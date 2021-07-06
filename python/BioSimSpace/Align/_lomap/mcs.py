# ******************
# MODULE DOCSTRING
# ******************

"""

LOMAP: Maximum Common Subgraph and scoring calculations
=====

Alchemical free energy calculations hold increasing promise as an aid to drug 
discovery efforts. However, applications of these techniques in discovery 
projects have been relatively few, partly because of the difficulty of planning 
and setting up calculations. The Lead Optimization Mapper (LOMAP) is an 
automated algorithm to plan efficient relative free energy calculations between 
potential ligands within a substantial of compounds.

"""

# *****************************************************************************
# Lomap2: A toolkit to plan alchemical relative binding affinity calculations
# Copyright 2015 - 2016  UC Irvine and the Authors
#
# Authors: Dr Gaetano Calabro' and Dr David Mobley
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, see http://www.gnu.org/licenses/
# *****************************************************************************


# ****************
# MODULE IMPORTS
# ****************


from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
from rdkit.Chem import Draw
from rdkit.Chem import rdmolops
from rdkit import DataStructs
from rdkit.Geometry.rdGeometry import Point3D
import sys
import math
from rdkit import RDLogger
import logging
import argparse

# *******************************
# Maximum Common Subgraph Class
# *******************************


__all__ = ['MCS']


class MCS(object):
    """

    This class is used to compute the Maximum Common Subgraph (MCS) between two
    RDkit molecule objects and to score their similarity by using defined rules 
    
    """

    def __init__(self, moli, molj, options=argparse.Namespace(time=20, verbose='info', max3d=1000, threed=False)):
        """
        Initialization function
    
        Parameters
        ----------

        moli : RDKit molecule object 
            the first molecule used to perform the MCS calculation
        molj : RDKit molecule object 
            the second molecule used to perform the MCS calculation
        options : argparse python object 
            the list of user options 
       
        """

        def substructure_centre(mol, mol_sub):
            """

            This function takes a molecule and a list of atom indices
            in that molecule and returns an RDKit Point3D representing
            the geometric centre of the atoms in the list

            """

            sum = Point3D()
            for i in mol_sub:
                sum += mol.GetConformer().GetAtomPosition(i)
            return sum / len(mol_sub)


        def best_substruct_match_to_mcs(moli,molj,by_rmsd=True):
            """

            This function looks over all of the substructure matches and returns the one
            with the best 3D correspondence (if by_rmsd is true), or the fewest number
            of atomic number mismatches (if by_rmsd is false)

            Note that the 3D correspondence does a translational centreing (but
            does not rotate).

            """

            # Sanity checking
            if not moli.HasSubstructMatch(self.mcs_mol):
                raise ValueError('RDkit MCS Subgraph first molecule search failed')

            if not molj.HasSubstructMatch(self.mcs_mol):
                raise ValueError('RDkit MCS Subgraph second molecule search failed')

            moli_sub = moli.GetSubstructMatches(self.mcs_mol,uniquify=False)
            molj_sub = molj.GetSubstructMatches(self.mcs_mol,uniquify=False)
            best_rmsd=1e8
            for mapi in moli_sub:
                for mapj in molj_sub:
                    # Compute the translation to bring molj's centre over moli
                    coord_delta = 0
                    if by_rmsd:
                        coord_delta = (substructure_centre(moli,mapi)
                                 - substructure_centre(molj,mapj))
                    rmsd=0
                    for pair in zip(mapi,mapj):
                        if by_rmsd:
                            rmsd += (moli.GetConformer().GetAtomPosition(pair[0]) 
                                   - molj.GetConformer().GetAtomPosition(pair[1])
                                   - coord_delta).LengthSq()
                        elif (moli.GetAtomWithIdx(pair[0]).GetAtomicNum() != 
                              molj.GetAtomWithIdx(pair[1]).GetAtomicNum()):
                            rmsd+=1
                    if rmsd < best_rmsd:
                        besti=mapi
                        bestj=mapj
                        best_rmsd=rmsd

            return (besti,bestj)

        def trim_mcs_mol(max_deviation=2.0):
            """

            This function is used to trim the MCS molecule to remove mismatched atoms i.e atoms
            where the topological mapping does not work in 3D coordinates.

            The sets of mapped atoms are translated to bring their geometric centres
            into alignment before trimming
           
            Parameters
            ----------

            max_deviation : the maximum difference in Angstroms between mapped atoms to allow

            """

            while True:
                (mapi,mapj) = best_substruct_match_to_mcs(self.__moli_noh,self.__molj_noh,by_rmsd=True)
                # Compute the translation to bring molj's centre over moli
                coord_delta = (substructure_centre(self.__moli_noh,mapi)
                             - substructure_centre(self.__molj_noh,mapj))
                worstatomidx=-1
                worstdist=0
                atomidx=0
                for pair in zip(mapi,mapj):
                    dist = (self.__moli_noh.GetConformer().GetAtomPosition(pair[0])
                          - self.__molj_noh.GetConformer().GetAtomPosition(pair[1])
                          - coord_delta).Length()
                    if dist > worstdist:
                        worstdist=dist
                        worstatomidx=atomidx
                    atomidx=atomidx+1

                if worstdist > max_deviation:
                    # Remove the furthest-away atom and try again
                    rwm = Chem.RWMol(self.mcs_mol)
                    rwm.RemoveAtom(worstatomidx)
                    if options.verbose == 'pedantic':
                       logging.info('Removing atom %d from MCS based on distance %f' %(worstatomidx,worstdist))
                    self.mcs_mol=Chem.Mol(rwm)
                else:
                    break

        def trim_mcs_fix_broken_rdkit_code():
            """

            Detect cases where the RDKit has generated an incorrect MCS (for alpha vs beta naphthyl, for
            instance), and delete some atoms to break the ring. The excess atoms will then be
            removed later by delete_broken_ring()
            
            Can be removed once RDKit is fixed

            Algorithm: find a bond in moli where the atoms are both in the MCS, they are bonded in
            moli, but are not bonded in the MCS.

            """

            to_remove = []
            for ai in self.moli.GetAtoms():
                if ai.HasProp('to_mcs'):    # is ai in the MCS?
                    aimcs = int(ai.GetProp('to_mcs'))
                    for bai in ai.GetNeighbors():
                        if bai.HasProp('to_mcs'):  # Atom bonded to ai is also in the MCS
                            baimcs = int(bai.GetProp('to_mcs'))
                            if (aimcs<baimcs):  # only do each bond once!
                                # Check if the corresponding MCS atoms are bonded
                                if not self.mcs_mol.GetBondBetweenAtoms(aimcs,baimcs):
                                    to_remove.append(aimcs)
                                    if options.verbose == 'pedantic':
                                       logging.info('Bond in first mol between atoms %d and %d not matched in MCS' %(ai.GetIdx(),bai.GetIdx()))
                                    

            if to_remove:
                # Delete atoms from the MCS, highest index first
                to_remove.sort(reverse=True)

                if options.verbose == 'pedantic':
                   logging.info('Removing %d atoms from MCS based on detection of broken RDKit ring bond matching' %(len(to_remove)))

                edit_mcs_mol = Chem.EditableMol(self.mcs_mol)
                for i in to_remove:
                    edit_mcs_mol.RemoveAtom(i)

                self.mcs_mol = edit_mcs_mol.GetMol()
                map_mcs_mol()   # Regenerate mappings

        def trim_mcs_chiral_atoms():
            """
                Remove all atoms in the MCS where there might be a chirality inversion i.e.
                (a) the corresponding atoms in the input molecules are both chiral, and
                (b) the parity of the atom mapping in the input molecules is reversed

                Calls map_mcs_mol as it uses the mappings generated there. 

            """

            def permutation_parity(perm):
                """ 
                    Returns the parity of the provided permutation tuple/array: even parity
                    is True, odd parity is False.
                """
                parity = True
                for i in range(len(perm)-1):
                    for j in range(i+1,len(perm)):
                        if (perm[i]<perm[j]): parity = not parity

                return parity

            def atom_mcs_chiral_parity(a):
                """
                    Take the neighbours of chiral atom a. Get the index of each of these atoms
                    in the MCS. Combine the parity of this list with the chirality flag for
                    a to determine the "MCS parity".
                """
                nbrs=[]
                for aj in a.GetNeighbors():
                    try:
                        nbrs.append(int(aj.GetProp('to_mcs')))
                    except Exception:
                        nbrs.append(1000)   # should not be more than one!

                if not permutation_parity(nbrs):
                    if a.GetChiralTag()==Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW: return Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW
                    if a.GetChiralTag()==Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW: return Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW
                return a.GetChiralTag()

                
                        
            def flag_inverted_atoms_in_mcs():
                """
                    Flag all atoms in the MCS where the chirality is inverted between
                    moli and molj. Atoms are flagged with CHI_TETRAHEDRAL_CW
                """

                # Generate atommappings as they are useful below
                map_mcs_mol()

                chiral_at_moli = [seq[0] for seq in Chem.FindMolChiralCenters(moli)]
                chiral_at_molj = [seq[0] for seq in Chem.FindMolChiralCenters(molj)]

                invertedatoms = []

                for i in chiral_at_moli:
                    # Is atom i in the MCS?
                    ai = moli.GetAtomWithIdx(i)
                    if (ai.HasProp('to_mcs')):
                        for j in chiral_at_molj:
                            # Is atom j in the MCS?
                            aj = molj.GetAtomWithIdx(j)
                            if (aj.HasProp('to_mcs')):
                                # Are they the same atom?
                                if (ai.GetProp('to_mcs') == aj.GetProp('to_mcs')):
                                    # OK, atoms are both chiral, and match the same MCS atom.
                                    # Take the list of neighbours for ai, and get their indices in 
                                    # the MCS. Use the parity of this index list together with the
                                    # chiral parity of ai to work out the "MCS parity". Do the same
                                    # for aj and check if the two are the same.
                                    # 
                                    # If not, flag with the CHI_TETRAHEDRAL_CW property.
                                    pi = atom_mcs_chiral_parity(ai)
                                    pj = atom_mcs_chiral_parity(aj)
                                    if (pi!=pj):
                                        invertedatoms.append(int(aj.GetProp('to_mcs')))

                for i in invertedatoms:
                    mcsat = self.mcs_mol.GetAtomWithIdx(i)
                    mcsat.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
                    if options.verbose == 'pedantic':
                       logging.info('Inverted chiral atom detected: %d' %(i))


            # Flag inverted atoms
            flag_inverted_atoms_in_mcs()

            # Trim inverted chiral Atoms. The algorithm is to delete the chiral centre,
            # fragment the molecule, and keep only the two largest fragments. Rinse and
            # repeat until no more flagged chiral centres remain.

            while True:
                mcs_chiral_set = set()
                atom_idx = -1;

                for atom in self.mcs_mol.GetAtoms():
                    # Note that any atom in the MCS which has inverted chirality between the input mols is
                    # flagged with CHI_TETRAHEDRAL_CW
                    if (atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW):
                        atom_idx=atom.GetIdx()
                        atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_UNSPECIFIED)  # Remove from consideration for the next loop
                        break

                if atom_idx == -1:  # Not found any more chiral atoms, so done
                    break

                # Move the chiral atom to the end (avoids indexing problems)
                newindexes = list(range(self.mcs_mol.GetNumAtoms()))
                newindexes.remove(atom_idx)
                newindexes.append(atom_idx)
                self.mcs_mol = Chem.RenumberAtoms(self.mcs_mol,newindexes)

                # Now we loop, deleting groups attached to the chiral atom, until the 
                # chiral atom has at most two heavy atom connections
                # Note that getAtoms()[-1] returns the first atom not the last if you 
                # don't convert it to a list. Grr.
                while list(self.mcs_mol.GetAtoms())[-1].GetDegree() > 2 :

                    # Delete the chiral atom in a temporary molecule, and fragment. Since the
                    # chiral atom was the last one, the indexes in the temporary molecule are the
                    # same as in self.mcs_mol
                    edit_mol = Chem.EditableMol(self.mcs_mol)
                    edit_mol.RemoveAtom(self.mcs_mol.GetNumAtoms()-1)
                    tmp_mol = edit_mol.GetMol()
                    fragments = Chem.rdmolops.GetMolFrags(tmp_mol)

                    # Get index of smallest fragments
                    min_idx = 0
                    lgt_min = 10000

                    for idx in range(0, len(fragments)):
                        lgt = len(fragments[idx])
                        if lgt < lgt_min:
                            lgt_min = lgt
                            min_idx = idx

                    # Get the atoms in this fragment and sort them so we delete the
                    # largest index first
                    min_frag = list(fragments[min_idx])
                    min_frag.sort(reverse=True)

                    if options.verbose == 'pedantic':
                       logging.info('Removing %d atoms to remove chiral inversion' %(len(min_frag)))
                    edit_mol = Chem.EditableMol(self.mcs_mol)
                    for idx in min_frag:
                        edit_mol.RemoveAtom(idx)
                    self.mcs_mol = edit_mol.GetMol()

            map_mcs_mol()   # Regenerate mappings after deletion
            # Done!

        def delete_broken_ring():
            """
                This function checks the MCS to see if there are any
                atoms which are in a ring in the parent molecules, but
                not in a ring in the MCS. This may occur if we have deleted
                some atoms from the MCS in 3D coordinate matching, for
                example.
            """

            to_remove = []
            for at in self.mcs_mol.GetAtoms():
                moli_idx = int(at.GetProp('to_moli'))
                moli_at = self.__moli_noh.GetAtomWithIdx(moli_idx)
                molj_idx = int(at.GetProp('to_molj'))
                molj_at = self.__moli_noh.GetAtomWithIdx(moli_idx)

                # Testing moli and molj is redundant due to the way that the
                # MCS is calculated, but I'd rather be paranoid here
                if (moli_at.IsInRing() and molj_at.IsInRing() and not at.IsInRing()):
                    to_remove.append(at.GetIdx())

            if to_remove:
                # Delete atoms from the MCS, highest index first
                to_remove.sort(reverse=True)

                if options.verbose == 'pedantic':
                   logging.info('Removing %d atoms from MCS to clear up partial rings' %(len(to_remove)))

                edit_mcs_mol = Chem.EditableMol(self.mcs_mol)
                for i in to_remove:
                    edit_mcs_mol.RemoveAtom(i)

                self.mcs_mol = edit_mcs_mol.GetMol()

                map_mcs_mol()   # Regenerate mappings after deletion


        def map_mcs_mol():
            """

            This function is used to define a map between the generated mcs, the
            molecules and vice versa
           
            """

            # Get self-mapping for the MCS
            mcsi_sub = tuple(range(self.mcs_mol.GetNumAtoms()))

            (moli_sub,molj_sub) = best_substruct_match_to_mcs(self.__moli_noh,self.__molj_noh,by_rmsd=self.options.threed)

            # mcs to moli
            map_mcs_mol_to_moli_sub = list(zip(mcsi_sub, moli_sub))

            # Clear all properties as we may call this function more than once
            for a in self.mcs_mol.GetAtoms():
                a.ClearProp('to_moli')
                a.ClearProp('to_molj')
            for a in self.moli.GetAtoms():
                a.ClearProp('to_mcs')
            for a in self.molj.GetAtoms():
                a.ClearProp('to_mcs')

            # An RDkit atomic property is defined to store the mapping to moli
            for idx in map_mcs_mol_to_moli_sub:
                self.mcs_mol.GetAtomWithIdx(idx[0]).SetProp('to_moli', str(idx[1]))
                self.moli.GetAtomWithIdx(idx[1]).SetProp('to_mcs', str(idx[0]))

            mcsj_sub = tuple(range(self.mcs_mol.GetNumAtoms()))

            # mcs to molj
            map_mcs_mol_to_molj_sub = list(zip(mcsj_sub, molj_sub))

            # Map between the two molecules
            self.__map_moli_molj = list(zip(moli_sub, molj_sub))

            # An RDkit atomic property is defined to store the mapping to molj
            for idx in map_mcs_mol_to_molj_sub:
                self.mcs_mol.GetAtomWithIdx(idx[0]).SetProp('to_molj', str(idx[1]))
                self.molj.GetAtomWithIdx(idx[1]).SetProp('to_mcs', str(idx[0]))

            # For each mcs atom we save its original index in a specified 
            # property. This could be very useful in the code development
            # when deletion or atom insertions are performed
            for at in self.mcs_mol.GetAtoms():
                at.SetProp('org_idx', str(at.GetIdx()))

            return

        def set_ring_counter(mol):

            """

            This function is used to attach to each molecule atom a ring counter
            rc. This parameter is used to asses if a ring has been broken or not
            during the MCS mapping
         
            Parameters
            ----------
            mol : RDKit Molecule obj
                the molecule used to define the atom ring counters
            """

            # set to zero the atom ring counters
            for at in mol.GetAtoms():
                at.SetProp('rc', '0')

            rginfo = mol.GetRingInfo()

            rgs = rginfo.AtomRings()

            # print rgs

            rgs_set = set([e for l in rgs for e in l])

            for idx in rgs_set:
                for r in rgs:
                    if idx in r:
                        val = int(mol.GetAtomWithIdx(idx).GetProp('rc'))
                        val = val + 1
                        mol.GetAtomWithIdx(idx).SetProp('rc', str(val))
            return

        # START of __init__ function
        # Set logging level and format
        logging.basicConfig(format='%(levelname)s:\t%(message)s', level=logging.INFO)

        self.options=options

        # Global beta setting for atom penalties
        self.beta = 0.1

        # Local pointers to the passed molecules
        self.moli = moli
        self.molj = molj

        # Sanitize input molecules
        Chem.SanitizeMol(self.moli)
        Chem.SanitizeMol(self.molj)

        # Set chirality flags from 3D coords if working in 3D
        if self.options.threed:
            Chem.rdmolops.AssignAtomChiralTagsFromStructure(self.moli,replaceExistingTags=True)
            Chem.rdmolops.AssignAtomChiralTagsFromStructure(self.molj,replaceExistingTags=True)

        if not options.verbose == 'pedantic':
            lg = RDLogger.logger()
            lg.setLevel(RDLogger.CRITICAL)

        # Local pointers to the passed molecules without hydrogens
        # These variables are defined as private
        try:
            self.__moli_noh = AllChem.RemoveHs(moli)
            self.__molj_noh = AllChem.RemoveHs(molj)
        except Exception:
            self.__moli_noh = AllChem.RemoveHs(moli, sanitize=False)
            self.__molj_noh = AllChem.RemoveHs(molj, sanitize=False)

            Chem.SanitizeMol(self.__moli_noh, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
            Chem.SanitizeMol(self.__molj_noh, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)

        # MCS calculation. In RDKit the MCS is a smart string. Ring atoms are 
        # always mapped in ring atoms. 
        # Don't add the mcs result as a member variable as it can't be pickled
        __mcs = rdFMCS.FindMCS([self.__moli_noh, self.__molj_noh],
                                    timeout=options.time,
                                    atomCompare=rdFMCS.AtomCompare.CompareAny,
                                    bondCompare=rdFMCS.BondCompare.CompareAny,
                                    matchValences=False,
                                    ringMatchesRingOnly=True,
                                    completeRingsOnly=True,
                                    matchChiralTag=False)

        # Note that we need matchChiralTag=False as we want to match chiral atoms with different
        # parities, we just need to trim the MCS to the largest possible match that doesn't have
        # a mismatched chiral centre in it (eg we want to match FC[C@](C)CO to FC[C@@](C)CO
        # using the MCS FCCCO - this includes the chiral atom but we delete the methyl group
        # that led it to be chiral. The trimming is done in trim_mcs_chiral_atoms()

        # Checking
        if __mcs.canceled:
            logging.warning('Timeout reached to find the MCS between the molecules')

        if __mcs.numAtoms == 0:
            raise ValueError('No MCS was found between the molecules')

        # The found MCS pattern (smart strings) is converted to a RDKit molecule
        self.mcs_mol_smarts = __mcs.smartsString
        self.mcs_mol = Chem.MolFromSmarts(__mcs.smartsString)

        # There's a symmetry-related bug here: if there was more than one MCS
        # of the same size and score, we'll get only one at random. We then try
        # to choose the mapping that matches 3D coords the best, but one of the 
        # not-considered MCSes that we never saw may give a better mapping. 
        # We can rescue some of this by converting all partial-query atoms to 
        # full query atoms
        testmol = Chem.MolFromSmarts("*")   # Create a "match anything" query atom for us to copy
        for a in self.mcs_mol.GetAtoms():
            if a.DescribeQuery().startswith("AtomOr"):  # Matches more than one element
                a.SetQuery(testmol.GetAtoms()[0])   # Set this atom to a copy of the "match anything" atom

        try:  # Try to sanitize the MCS molecule
            Chem.SanitizeMol(self.mcs_mol)
        except Exception:  # if not, try to recover the atom aromaticity which is
            # important for the ring counter
            sanitFail = Chem.SanitizeMol(self.mcs_mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY,
                                         catchErrors=True)
            if sanitFail:  # if not, the MCS is skipped
                raise ValueError('Sanitization Failed...')

        # Trim the MCS to remove atoms with too-large real-space deviations
        if self.options.max3d>0 :
            try:
                trim_mcs_mol(max_deviation=self.options.max3d)
            except Exception as e:
                raise ValueError(str(e))

        # Trim the MCS further to remove chirality mismatches
        trim_mcs_chiral_atoms()

        # Check to see if we've hit the RDKit incorrect-MCS bug
        trim_mcs_fix_broken_rdkit_code()

        # Cleanup any partial rings remaining
        delete_broken_ring()

        # Mapping between the found MCS molecule and moli,  molj
        try:
            map_mcs_mol()
        except Exception as e:
            raise ValueError(str(e))

        # Set the ring counters for each molecule
        set_ring_counter(self.__moli_noh)
        set_ring_counter(self.__molj_noh)
        set_ring_counter(self.mcs_mol)

        # for at in self.mcs_mol.GetAtoms():
        #     print 'at = %d rc = %d' % (at.GetIdx(), int(at.GetProp('rc')))

        if not options.verbose == 'pedantic':
            lg.setLevel(RDLogger.WARNING)

        return

    def get_map(self):
        """

        This function is used to return a list of pairs of atom indexes generated
        by the mapping between the two molecules used to calculate the MCS. 
        The calculated mapping is performed without considering hydrogens 

        Returns
        -------
        pair of indexes related to the atom mapping 

        """

        return self.__map_moli_molj

    ############ MCS BASED RULES ############

    def mcsr(self):

        """
        This rule computes the similarity between the two passed molecules 
        used to compute the MCS
        
        Returns
        -------
        scr_mcsr : float
            the rule score

             
        """

        # The number of heavy atoms in each molecule
        nha_moli = self.moli.GetNumHeavyAtoms()
        nha_molj = self.molj.GetNumHeavyAtoms()
        # Note that the mcs_mol (a) doesn't contain hydrogens, and (b) does contain
        # wildcard atoms, which don't count as 'heavy'. Use the total atom count instead.
        nha_mcs_mol = self.mcs_mol.GetNumAtoms()

        # score
        scr_mcsr = math.exp(-self.beta * (nha_moli + nha_molj - 2 * nha_mcs_mol))

        logging.info('MCSR from MCS size %d, molecule sizes %d,%d is %f' %(nha_mcs_mol,nha_moli,nha_molj,scr_mcsr))

        return scr_mcsr

    # MNACR rule
    def mncar(self, ths=4):

        """
        This rule cut the similarity score between two molecules if they do
        not share the selected number of atoms 

        
        Parameters
        ----------
        ths : float
            the minumum number of atoms to share
        
        Returns
        -------
        scr_mncar : float
            the rule score     
        """

        # This rule has been modified from the rule desribed in the Lomap paper
        # to match the LOMAP first implementation provided by schrodinger

        nha_mcs_mol = self.mcs_mol.GetNumHeavyAtoms()
        nha_moli = self.moli.GetNumHeavyAtoms()
        nha_molj = self.molj.GetNumHeavyAtoms()

        scr_mncar = float((nha_mcs_mol >= ths) or (nha_moli < ths + 3) or (nha_molj < ths + 3))

        return scr_mncar

    # TMCRS rule (Trim rule) 
    # MDM Note: we don't use this as we don't have the same limitation on partial ring
    # deletion as Schrodinger
    # NB removed the chirality check - the MCS is now trimmed to remive chirality
    def tmcsr(self, strict_flag=True):
        return 1.0

    # AtomicNumber rule 
    def atomic_number_rule(self):

        """
        This rule checks how many elements have been changed in the MCS 
        and a score based on the fraction of MCS matches that are the same atomic number.
        When used with beta=0.1 and multiplied by mcsr, this is equivalent to counting
        mismatched atoms at only half weight.

        This has been extended to modify the amount of mismatch according to the 
        atoms being mapped. 
             
        """

        # A value of 0.5 is the same behaviour as before, a value of 1 means that the 
        # atoms are perfectly equivalent, a value of 0 means that the atoms are perfectly
        # non-equivalent (i.e the penalty should basically remove this atom pair from the
        # MCS). The default for pairs not in this data structure is 0.5. 
        # 
        # Note that we don't need the symmetry equivalent values: we will use the large of 
        # [i][j] and [j][i]
        transform_difficulty={ 
          # H to element - not sure this has any effect currently
          1: { 9: 0.5, 17: 0.25, 35: 0, 53: -0.5 },
          # O to element - methoxy to Cl/Br is easier than expected
          8: { 17: 0.85, 35: 0.85 },
          # F to element 
          9: { 17: 0.5, 35: 0.25, 53: 0 },
          # Cl to element 
          17: { 35: 0.85, 53: 0.65 },
          # Br to element
          35: { 53: 0.85 },
        }
        nmismatch=0
        for at in self.mcs_mol.GetAtoms():
            moli_idx = int(at.GetProp('to_moli'))
            molj_idx = int(at.GetProp('to_molj'))
            moli_a = self.__moli_noh.GetAtoms()[moli_idx]
            molj_a = self.__molj_noh.GetAtoms()[molj_idx]

            if moli_a.GetAtomicNum() != molj_a.GetAtomicNum():
                ij=-1
                ji=-1
                try:
                    ij=transform_difficulty[moli_a.GetAtomicNum()][molj_a.GetAtomicNum()]
                except KeyError:
                    pass
                try:
                    ji=transform_difficulty[molj_a.GetAtomicNum()][moli_a.GetAtomicNum()]
                except KeyError:
                    pass
                diff = max(ij,ji)
                if (diff==-1):
                    diff=0.5    # default for elements not found
                    
                nmismatch+=(1-diff)

        an_score =  math.exp(-1 * self.beta * nmismatch)
        logging.info('atomic number score from %d mismatches is %f' %(nmismatch,an_score))
        return an_score

    # Hybridization rule 
    def hybridization_rule(self, penalty_weight = 1.5):

        """
        This rule checks how many atoms have changed hybridization state.
        The penalty weight means how many "atoms" different a hybridization state change
        is: 1 means that the atom is effectively removed from the MCS for scoring purposes,
        0 means that hybridization changes are free.
        When used with beta=0.1 and multiplied by mcsr, this is equivalent to counting
        mismatched atoms at a weight of (1-penalty_weight)

        """

        def atom_hybridization(a):
            """
            RDKit has an un-useful hybridization definition. Instead, just look at the number
            of multiple bonds from an atom
            """
            if a.GetIsAromatic(): return 2
            xs=0
            for b in a.GetBonds():
                if b.GetBondType()==Chem.rdchem.BondType.AROMATIC: return 2
                if b.GetBondType()==Chem.rdchem.BondType.DOUBLE: xs += 1
                if b.GetBondType()==Chem.rdchem.BondType.TRIPLE: xs += 2
                if b.GetBondType()==Chem.rdchem.BondType.ONEANDAHALF: xs += 0.5

            # O- is sp2 to avoid problems with carboxylate etc
            if a.GetAtomicNum()==8 and a.GetFormalCharge()<0: return 2 # sp2

            if xs==0: return 3 # sp3
            if xs>1.1: return 1 # sp
            return 2  # sp2

        nmismatch=0
        for at in self.mcs_mol.GetAtoms():
            moli_idx = int(at.GetProp('to_moli'))
            molj_idx = int(at.GetProp('to_molj'))
            moli_a = self.__moli_noh.GetAtoms()[moli_idx]
            molj_a = self.__molj_noh.GetAtoms()[molj_idx]

            hybi = atom_hybridization(moli_a) 
            hybj = atom_hybridization(molj_a) 
            mismatch= hybi != hybj

            # Allow Nsp3 to match Nsp2, otherwise guanidines etc become painful
            if moli_a.GetAtomicNum()==7 and molj_a.GetAtomicNum()==7 and (hybi in [2,3]) and hybj in [2,3]: mismatch=False

            if mismatch:
                nmismatch+=1
                logging.info("Hybridization mismatch %d %s %d vs %d %s %d",moli_a.GetIdx(),moli_a.GetSymbol(),hybi,molj_a.GetIdx(),molj_a.GetSymbol(),hybj)

        hyb_score =  math.exp(-1 * self.beta * nmismatch * penalty_weight)
        logging.info('hybridization score from %d mismatches is %f' %(nmismatch,hyb_score))
        return hyb_score


    # Sulfonamides rule
    def sulfonamides_rule(self, penalty=4):

        """
        This rule checks to see if we are growing a complete sulfonamide, and 
        returns 0 if we are. This means that if this rule is used we effectively disallow
        this transition. Testing has shown that growing -SO2NH2 from scratch performs
        very badly.

        Parameters
        ----------
        penalty : the number of atom mismatches that failing this rule will lower the score by
             
        """

        def adds_sulfonamide(mol):
            """
            Returns true if the removal of the MCS from the provided molecule
            leaves a sulfonamide
            """

            if not mol.HasSubstructMatch(self.mcs_mol):
                raise ValueError('RDkit MCS Subgraph molecule search failed in sulfonamide check')

            
            rwm=rdmolops.DeleteSubstructs(mol, self.mcs_mol)
            return rwm.HasSubstructMatch(Chem.MolFromSmarts('S(=O)(=O)N'))

        fail = 1 if (adds_sulfonamide(self.__moli_noh)) else 0
        fail = 1 if (adds_sulfonamide(self.__molj_noh)) else fail
        sulf_score =  math.exp(-1 * self.beta * fail * penalty)
        logging.info('sulfonamide score is %f' %(sulf_score))
        return sulf_score

    # Heterocycles rule
    def heterocycles_rule(self, penalty=4):

        """
        This rule checks to see if we are growing a heterocycle from a hydrogen, and 
        returns <1 if we are. This means that if this rule is used we penalise
        this transition. Testing has shown that growing a pyridine or other heterocycle
        is unlikely to work (better to grow phenyl then mutate)

        Parameters
        ----------
        penalty : the number of atom mismatches that failing this rule will lower the score by
             
             
        """

        def adds_heterocycle(mol):
            """
            Returns true if the removal of the MCS from the provided molecule
            leaves a sulfonamide
            """

            if not mol.HasSubstructMatch(self.mcs_mol):
                raise ValueError('RDkit MCS Subgraph molecule search failed in sulfonamide check')

            
            rwm=rdmolops.DeleteSubstructs(mol, self.mcs_mol)
            # Only picking up N/C containing heterocycles - odd cases like pyran derivatives are not caught
            grow6mheterocycle =  rwm.HasSubstructMatch(Chem.MolFromSmarts('[n]1[c,n][c,n][c,n][c,n][c,n]1'))

            # Note that growing pyrrole, furan or thiophene is allowed
            grow5mheterocycle =  rwm.HasSubstructMatch(Chem.MolFromSmarts('[o,n,s]1[n][c,n][c,n][c,n]1'))
            grow5mheterocycle |=  rwm.HasSubstructMatch(Chem.MolFromSmarts('[o,n,s]1[c,n][n][c,n][c,n]1'))
            return (grow6mheterocycle | grow5mheterocycle)


        fail = 1 if (adds_heterocycle(self.__moli_noh)) else 0
        fail = 1 if (adds_heterocycle(self.__molj_noh)) else fail
        het_score = math.exp(-1 * self.beta * fail * penalty)
        logging.info('heterocycle score is %f' %(het_score))
        return het_score

    def transmuting_methyl_into_ring_rule(self, penalty=6):

        """
         Rule to prevent turning a methyl into a ring atom and similar transformations
         (you can grow a ring, but you can't transmute into one)

        Parameters
        ----------
        penalty : the number of atom mismatches that failing this rule will lower the score by
             

        """
        moli=self.__moli_noh
        molj=self.__molj_noh

        # Get list of bonds in mol i and j that go from the MCS to a non-MCS atom,
        # arranged in tuples with the index of the MCS atom
        moli_sub = moli.GetSubstructMatch(self.mcs_mol)
        molj_sub = molj.GetSubstructMatch(self.mcs_mol)

        is_bad=False

        for i in range(0,len(moli_sub)):
            edge_bondsi = [ b.GetBeginAtomIdx() for b in moli.GetBonds() if (b.GetEndAtomIdx()==moli_sub[i] and not b.GetBeginAtomIdx() in moli_sub) ]
            edge_bondsi += [ b.GetEndAtomIdx() for b in moli.GetBonds() if (b.GetBeginAtomIdx()==moli_sub[i] and not b.GetEndAtomIdx() in moli_sub) ]
            edge_bondsj = [ b.GetBeginAtomIdx() for b in molj.GetBonds() if (b.GetEndAtomIdx()==molj_sub[i] and not b.GetBeginAtomIdx() in molj_sub) ]
            edge_bondsj += [ b.GetEndAtomIdx() for b in molj.GetBonds() if (b.GetBeginAtomIdx()==molj_sub[i] and not b.GetEndAtomIdx() in molj_sub) ]
            #print("Atom",i,"index",moli_sub[i],"edge atoms on mol 1 are",edge_bondsi);
            #print("Atom",i,"index",molj_sub[i],"edge atoms on mol 2 are",edge_bondsj);

            for edgeAtom_i in edge_bondsi:
                for edgeAtom_j in edge_bondsj:
                    if (moli.GetAtomWithIdx(edgeAtom_i).IsInRing() ^ molj.GetAtomWithIdx(edgeAtom_j).IsInRing()):
                        is_bad=True

        mescore = math.exp(-1 * self.beta * penalty) if is_bad else 1
        logging.info('methyl-to-ring transformation score is %f' %(mescore))
        return mescore

    def transmuting_ring_sizes_rule(self):

        """
         Rule to prevent turning a ring atom into a ring atom with a different ring size
         (you can grow a ring, but you can't turn a cyclopentyl into a cyclohexyl)

         Hard rule: sets score to near zero if violated

        """
        moli=self.__moli_noh
        molj=self.__molj_noh

        # Get list of bonds in mol i and j that go from the MCS to a non-MCS atom,
        # arranged in tuples with the index of the MCS atom
        moli_sub = moli.GetSubstructMatch(self.mcs_mol)
        molj_sub = molj.GetSubstructMatch(self.mcs_mol)

        is_bad=False

        for i in range(0,len(moli_sub)):
            edge_bondsi = [ b.GetBeginAtomIdx() for b in moli.GetBonds() if (b.GetEndAtomIdx()==moli_sub[i] and not b.GetBeginAtomIdx() in moli_sub) ]
            edge_bondsi += [ b.GetEndAtomIdx() for b in moli.GetBonds() if (b.GetBeginAtomIdx()==moli_sub[i] and not b.GetEndAtomIdx() in moli_sub) ]
            edge_bondsj = [ b.GetBeginAtomIdx() for b in molj.GetBonds() if (b.GetEndAtomIdx()==molj_sub[i] and not b.GetBeginAtomIdx() in molj_sub) ]
            edge_bondsj += [ b.GetEndAtomIdx() for b in molj.GetBonds() if (b.GetBeginAtomIdx()==molj_sub[i] and not b.GetEndAtomIdx() in molj_sub) ]
            #print("Atom",i,"index",moli_sub[i],"edge atoms on mol 1 are",edge_bondsi);
            #print("Atom",i,"index",molj_sub[i],"edge atoms on mol 2 are",edge_bondsj);

            for edgeAtom_i in edge_bondsi:
                for edgeAtom_j in edge_bondsj:
                    #print("Checking ring for atom",edgeAtom_i,edgeAtom_j,moli.GetAtomWithIdx(edgeAtom_i).IsInRing(),molj.GetAtomWithIdx(edgeAtom_j).IsInRing())
                    if (moli.GetAtomWithIdx(edgeAtom_i).IsInRing() and molj.GetAtomWithIdx(edgeAtom_j).IsInRing()):
                        for ring_size in range(3,8):
                            if (moli.GetAtomWithIdx(edgeAtom_i).IsInRingSize(ring_size) ^ molj.GetAtomWithIdx(edgeAtom_j).IsInRingSize(ring_size)):
                                logging.info('tRansforming ring sizes score is 0 based on atom %d in moli and %d in molj' %(edgeAtom_i,edgeAtom_j))
                                is_bad=True
                            if (moli.GetAtomWithIdx(edgeAtom_i).IsInRingSize(ring_size) or molj.GetAtomWithIdx(edgeAtom_j).IsInRingSize(ring_size)):
                                break

        return 0.1 if is_bad else 1

    def heavy_atom_mcs_map(self):
        '''
        Returns a list of tuples mapping atoms from moli to molj
        Heavy atoms only, returned sorted by first index
        '''
        maplist=[]
        for at in self.mcs_mol.GetAtoms():
            moli_idx = int(at.GetProp('to_moli'))
            molj_idx = int(at.GetProp('to_molj'))
            maplist.append((moli_idx,molj_idx))
        maplist.sort()
        return maplist

    def heavy_atom_match_list(self):
        '''
        Returns a string listing the MCS match between the two molecules as 
          atom_m1:atom_m2,atom_m1:atom_m2,...
        Heavy atoms only
        '''
        return ",".join([str(i)+":"+str(j) for (i,j) in self.heavy_atom_mcs_map()])

    def all_atom_match_list(self):
        '''
        Returns a string listing the MCS match between the two molecules as 
          atom_m1:atom_m2,atom_m1:atom_m2,...
        All atoms including hydrogens. The string is sorted by first index.
        We need to be careful that this function is symmetric, and that hydrogens
        are mapped correctly.
        '''

        def get_attached_atoms_not_in_mcs(mol,i):
            ''' Get atoms attached to atom i which are not in the MCS '''
            attached=[]
            for b in mol.GetBonds():
                if b.GetEndAtomIdx()==i or b.GetBeginAtomIdx()==i:
                    j=b.GetEndAtomIdx()
                    if (j==i):
                        j=b.GetBeginAtomIdx()
                    # OK, so j is the atom at the other end of the bond atom atom i. Is it in the MCS?
                    inMCS = mol.GetAtomWithIdx(j).HasProp('to_mcs')
                    if not inMCS:
                        attached.append(j)
            return attached


        moli=self.moli
        molj=self.molj

        maplist=self.heavy_atom_mcs_map()

        # OK, this is painful, as the MCS only includes heavies. We now need to match up
        # hydrogens hanging off the MCS

        # Iterate over all atoms in the MCS
        for at in self.mcs_mol.GetAtoms():
            moli_idx = int(at.GetProp('to_moli'))
            molj_idx = int(at.GetProp('to_molj'))
            attached_i = get_attached_atoms_not_in_mcs(moli,moli_idx)
            attached_j = get_attached_atoms_not_in_mcs(molj,molj_idx)

            # Now, we need to match these up, with the caveat that we *must* not match
            # a heavy to a heavy (as if we were allowed to match these, then they would be
            # in the MCS! 
            #
            # In 3D mode, ensure that maps happen tothe closest atom in 3D coordinates - 
            # this gets mappings to prochiral hydrogens correct (SFT-15791)

            # Match H to H first
            while attached_i and attached_j:
                
                hidx_i=-1
                hidx_j=-1
                best_dist=10000
                for ai in attached_i:
                    if moli.GetAtomWithIdx(ai).GetAtomicNum()==1:
                        for aj in attached_j:
                            if molj.GetAtomWithIdx(aj).GetAtomicNum()==1:
                                dist = (moli.GetConformer().GetAtomPosition(ai)
                                      - molj.GetConformer().GetAtomPosition(aj)).Length()
                                if (dist < best_dist or not self.options.threed):
                                    hidx_i=ai
                                    hidx_j=aj
                                    best_dist=dist
                if (hidx_i<0):
                    # OK, no hydrogen-hydrogen matches left. Try to match a hydrogen to a non-hydrogen
                    for ai in attached_i:
                        for aj in attached_j:
                            if moli.GetAtomWithIdx(ai).GetAtomicNum()==1 or molj.GetAtomWithIdx(aj).GetAtomicNum()==1:
                                dist = (moli.GetConformer().GetAtomPosition(ai)
                                      - molj.GetConformer().GetAtomPosition(aj)).Length()
                                if (dist < best_dist or not self.options.threed):
                                    hidx_i=ai
                                    hidx_j=aj
                                    best_dist=dist

                if (hidx_i>=0):
                    # Found a mappable pair: add and try again
                    maplist.append((hidx_i,hidx_j))
                    attached_i.remove(hidx_i)
                    attached_j.remove(hidx_j)
                else:
                    break   # No mappable pairs left

        maplist.sort()
        return ",".join([str(i)+":"+str(j) for (i,j) in maplist])

"""
Table of #atoms-changed to score for beta=0.1

0   1.0
1   0.905
2   0.819
3   0.741
4   0.670
5   0.607
6   0.549
7   0.500
8   0.449
9   0.407
10  0.369
11  0.333
12  0.301
13  0.273
14  0.247
15  0.223
16  0.202
"""
              
if "__main__" == __name__:

    # MCS is wrong unless you convert down to SMILES and back up
    #mola = Chem.MolFromMolFile('../p38_3flz.sdf', sanitize=False, removeHs=False)
    #molb = Chem.MolFromMolFile('../p38_2i.sdf', sanitize=False, removeHs=False)

    # Chirality problem with 4- vs 3-heavy atom attachment chiral centres
    #mola = Chem.MolFromMolFile('../bace_mk1.sdf', sanitize=False, removeHs=False)
    #molb = Chem.MolFromMolFile('../bace_cat_13d.sdf', sanitize=False, removeHs=False)

    # Me -> OMe wrong as O is labeled SP2
    #mola = Chem.MolFromMolFile('../jnk_18630.sdf', sanitize=False, removeHs=False)
    #molb = Chem.MolFromMolFile('../jnk_18633.sdf', sanitize=False, removeHs=False)
    # Cl -> OMe wrong as O is labelled SP2
    #mola = Chem.MolFromMolFile('../jnk_18626.sdf', sanitize=False, removeHs=False)
    #molb = Chem.MolFromMolFile('../jnk_18632.sdf', sanitize=False, removeHs=False)

    # These two have an alpha/beta ring attachment change but get an MCS that includes all atoms
    #mola = Chem.MolFromMolFile('../mcl_46.sdf', sanitize=False, removeHs=False)
    #molb = Chem.MolFromMolFile('../mcl_48.sdf', sanitize=False, removeHs=False)
    
    # Hybridization issues with amide N being labeled SP2
    #mola = Chem.MolFromMolFile('../tyk_42.sdf', sanitize=False, removeHs=False)
    #molb = Chem.MolFromMolFile('../tyk_54.sdf', sanitize=False, removeHs=False)

    # Chirality testing
    #mola = Chem.MolFromMolFile('../test/chiral/bace_cat_13d.sdf', sanitize=False, removeHs=False)
    #molb = Chem.MolFromMolFile('../test/chiral/bace_cat_13d_perm2.sdf', sanitize=False, removeHs=False)

    # Testing bridge addition: need to map heavy atoms to the correct prochiral hydrogen
    mola = Chem.MolFromMolFile('../test/chiral/tpbs2_lig1.sdf', sanitize=False, removeHs=False)
    molb = Chem.MolFromMolFile('../test/chiral/tpbs2_lig2.sdf', sanitize=False, removeHs=False)

    #mola = Chem.MolFromSmiles(Chem.MolToSmiles(mola))
    #molb = Chem.MolFromSmiles(Chem.MolToSmiles(molb))
    print("Mola: ",Chem.MolToSmiles(mola))
    print("Molb: ",Chem.MolToSmiles(molb))

    # MCS calculation
    try:
        MC = MCS(mola, molb, argparse.Namespace(time=20, verbose='pedantic', max3d=5, threed=True))
        #MC = MCS(mola, molb, argparse.Namespace(time=20, verbose='info', max3d=0, threed=False))
        #MC = MCS(mola, molb)
    except Exception:
        raise ValueError('NO MCS FOUND......')

    # # Rules calculations
    mcsr = MC.mcsr()
    strict = MC.tmcsr(strict_flag=True)
    loose = MC.tmcsr(strict_flag=False)
    mncar = MC.mncar()
    atnum = MC.atomic_number_rule()
    hybrid = MC.hybridization_rule()
    sulf = MC.sulfonamides_rule()
    het = MC.heterocycles_rule()
    growring = MC.transmuting_methyl_into_ring_rule()
    changering = MC.transmuting_ring_sizes_rule()

    print('TMCRS STRICT = %f , TMCRS LOOSE = %f' % (strict, loose))
    print('MCSR = ', mcsr)
    print('MNCAR = ', mncar)
    print('ATNUM = ', atnum)

    tmp = mcsr * mncar

    print('Total Strict = %f , Total Loose = %f' % (tmp * strict, tmp * loose))

    print('MCS is ',MC.mcs_mol.GetNumHeavyAtoms(),' ',Chem.MolToSmiles(MC.mcs_mol))
    for at in MC.mcs_mol.GetAtoms():
        moli_idx = int(at.GetProp('to_moli'))
        molj_idx = int(at.GetProp('to_molj'))
        moli_a = mola.GetAtoms()[moli_idx]
        molj_a = molb.GetAtoms()[molj_idx]
        print("MCS match: ",moli_idx,moli_a.GetAtomicNum(),molj_idx,molj_a.GetAtomicNum())

    print("hybridization change:",hybrid)
    print("sulfonamides:",sulf)
    print("heterocycles:",het)
    print("growring:",growring)
    print("transmuting_ring_sizes_rule:",changering)
    score = mncar * mcsr * atnum * hybrid
    score *= sulf * het * growring
    score *= changering
    print("FINAL SCORE:",score)
    print("Match list:",MC.heavy_atom_match_list())
    print("Match list:",MC.all_atom_match_list())

