#!/usr/bin/env python
# coding: utf-8

# Author: Julien Michel
# 
# email: julien.michel@ed.ac.uk

# # PrepareFEP
# Loads a pair of input files, perform mapping between the first molecule of each input. Write down input files for a SOMD FEP calculation.

import os
import zipfile
from Sire.Mol import AtomIdx
import BioSimSpace as BSS

import os
import re
import sys
import argparse
from Sire.Base import *
from datetime import datetime
# Make sure that the OPENMM_PLUGIN_DIR enviroment variable is set correctly.
os.environ["OPENMM_PLUGIN_DIR"] = getLibDir() + "/plugins"

from Sire.IO import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.System import *
from Sire.Move import *
from Sire.MM import *
from Sire.FF import *
from Sire.Units import *
from Sire.Vol import *
from Sire.Maths import *
from Sire.Qt import *
from Sire.ID import *
from Sire.Config import *
from Sire.Analysis import *
from Sire.Tools.DCDFile import *
from Sire.Tools import Parameter, resolveParameters
import Sire.Stream
import time
import numpy as np


node = BSS.Gateway.Node("A node to generate input files for a SOMD relative free energy calculation.")
node.addInput("pdb1", BSS.Gateway.FileSet(help="A coordinates file for the first molecule"))
node.addInput("pdb2", BSS.Gateway.FileSet(help="A coordinates file for the second molecule"))
node.addInput("xml1", BSS.Gateway.FileSet(help="A parameters file for the first molecule"))
node.addInput("xml2", BSS.Gateway.FileSet(help="A parameters file for the second molecule"))
node.addInput("prematch", BSS.Gateway.String(help="list of atom indices that are matched between input2 and input1. Syntax is of the format 1-3,4-8,9-11... Ignored if a mapping is provided", default=""))
node.addInput("mapping", BSS.Gateway.File(help="csv file that contains atom indices in input1 mapped ot atom indices in input2", optional=True))
node.addInput("timeout", BSS.Gateway.Time(help="The timeout for the maximum common substructure search", default=10*BSS.Units.Time.second))
node.addInput("allow_ring_breaking", BSS.Gateway.Boolean(help="Whether to allow opening/closing of rings during merge", default=False))
node.addInput("allow_ring_size_change", BSS.Gateway.Boolean(help="Whether to allow ring size changes during merge", default=False))
node.addInput("output", BSS.Gateway.String(help="The root name for the files describing the perturbation input1->input2."))


# In[ ]:


node.addOutput("nodeoutput", BSS.Gateway.FileSet(help="SOMD input files for a perturbation of input1->input2."))


# In[ ]:


node.showControls()

def writeLog(ligA, ligB, mapping):
    """ Human readable report on atoms used for the mapping."""
    atoms_in_A = list(mapping.keys())
    stream = open('somd.mapping','w')
    atAdone = []
    atBdone= []
    for atAidx in atoms_in_A:
        atA = ligA._sire_object.select(AtomIdx(atAidx))
        atB = ligB._sire_object.select(AtomIdx(mapping[atAidx]))
        stream.write("%s %s --> %s %s\n" % (atA.index(), atA.name(),atB.index(), atB.name()))
        atAdone.append(atA)
        atBdone.append(atB)
    for atom in ligA._sire_object.atoms():
        if atom in atAdone:
            continue
        stream.write("%s %s --> dummy\n" % (atom.index(), atom.name()))
    for atom in ligB._sire_object.atoms():
        if atom in atBdone:
            continue
        stream.write("dummy --> %s %s\n" % (atom.index(), atom.name()))
    stream.close()


# In[ ]:


def loadMapping(mapping_file):
    """Parse a text file that specifies mappings between atomic indices in input1 --> atoms in input2"""
    stream = open(mapping_file,'r')
    buffer = stream.readlines()
    stream.close()
    mapping = {}
    for line in buffer:
        if line.startswith("#"):
            continue
        elems = line.split(",")
        idx1 = int(elems[0])
        idx2 = int(elems[1])
        mapping[idx1] = idx2
    
    return mapping


#########################################
#       Config file parameters          #  
#########################################



combining_rules = Parameter("combining rules", "geometric",
                            """Combining rules to use for the non-bonded interactions.""")

cutoff_type = Parameter("cutoff type", "nocutoff", """The cutoff method to use during the simulation.""")

cutoff_dist = Parameter("cutoff distance", 12 * angstrom,
                        """The cutoff distance to use for the non-bonded interactions.""")

use_restraints = Parameter("use restraints", False, """Whether or not to use harmonic restaints on the solute atoms.""")


def createSystem(molecules):
    #print("Applying flexibility and zmatrix templates...")
    print("Creating the system...")

    moleculeNumbers = molecules.molNums()
    moleculeList = []

    for moleculeNumber in moleculeNumbers:
        molecule = molecules.molecule(moleculeNumber)[0].molecule()
        moleculeList.append(molecule)

    molecules = MoleculeGroup("molecules")
    ions = MoleculeGroup("ions")

    for molecule in moleculeList:
        natoms = molecule.nAtoms()
        if natoms == 1:
            ions.add(molecule)
        else:
            molecules.add(molecule)

    all = MoleculeGroup("all")
    all.add(molecules)
    all.add(ions)

    # Add these groups to the System
    system = System()

    system.add(all)
    system.add(molecules)
    system.add(ions)

    return system


def setupForcefields(system, space):

    print("Creating force fields... ")

    all = system[MGName("all")]
    molecules = system[MGName("molecules")]
    ions = system[MGName("ions")]

    # - first solvent-solvent coulomb/LJ (CLJ) energy
    internonbondedff = InterCLJFF("molecules:molecules")
    if (cutoff_type.val != "nocutoff"):
        internonbondedff.setUseReactionField(True)
        internonbondedff.setReactionFieldDielectric(rf_dielectric.val)
    internonbondedff.add(molecules)

    inter_ions_nonbondedff = InterCLJFF("ions:ions")
    if (cutoff_type.val != "nocutoff"):
        inter_ions_nonbondedff.setUseReactionField(True)
        inter_ions_nonbondedff.setReactionFieldDielectric(rf_dielectric.val)

    inter_ions_nonbondedff.add(ions)

    inter_ions_molecules_nonbondedff = InterGroupCLJFF("ions:molecules")
    if (cutoff_type.val != "nocutoff"):
        inter_ions_molecules_nonbondedff.setUseReactionField(True)
        inter_ions_molecules_nonbondedff.setReactionFieldDielectric(rf_dielectric.val)

    inter_ions_molecules_nonbondedff.add(ions, MGIdx(0))
    inter_ions_molecules_nonbondedff.add(molecules, MGIdx(1))

    # Now solute bond, angle, dihedral energy
    intrabondedff = InternalFF("molecules-intrabonded")
    intrabondedff.add(molecules)

    # Now solute intramolecular CLJ energy
    intranonbondedff = IntraCLJFF("molecules-intranonbonded")

    if (cutoff_type.val != "nocutoff"):
        intranonbondedff.setUseReactionField(True)
        intranonbondedff.setReactionFieldDielectric(rf_dielectric.val)

    intranonbondedff.add(molecules)

    # solute restraint energy
    #
    # We restrain atoms based ont he contents of the property "restrainedatoms"
    #
    restraintff = RestraintFF("restraint")

    if use_restraints.val:
        molnums = molecules.molecules().molNums()

        for molnum in molnums:
            mol = molecules.molecule(molnum)[0].molecule()
            try:
                mol_restrained_atoms = propertyToAtomNumVectorList(mol.property("restrainedatoms"))
            except UserWarning as error:
                error_type = re.search(r"(Sire\w*::\w*)", str(error)).group(0)
                if error_type == "SireBase::missing_property":
                    continue
                else:
                    raise error

            for restrained_line in mol_restrained_atoms:
                atnum = restrained_line[0]
                restraint_atom = mol.select(atnum)
                restraint_coords = restrained_line[1]
                restraint_k = restrained_line[2] * kcal_per_mol / (angstrom * angstrom)

                restraint = DistanceRestraint.harmonic(restraint_atom, restraint_coords, restraint_k)

                restraintff.add(restraint)

    # Here is the list of all forcefields
    forcefields = [internonbondedff, intrabondedff, intranonbondedff, inter_ions_nonbondedff,
                   inter_ions_molecules_nonbondedff, restraintff]

    for forcefield in forcefields:
        system.add(forcefield)

    system.setProperty("space", space)
    system.setProperty("switchingFunction", CHARMMSwitchingFunction(cutoff_dist.val))
    system.setProperty("combiningRules", VariantProperty(combining_rules.val))

    total_nrg = internonbondedff.components().total() + \
                intranonbondedff.components().total() + intrabondedff.components().total() + \
                inter_ions_nonbondedff.components().total() + inter_ions_molecules_nonbondedff.components().total() + \
                restraintff.components().total()

    e_total = system.totalComponent()

    system.setComponent(e_total, total_nrg)

    # Add a monitor that calculates the average total energy and average energy
    # deltas - we will collect both a mean average and an zwanzig average
    system.add("total_energy", MonitorComponent(e_total, Average()))

    return system


def vsiteListToProperty(list):
    prop = Properties()
    i = 0
    for entry in list:
        for key, value in entry.items():
            prop.setProperty("%s(%d)" % (key,i), VariantProperty(value))
        i += 1
    prop.setProperty("nvirtualsites",VariantProperty(i))
    return prop


def readXmlParameters(pdbfile, xmlfile):
# 1) Read a pdb file describing the system to simulate
    p = PDB2(pdbfile)
    s = p.toSystem()
    molecules = s.molecules()
    #print (molecules)
    with open (pdbfile, "r") as f:
        first_line = f.readline()
        if first_line.split()[0] == "CRYST1" :
            print (first_line)
            pbc_x = float(line.split()[1])
            pbc_y = float(line.split()[2])
            pbc_z = float(line.split()[3])
            space = PeriodicBox(Vector(pbc_x, pbc_y, pbc_z))
        else:
            space = Cartesian()
    print("space:", space)
#     system = System()

     # 2) Now we read the xml file, and store parameters for each molecule


    import xml.dom.minidom as minidom
    xmldoc = minidom.parse(xmlfile)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: TYPE ~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_type = xmldoc.getElementsByTagName('Type')
    dicts_type = []
    for items in itemlist_type:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_type.append(d)
    dicts_tp =  str(dicts_type).split()
    #print (dicts_tp)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: ATOM ~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_atom = xmldoc.getElementsByTagName('Atom')
    dicts_atom = []
    for items in itemlist_atom:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_atom.append(d)
    dicts_at =  str(dicts_atom).split()
    #print (dicts_at)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: BOND ~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_bond = xmldoc.getElementsByTagName('Bond')
    dicts_bond = []
    for items in itemlist_bond:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_bond.append(d)
    dicts_b =  str(dicts_bond).split()
    #print (dicts_b)

    nbond = itemlist_bond.length

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: ANGLE ~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_angle = xmldoc.getElementsByTagName('Angle')
    dicts_angle = []
    for items in itemlist_angle:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_angle.append(d)
    dicts_ang =  str(dicts_angle).split()
    print (dicts_angle)

    nAngles= itemlist_angle.length
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: PROPER ~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_proper = xmldoc.getElementsByTagName('Proper')
    dicts_proper = []
    for items in itemlist_proper:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_proper.append(d)
    dicts_pr =  str(dicts_proper).split()
    #print (dicts_pr)

    nProper = itemlist_proper.length

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: IMPROPER ~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_improper = xmldoc.getElementsByTagName('Improper')
    dicts_improper = []
    for items in itemlist_improper:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_improper.append(d)
    dicts_impr =  str(dicts_improper).split()
    #print (dicts_impr)
    nImproper = itemlist_improper.length

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: VIRTUAL SITES ~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_VirtualSite = xmldoc.getElementsByTagName('VirtualSite')
    dicts_virtualsite = []
    for items in itemlist_VirtualSite:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_virtualsite.append(d)
    #dicts_vs =  str(dicts_virtualsite).split()
    #print (dicts_vs)
    nVirtualSites = itemlist_VirtualSite.length 


    v_site_CLJ = []
    for i in range(0, int(len(dicts_atom))):
        if dicts_atom[i]['type'][0] == 'v':
            v_site_CLJ = dicts_atom[i]
            dicts_virtualsite.append(v_site_CLJ)

    for i in range(0, len(itemlist_VirtualSite)):
        dicts_virtualsite[i].update(dicts_virtualsite[i+len(itemlist_VirtualSite)])
        dicts_virtualsite[i].update(dicts_virtualsite[i+2*len(itemlist_VirtualSite)]) 

    dict_vs = []
    for i in range(0, len(itemlist_VirtualSite)):
        dicts_virtualsite[i]
        dict_vs.append(dicts_virtualsite[i])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: RESIDUE ~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_residue = xmldoc.getElementsByTagName('Residue')
    dicts_residue = []
    for items in itemlist_residue:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_residue.append(d)
    dicts_res =  str(dicts_residue).split()
    #print (dicts_res)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~ TAG NAME: NON BONDED FORCE ~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_nonbond = xmldoc.getElementsByTagName('NonbondedForce')
    dicts_nonb = []
    for items in itemlist_nonbond:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_nonb.append(d)
    dicts_nb =  str(dicts_nonb).split()
    #print (dicts_nb)
    nNonBonded = itemlist_nonbond.length

    # 3) Now we create an Amberparameters object for each molecule
    molnums = molecules.molNums()

    newmolecules = Molecules()
    for molnum in molnums:
        mol = molecules.at(molnum)
        #print (mol) 
        


        # Add potential virtual site parameters
        if len(dicts_virtualsite) > 0:
            mol = mol.edit().setProperty("virtual-sites", vsiteListToProperty(dict_vs)).commit()

            
        # We populate the Amberparameters object with a list of bond, angle, dihedrals
        # We look up parameters from the contents of the xml file
        # We also have to set the atomic parameters (q, sigma, epsilon)

        editmol = mol.edit()
        mol_params = AmberParameters(editmol) #SireMol::AmberParameters()
        atoms = editmol.atoms()
        # We update atom parameters see setAtomParameters in SireIO/amber.cpp l2122 
        natoms = editmol.nAtoms()
        #print("number of atoms is %s" %natoms)
        
    #natoms don't include the virtual sites! 

    # Loop over each molecule in the molecules object


        opls=[]
        for i in range (0, int(len(dicts_atom)/2)): 
            opl={} 
            opl = dicts_atom[i]['type'] 
            opls.append(opl) 

        name=[]
        for i in range (0, int(len(dicts_atom)/2)): 
            nm={} 
            nm = dicts_atom[i]['name'] 
            name.append(nm) 


        two=[] 
        #print(len(name)) 
        for i in range(0, len(name)): 
            t=(opls[i],name[i]) 
            two.append(t) 

        import numpy as np 

        atom_sorted = []
        for j in range(0, len(two)): 
            for i in range(int(len(dicts_atom)/2), len(dicts_atom)):   
                if dicts_atom[i]['type'] == two[j][0]: 
                    dic_a = {}
                    dic_a = dicts_atom[i]
                    atom_sorted.append(dic_a)
      
        type_sorted = []
        for j in range(0, len(two)): 
            for i in range(0, int(len(dicts_type))):   
                if dicts_type[i]['name'] == two[j][0]: 
                    dic_t = {}
                    dic_t = dicts_type[i]
                    type_sorted.append(dic_t)
        print(" ")
        print("There are ",natoms," atoms in this molecule. ")
        print("*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*")

        for atom in atoms: 
            editatom = editmol.atom(atom.index())

            i = int(str(atom.number()).split('(')[1].replace(")" , " ")) 

            editatom.setProperty("charge", float(atom_sorted[i-1]['charge']) * mod_electron)
            editatom.setProperty("mass", float(type_sorted[i-1]['mass']) * g_per_mol) 
            editatom.setProperty("LJ", LJParameter( float(atom_sorted[i-1]['sigma'])*10 * angstrom , float(atom_sorted[i-1]['epsilon'])/4.184 * kcal_per_mol))
            editatom.setProperty("ambertype", dicts_atom[i-1]['type'])
           
            editmol = editatom.molecule()

        # Now we create a connectivity see setConnectivity in SireIO/amber.cpp l2144
        # XML data tells us how atoms are bonded in the molecule (Bond 'from' and 'to')
        
        if natoms > 1:
            print("Set up connectivity")

            con = []
            for i in range(0,int(nbond/2)):
                if natoms > 1:   
                    connect_prop= {}
                    connect_prop = dicts_bond[i]['from'], dicts_bond[i]['to']
                con.append(connect_prop)
            

            conn = Connectivity(editmol.info()).edit()

            for j in range(0,len(con)):
                conn.connect(atoms[int(con[j][0]) ].index(), atoms[int(con[j][1]) ].index()) 
                   
            
            editmol.setProperty("connectivity", conn.commit()).commit()
            mol = editmol.setProperty("connectivity", conn.commit()).commit()
            system.update(mol)

             # Now we add bond parameters to the Sire molecule. We also update amberparameters see SireIO/amber.cpp l2154

            internalff = InternalFF()
                    
            bondfuncs = TwoAtomFunctions(mol)
            r = internalff.symbols().bond().r()

            for j in range(0,len(con)):
                bondfuncs.set(atoms[int(con[j][0]) ].index(), atoms[int(con[j][1]) ].index(), float(dicts_bond[j+len(con)]['k'])/(2*100*4.184)* (float(dicts_bond[j+len(con)]['length'])*10 - r) **2  )
                bond_id = BondID(atoms[int(con[j][0])].index(), atoms[int(con[j][1])].index())
                mol_params.add(bond_id, float(dicts_bond[j+len(con)]['k'])/(2*100*4.184), float(dicts_bond[j+len(con)]['length'])*10 )   
                editmol.setProperty("bonds", bondfuncs).commit()
                molecule = editmol.commit()
             

            mol_params.getAllBonds() 

            editmol.setProperty("amberparameters", mol_params).commit() # Weird, should work - investigate ? 
            molecule = editmol.commit()

        # Now we add angle parameters to the Sire molecule. We also update amberparameters see SireIO/amber.cpp L2172
        if natoms > 2:
            print("Set up angles")

            anglefuncs = ThreeAtomFunctions(mol)
            at1 = []
            for i in range(0, nAngles):
                a1 = {}
                to_str1 = str(re.findall(r"\d+",str(dicts_angle[i]['class1'])))
                if dicts_atom[i]['type'][0] == 'o': #if opls_
                    a1 = int(to_str1.replace("[","").replace("]","").replace("'","") )-800

                    at1.append(a1)
          
            print(at1)
            at2 = []
            for i in range(0, nAngles):
                a2 = {}
                to_str2 = str(re.findall(r"\d+",str(dicts_angle[i]['class2'])))
                if dicts_atom[i]['type'][0] == 'o': #if opls_
                    a2 = int(to_str2.replace("[","").replace("]","").replace("'","") )-800
                    at2.append(a2)
         

            at3 = []
            for i in range(0, nAngles):
                a3 = {}
                to_str3 = str(re.findall(r"\d+",str(dicts_angle[i]['class3'])))
                if dicts_atom[i]['type'][0] == 'o': #if opls_
                    a3 = int(to_str3.replace("[","").replace("]","").replace("'","") )-800

                    at3.append(a3)
            print ("lengths:", len(at1), len(at2), len(at3))
            print("number of Angles =",nAngles)
            theta = internalff.symbols().angle().theta()
            for j in range(0,nAngles- nVirtualSites):

                anglefuncs.set( atoms[at1[j]].index(), atoms[at2[j]].index(), atoms[at3[j]].index(), float(dicts_angle[j]['k'])/(2*4.184) * ( (float(dicts_angle[j]['angle']) - theta )**2 ))
                angle_id = AngleID( atoms[int(at1[j])].index(), atoms[int(at2[j])].index(), atoms[int(at3[j])].index())
                mol_params.add(angle_id, float(dicts_angle[j]['k'])/(2*4.184), float(dicts_angle[j]['angle']) ) 
            
        # Now we add dihedral parameters to the Sire molecule. We also update amberparameters see SireIO/amber.cpp L2190

        if natoms > 3:
            print("Set up dihedrals")
            di1 = []

            for i in range(0, nProper):
                d1 = {}
                to_str1 = str(re.findall(r"\d+",str(dicts_proper[i]['class1'])))
                if dicts_atom[0]['type'][0] == 'o':#if opls_
                    d1 = int(to_str1.replace("[","").replace("]","").replace("'","") )-800
                elif dicts_atom[i]['type'][0] == 'Q':  #if QUBE_
                    d1 = int(to_str1.replace("[","").replace("]","").replace("'","") )
               
                di1.append(d1)


            di2 = []
            for i in range(0, nProper):
                d2 = {}
                to_str2 = str(re.findall(r"\d+",str(dicts_proper[i]['class2'])))
                if dicts_atom[0]['type'][0] == 'o':#if opls_
                    d2 = int(to_str2.replace("[","").replace("]","").replace("'","") )-800

                    di2.append(d2)

            di3 = []
            for i in range(0, nProper):
                d3 = {}
                to_str3 = str(re.findall(r"\d+",str(dicts_proper[i]['class3'])))
                if dicts_atom[0]['type'][0] == 'o':#if opls_
                    d3 = int(to_str3.replace("[","").replace("]","").replace("'","") )-800

                    di3.append(d3)
          

            di4 = []
            for i in range(0, nProper):
                d4 = {}
                to_str4 = str(re.findall(r"\d+",str(dicts_proper[i]['class4'])))
                if dicts_atom[0]['type'][0] == 'o':#if opls_
                    d4 = int(to_str4.replace("[","").replace("]","").replace("'","") )-800

                    di4.append(d4)
 

            dihedralfuncs = FourAtomFunctions(mol)
    
            phi = internalff.symbols().dihedral().phi()
            for i in range(0,nProper):  
                if atoms[int(di1[i])].index() != atoms[int(di4[i])].index():
                    dihedral_id = DihedralID( atoms[int(di1[i])].index(), atoms[int(di2[i])].index(), atoms[int(di3[i])].index(), atoms[int(di4[i])].index()) 
                    dih1= float(dicts_proper[i]['k1'])/4.184*(1+Cos(int(dicts_proper[i]['periodicity1'])* phi- float(dicts_proper[i]['phase1'])))
                    dih2= float(dicts_proper[i]['k2'])/4.184*(1+Cos(int(dicts_proper[i]['periodicity2'])* phi- float(dicts_proper[i]['phase2'])))
                    dih3= float(dicts_proper[i]['k3'])/4.184*(1+Cos(int(dicts_proper[i]['periodicity3'])* phi- float(dicts_proper[i]['phase3'])))
                    dih4= float(dicts_proper[i]['k4'])/4.184*(1+Cos(int(dicts_proper[i]['periodicity4'])* phi- float(dicts_proper[i]['phase4'])))
                    dih_fun = dih1 + dih2 +dih3 +dih4
                    dihedralfuncs.set(dihedral_id, dih_fun)

                    for t in range(1,5):
                        mol_params.add(dihedral_id, float(dicts_proper[i]['k%s'%t])/4.184, int(dicts_proper[i]['periodicity%s'%t]), float(dicts_proper[i]['phase%s'%t]) ) 
            

            print("Set up impropers")

            di_im1 = []
            for i in range(0, nImproper):
                d1 = {}
                to_str1 = str(re.findall(r"\d+",str(dicts_improper[i]['class1'])))
                if dicts_atom[0]['type'][0] == 'o':#if opls_
                    d1 = int(to_str1.replace("[","").replace("]","").replace("'","") )-800

                    di_im1.append(d1)


            di_im2 = []
            for i in range(0, nImproper):
                d2 = {}
                to_str2 = str(re.findall(r"\d+",str(dicts_improper[i]['class2'])))
                if dicts_atom[0]['type'][0] == 'o':#if opls_
                    d2 = int(to_str2.replace("[","").replace("]","").replace("'","") )-800
                elif dicts_atom[i]['type'][0] == 'Q':
                    d2 = int(to_str2.replace("[","").replace("]","").replace("'","") )

                
                di_im2.append(d2)

            di_im3 = []
            for i in range(0, nImproper):
                d3 = {}
                to_str3 = str(re.findall(r"\d+",str(dicts_improper[i]['class3'])))
                if dicts_atom[0]['type'][0] == 'o':#if opls_
                    d3 = int(to_str3.replace("[","").replace("]","").replace("'","") )-800
                elif dicts_atom[i]['type'][0] == 'Q':
                    d3 = int(to_str3.replace("[","").replace("]","").replace("'","") )
                
                di_im3.append(d3)

            di_im4 = []
            for i in range(0, nImproper):
                d4 = {}
                to_str4 = str(re.findall(r"\d+",str(dicts_improper[i]['class4'])))
                if dicts_atom[0]['type'][0] == 'o':#if opls_
                    d4 = int(to_str4.replace("[","").replace("]","").replace("'","") )-800
                elif dicts_atom[i]['type'][0] == 'Q':
                    d4 = int(to_str4.replace("[","").replace("]","").replace("'","") )
                
                di_im4.append(d4)
            
            improperfuncs = FourAtomFunctions(mol)

            phi_im = internalff.symbols().improper().phi()


            for i in range(0,nImproper):  
                improper_id = ImproperID( atoms[int(di_im2[i])].index(), atoms[int(di_im3[i])].index(), atoms[int(di_im1[i])].index(), atoms[int(di_im4[i])].index()) 
                imp1= float(dicts_improper[i]['k1'])*(1/4.184)*(1+Cos(int(dicts_improper[i]['periodicity1'])* phi_im - float(dicts_improper[i]['phase1'])))
                imp2= float(dicts_improper[i]['k2'])*(1/4.184)*(1+Cos(int(dicts_improper[i]['periodicity2'])* phi_im - float(dicts_improper[i]['phase2'])))
                imp3= float(dicts_improper[i]['k3'])*(1/4.184)*(1+Cos(int(dicts_improper[i]['periodicity3'])* phi_im - float(dicts_improper[i]['phase3'])))
                imp4= float(dicts_improper[i]['k4'])*(1/4.184)*(1+Cos(int(dicts_improper[i]['periodicity4'])* phi_im - float(dicts_improper[i]['phase4'])))
                imp_fun = imp1 + imp2 +imp3 +imp4
                improperfuncs.set(improper_id, imp_fun)
                #print(improperfuncs.potentials())

                for t in range(1,5):
                    mol_params.add(improper_id, float(dicts_improper[i]['k%s'%t])*(1/4.184), int(dicts_improper[i]['periodicity%s'%t]), float(dicts_improper[i]['phase%s'%t]) ) 


            mol = editmol.setProperty("bond", bondfuncs).commit()
            mol = editmol.setProperty("angle" , anglefuncs).commit()
            mol = editmol.setProperty("dihedral" , dihedralfuncs).commit()
            mol = editmol.setProperty("improper" , improperfuncs).commit()
            system.update(mol)

        # Now we work out non bonded pairs see SireIO/amber.cpp L2213


            print("Set up nbpairs")
            print("*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*")
            ## Define the bonded pairs in a list that is called are12
            #print("Now calculating 1-2 intercactions")
            are12 = []
            for i in range(0, natoms): 
                for j in range (0, natoms): 
                    if conn.areBonded(atoms[i].index(), atoms[j].index()) == True:
                        #ij = {}
                        ij= (i,j)
                        are12.append(ij)
            are12_bckup = are12[:]


            #print("Now calculating 1-3 intercactions")
            are13 = []
            for i in range(0, natoms): 
                for j in range (0, natoms): 
                    if conn.areAngled(atoms[i].index(), atoms[j].index()) == True:
                        ij = {}
                        ij= (i,j)
                        are13.append(ij)
            are13_bckup = are13[:]

           # print("Now calculating 1-4 intercactions")
            are14 = []
            for i in range(0, natoms): 
                for j in range (0, natoms):
                   
                    if conn.areDihedraled(atoms[i].index(), atoms[j].index()) == True and conn.areAngled(atoms[i].index(), atoms[j].index()) == False:
                        ij = {}
                        ij= (i,j)
                        are14.append(ij)
            are14_bckup = are14[:]

           # print("Now calculating the non-bonded intercactions")
            bonded_pairs_list = are12_bckup + are13_bckup + are14_bckup    
            nb_pair_list =[]

            for i in range(0, natoms): 
                #print("i=",i)
                for j in range (0, natoms):
                    if i != j and (i,j) not in bonded_pairs_list:
                        nb_pair_list.append((i,j))
            are_nb_bckup = nb_pair_list[:]

            nbpairs = CLJNBPairs(editmol.info(), CLJScaleFactor(0,0))
            #print("Now setting 1-2 intercactions")
            for i in range(0, len(are12)):
                scale_factor1 = 0
                scale_factor2 = 0
                nbpairs.set(atoms.index( int(are12[i][0])), atoms.index(int(are12[i][1])), CLJScaleFactor(scale_factor1,scale_factor2))

            #print("Now setting 1-3 intercactions")
            for i in range(0, len(are13)):
                scale_factor1 = 0
                scale_factor2 = 0
                nbpairs.set(atoms.index( int(are13[i][0])), atoms.index(int(are13[i][1])), CLJScaleFactor(scale_factor1,scale_factor2))

           # print("Now setting 1-4 intercactions")              
            for i in range(0, len(are14)): 
                scale_factor1 = 1/2
                scale_factor2 = 1/2
                nbpairs.set(atoms.index( int(are14[i][0])), atoms.index(int(are14[i][1])), CLJScaleFactor(scale_factor1,scale_factor2))
                mol_params.add14Pair(BondID(atoms.index( int(are14[i][0])), atoms.index( int(are14[i][1]))),scale_factor1 , scale_factor2)

           # print("Now setting non-bonded intercactions")  
            #print("*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*")              
            for i in range(0, len(nb_pair_list)):  
                scale_factor1 = 1
                scale_factor2 = 1
                nbpairs.set(atoms.index( int(nb_pair_list[i][0])), atoms.index(int(nb_pair_list[i][1])), CLJScaleFactor(scale_factor1,scale_factor2))

                # print("~~~~~~~~~~~~~~~~~~`")
                
            mol = editmol.setProperty("intrascale" , nbpairs).commit()
            system.update(mol)


        #print("Setup name of qube FF")
        # mol = mol.edit().setProperty("forcefield", ffToProperty("qube")).commit()
        # system.update(mol)

        molecule = editmol.commit()
        newmolecules.add(molecule)



    return (newmolecules, space)


# In[ ]:

do_mapping = True
custom_mapping = node.getInput("mapping")
#print (custom_mapping)
if custom_mapping is not None:
    do_mapping = False
    mapping = loadMapping(custom_mapping)
    #print (mapping)




# Optional input, dictionary of Atom indices that should be matched in the search. 
prematch = {}
prematchstring = node.getInput("prematch")
if len(prematchstring) > 0: 
    entries = prematchstring.split(",")
    for entry in entries:
        idxA, idxB = entry.split("-")
        prematch[int(idxA)] = int(idxB)
#print (prematch)




# Load system 1

xmlfile = node.getInput("xml1")[0]
pdbfile = node.getInput("pdb1")[0]
(molecules, space) = readXmlParameters(pdbfile, xmlfile) 




# Load system 2

xmlfile1 = node.getInput("xml2")[0]
pdbfile1 = node.getInput("pdb2")[0]
(molecules1, space1) = readXmlParameters(pdbfile1, xmlfile1) 


# We assume the molecules to perturb are the first molecules in each system
mol1 = molecules.first()
mol2 = molecules1.first()

#Convert the input to BioSimSpace molecules
bss_mol1 = BSS._SireWrappers.Molecule(mol1)
bss_mol2 = BSS._SireWrappers.Molecule(mol2)

if do_mapping:
    # Return a maximum of 10 matches, scored by RMSD and sorted from best to worst.
    mappings, scores = BSS.Align.matchAtoms(bss_mol1, bss_mol2, matches=10, prematch=prematch, return_scores=True, scoring_function="RMSDalign", timeout=node.getInput("timeout"))
    # We retain the top mapping
    mapping = mappings[0]
    #print (len(mappings))
    #print (mappings)



#print (mapping)
#for x in range(0,len(mappings)):
#    print (mappings[x], scores[x])


inverted_mapping = dict([[v,k] for k,v in mapping.items()])
#print (inverted_mapping)


# Align lig2 to lig1 based on the best mapping (inverted). The molecule is aligned based
# on a root mean squared displacement fit to find the optimal translation vector
# (as opposed to merely taking the difference of centroids).
lig1 = bss_mol1
lig2 = BSS.Align.rmsdAlign(bss_mol2, lig1, inverted_mapping)  

# Merge the two ligands based on the mapping.
merged = BSS.Align.merge(lig1, lig2, mapping, allow_ring_breaking=node.getInput("allow_ring_breaking"), allow_ring_size_change=node.getInput("allow_ring_size_change"))

# Create a composite system
system1 = BSS._SireWrappers.System(BSS._SireWrappers.Molecules(molecules))
system1.removeMolecules(bss_mol1) 
system1.addMolecules(merged)


# Log the mapping used
writeLog(lig1, lig2, mapping)
BSS.IO.saveMolecules("merged_at_lam0.pdb", system1, "PDB", { "coordinates" : "coordinates0" , "element": "element0" })
# Generate package specific input
protocol = BSS.Protocol.FreeEnergy(runtime = 2*BSS.Units.Time.femtosecond, num_lam=3)
process = BSS.Process.Somd(system1, protocol)
process.getOutput()
with zipfile.ZipFile("somd_output.zip", "r") as zip_hnd:
    zip_hnd.extractall(".")


# In[ ]:


root = node.getInput("output")
mergedpdb = "%s.mergeat0.pdb" % root
pert = "%s.pert" % root
prm7 = "%s.prm7" % root
rst7 = "%s.rst7" % root
mapping_str = "%s.mapping" % root


os.replace("merged_at_lam0.pdb", mergedpdb)
os.replace("somd.pert", pert)
os.replace("somd.prm7", prm7)
os.replace("somd.rst7", rst7)
os.replace("somd.mapping", mapping_str)
try:
    os.remove("somd_output.zip")
    os.remove("somd.cfg")
    os.remove("somd.err")
    os.remove("somd.out")
except Exception:
    pass

node.setOutput("nodeoutput",[mergedpdb, pert, prm7, rst7, mapping_str])


# In[ ]:


node.validate()
