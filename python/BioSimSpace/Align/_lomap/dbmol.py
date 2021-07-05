# ******************
# MODULE DOCSTRING
# ******************

"""

LOMAP
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

import argparse
import glob
import logging
import math
import multiprocessing
import os
import pickle
from ._version import get_versions

import networkx as nx
import numpy as np
from . import graphgen
from . import mcs
from rdkit import Chem
from rdkit import DataStructs

__all__ = ['DBMolecules', 'SMatrix', 'Molecule']


# *************************
# Molecule Database Class
# *************************

class DBMolecules(object):
    """
    
    This class is used as a container for all the Molecules  
    
    """

    # Initialization function
    def __init__(self, directory, parallel=1, verbose='off',
                 time=20, ecrscore=0.0, threed=False, max3d=1000.0, output=False,
                 name='out', output_no_images=False, output_no_graph=False, display=False,
                 allow_tree=False, max=6, cutoff=0.4, radial=False, hub=None, fast=False, 
                 links_file=None, known_actives_file=None, max_dist_from_actives=2):

        """
        Initialization of  the Molecule Database Class
    
        Parameters
        ----------
        directory : str 
           the mol2/sdf directory file name
        parallel : int
           the number of cores used to generate the similarity score matrices
        verbose : bool
           verbose mode
        time : int
           the maximum time in seconds used to perform the MCS search
        ecrscore: float
           the electrostatic score to be used (if != 0) if two molecule have diffrent charges
        threed: bool
           If true, symmetry-equivalent MCSes are filtered to prefer the one with the best real-space alignment
        max3d: float
           The MCS is filtered to remove atoms which are further apart than this threshold. The default of 1000 is 
           effectively "no filter"
        output : bool
           a flag used to generate or not the output files
        name : str
           the file name prefix used to produce the output files
        output_no_images : bool
           a flag used to disable the generation of the output image files
        output_no_graph : bool
           a flag used to disable the generation of the output graph (.dot) file
        display : bool
           a flag used to display or not a network made by using matplotlib
        allow_tree: bool
           if set, then the final graph does not need a cycle covering and will be a tree
        max : int
           the maximum diameter of the resulting graph 
        cutoff : float
           the Minimum Similarity Score (MSS) used to build the graph
        links_file : str
           the name of a file containing links to seed the graph with
        known_actives_file : str
           the name of a file containing mols whose activity is known
        max_dist_from_actives : int
            The maximum number of links from any molecule to an active
            

        """

        # Set the Logging 
        if verbose == 'off':
            logging.basicConfig(format='%(message)s', level=logging.CRITICAL)

        if verbose == 'info':
            logging.basicConfig(format='%(message)s', level=logging.INFO)

        if verbose == 'pedantic':
            logging.basicConfig(format='%(message)s', level=logging.DEBUG)
            # logging.basicConfig(format='%(levelname)s:\t%(message)s', level=logging.DEBUG)

        if __name__ == '__main__':
            self.options = parser.parse_args()

        else:

            if not isinstance(output, bool):
                raise TypeError('The output flag is not a bool type')

            if not isinstance(output_no_images, bool):
                raise TypeError('The output_no_images flag is not a bool type')

            if not isinstance(output_no_graph, bool):
                raise TypeError('The output_no_graph flag is not a bool type')

            if not isinstance(display, bool):
                raise TypeError('The display flag is not a bool type')

            if not isinstance(radial, bool):
                raise TypeError('The radial flag is not a bool type')
            output_str = ''
            output_no_images_str = ''
            output_no_graph_str = ''
            display_str = ''
            radial_str = ''
            fast_str = ''
            threed_str = ''
            links_file_str = ''
            known_actives_file_str = ''
            allow_tree_str = ''

            if output:
                output_str = '--output'

            if output_no_images:
                output_no_images_str = '--output-no-images'

            if output_no_graph:
                output_no_graph_str = '--output-no-graph'

            if display:
                display_str = '--display'

            if radial:
                radial_str = '--radial'

            if fast:
                fast_str = '--fast'

            if threed:
                threed_str = '--threed'

            if allow_tree:
                allow_tree_str = '--allow-tree'

            if links_file:
                links_file_str = f'--links-file {links_file}'

            if known_actives_file:
                known_actives_file_str = f'--known-actives-file {known_actives_file}'

            names_str = '%s --parallel %s --verbose %s --time %s --ecrscore %s --max3d %s --name %s --max %s --max-dist-from-actives %s --cutoff %s --hub %s %s %s %s %s %s %s %s %s %s %s' \
                        % (
                        directory, parallel, verbose, time, ecrscore, max3d, name, max, max_dist_from_actives, cutoff, hub, output_str, display_str, output_no_images_str, output_no_graph_str,
                        radial_str, fast_str, threed_str, allow_tree_str, links_file_str, known_actives_file_str)

            #print("ARGS:",names_str)
            self.options = parser.parse_args(names_str.split())

        # Internal list container used to store the loaded molecule objects
        self.__list = self.read_molecule_files()

        # Dictionary which holds the mapping between the generated molecule IDs and molecule file names
        self.dic_mapping = {}
        self.inv_dic_mapping = {}

        # Hold the MCS index map strings for each molecule pair. Indexed by a tuple of molecule IDs (lowest first)
        self.mcs_map_store = {}

        # Pre-specified links between molecules - a map of molecule index tuples to score. 
        # A value < -1 means "recompute, but force the link to be included"
        # A negative value (>= -1) means "Use the absolute value as the score, but force the link to be included"
        # A positive value means "Use this value as the score, but treat the link as normal in the graph calculation"
        self.prespecified_links = {}

        # List of which molecules are "known actives". Note that all pairs of known actives 
        # are automatically added as prespecified links with a score of -1 (i.e force score to
        # 1 and force link to be included)
        self.known_actives = []

        for mol in self.__list:
            self.dic_mapping[mol.getID()] = mol.getName()
            self.inv_dic_mapping[mol.getName()] = mol.getID()

        if self.options.links_file and len(self.options.links_file)>0:
            self.parse_links_file(self.options.links_file)

        if self.options.known_actives_file and len(self.options.known_actives_file)>0:
            self.parse_known_actives_file(self.options.known_actives_file)

        # Index used to perform index selection by using __iter__ function
        self.__ci = 0

        # Symmetric matrices used to store the mcs scoring. The matrices are subclasses of numpy 
        self.strict_mtx = SMatrix(shape=(0,))
        self.loose_mtx = SMatrix(shape=(0,))

        # Empty pointer to the networkx graph 
        self.Graph = nx.Graph()

    def __iter__(self):
        """
        Index generator
        """
        return self

    def next(self):  # Python 3: def __next__(self)
        """
        Select the molecule during an iteration
        """

        if self.__ci > len(self.__list) - 1:
            self.__ci = 0
            raise StopIteration
        else:
            self.__ci = self.__ci + 1
            return self.__list[self.__ci - 1]

    def __getitem__(self, index):
        """
        Slicing and index selection function
        """

        return self.__list[index]

    def __setitem__(self, index, molecule):
        """
        Index setting function
        
        Parameters
        ----------
        index : int 
           the molecule index
        molecule : Molecule obj
           the molecule to assign to the molecule database by selecting the index: 
           DB[index] = molecule
        
        """

        if not isinstance(molecule, Molecule):
            raise ValueError('The passed molecule is not a Molecule object')

        self.__list[index] = molecule

    def __add__(self, molecule):

        """
        Add a new molecule to the molecule database    
        
        Parameters
        ----------
        molecule : Molecule obj 
           the molecule to append into the molecule database
        """

        if not isinstance(molecule, Molecule):
            raise ValueError('The passed molecule is not a Molecule object')

            self.__list.append(molecule)

    def nums(self):
        """
        This function recovers the total number of molecules currently stored in
        the molecule database
        """
        return len(self.__list)

    def read_molecule_files(self):
        """
        Read in all the mol2 or SDF files

        Returns
        -------
        molid_list : list of Molecule objects
           the container list of all the allocated Molecule objects

        """

        # This list is used as container to handle all the molecules read in by using RdKit.
        # All the molecules are instances of  Molecule class
        molid_list = []

        # List of molecule that failed to load in
        mol_error_list_fn = []

        logging.info(30 * '-')

        # The .mol2 and .sdf file formats are the only supported so far
        mol_fnames = glob.glob(self.options.directory + "/*.mol2")
        mol_fnames += glob.glob(self.options.directory + "/*.sdf")

        mol_fnames.sort()

        if len(mol_fnames) < 2:
            raise IOError('The directory %s must contain at least two mol2/sdf files' % self.options.directory)

        print_cnt = 0
        mol_id_cnt = 0

        for fname in mol_fnames:
            # The RDkit molecule object reads in as mol2/sdf file. The molecule is not sanitized and 
            # all the hydrogens are kept in place - we are assuming 3D input, correctly charged
            # and prepared in the protein active site
            if fname.endswith(".mol2"):
                rdkit_mol = Chem.MolFromMol2File(fname, sanitize=False, removeHs=False)
            else:
                rdkit_mol = Chem.MolFromMolFile(fname, sanitize=False, removeHs=False)

            # Reading problems
            if rdkit_mol == None:
                logging.warning('Error reading the file: %s' % os.path.basename(fname))
                mol_error_list_fn.append(os.path.basename(fname))
                continue

            # The Rdkit molecule is stored in a Molecule object
            mol = Molecule(rdkit_mol, mol_id_cnt, os.path.basename(fname))
            mol_id_cnt += 1

            # Cosmetic printing and status
            if print_cnt < 15 or print_cnt == (len(mol_fnames) - 1):
                logging.info('ID %s\t%s' % (mol.getID(), os.path.basename(fname)))

            if print_cnt == 15:
                logging.info('ID %s\t%s' % (mol.getID(), os.path.basename(fname)))
                logging.info(3 * '\t.\t.\n')

            print_cnt += 1

            molid_list.append(mol)

        logging.info(30 * '-')

        logging.info('Finish reading input files. %d structures in total....skipped %d\n' % (
        len(molid_list), len(mol_error_list_fn)))

        if mol_error_list_fn:
            logging.warning('Skipped molecules:')
            logging.warning(30 * '-')
            for fn in mol_error_list_fn:
                logging.warning('%s' % fn)
            print(30 * '-')

        return molid_list

    def parse_links_file(self, links_file):
        try:
            with open(links_file,"r") as lf:
                for line in lf:
                    mols = line.split();
                    if (len(mols)<2 or len(mols)>4):
                        raise IOError('Syntax error in links file parsing line:'+line) 
                    indexa = self.inv_dic_mapping[mols[0]]
                    indexb = self.inv_dic_mapping[mols[1]]
                    score = -2
                    if (len(mols)>2):
                        score = float(mols[2])
                    if (len(mols)>3):
                        if (mols[3] != "force"):
                            raise IOError('Syntax error parsing fourth argument in links file on line:'+line) 
                        score = -score
                    self.prespecified_links[(indexa,indexb)]=score
                    self.prespecified_links[(indexb,indexa)]=score
                    print("Added prespecified link for mols",mols,"->",(indexa,indexb),"score",score)
        except KeyError as e:
            raise IOError('Filename within the links file "'+links_file+'" not found: '+str(e)) from None

    def parse_known_actives_file(self, actives_file):
        try:
            with open(actives_file,"r") as lf:
                for line in lf:
                    mols = line.split();
                    indexa = self.inv_dic_mapping[mols[0]]
                    self.known_actives.append(indexa)
                    self.__list[indexa].setActive(True)
                    print("Added known activity for mol",mols[0],"->",indexa)
        except KeyError as e:
            raise IOError('Filename within the actives file "'+actives_file+'" not found: '+str(e)) from None
        # Add all combinations of these to the set of prespecified links
        for t in [(x,y) for x in self.known_actives for y in self.known_actives]:
            print("Added prespecified link for ",t)
            self.prespecified_links[t]=-1

    def set_MCSmap(self,i,j,MCmap):
        if (i<j):
            idx=(i,j)
        else:
            idx=(j,i)
        self.mcs_map_store[idx]=MCmap

    def get_MCSmap(self,i,j):
        if (i<j):
            idx=(i,j)
        else:
            idx=(j,i)
        if idx in self.mcs_map_store:
            return self.mcs_map_store[idx]
        return None

    def compute_mtx(self, a, b, strict_mtx, loose_mtx, true_strict_mtx, MCS_map):
        """
        Compute a chunk of the similarity score matrices. The chunk is selected
        by the start index a and the final index b. The matrices are indeed 
        treated as linear array

        Parameters
        ----------
        a : int 
           the start index of the chunk 
        b : int
           the final index of the chunk
        strict_mtx: python multiprocessing array

           strict similarity score matrix. This array is used as shared memory
           array managed by the different allocated processes. Each process 
           operates on a separate chunk selected by the indexes a and b

        loose_mtx: python multiprocessing array
           loose similarity score matrix. This array is used as shared memory 
           array managed by the different allocated processes. Each process 
           operates on a separate chunk selected by the indexes a and b

        true_strict_mtx: python multiprocessing array
           Holds the strict score *before* that is potentially
           modified by the prespecified link function (which sets the link score to 1.0).

        MCS_map: dict (multiprocessing)
            Holds a dict of (index tuple) -> string with the strings being the 
            MCS atom index map between the two molecules

        """

        # name = multiprocessing.current_process().name
        # print(name)
        # print('a = %d, b = %d' % (a,b))
        # print('\n') 

        def formal_charge(mol):
            total_charge_mol = 0.0

            try:
                # Assume mol2
                total_charge_mol=sum([float(a.GetProp('_TriposPartialCharge')) for a in mol.GetAtoms()])
            except:
                # wasn't mol2, so assume SDF with correct formal charge props for mols
                total_charge_mol=sum([a.GetFormalCharge() for a in mol.GetAtoms()])

            return total_charge_mol



        def ecr(mol_i, mol_j):
            """
            This function computes the similarity score between the passed molecules
            by using the EleCtrostatic Rule (ECR)

            Parameters
            ----------
            mol_i : Rdkit molecule object 
               the first molecules used to calculate the ECR rule  
            mol_j : Rdkit molecule object 
               the second molecules used to calculate the ECR rule 

            Returns
            -------
            scr_ecr: float
                the calculated similarity score (1 if mol_i and mol_j have the
                same total charges, 0  otherwire)

            """

            total_charge_mol_i = formal_charge(mol_i) 
            total_charge_mol_j = formal_charge(mol_j) 

            if abs(total_charge_mol_j - total_charge_mol_i) < 1e-3:
                scr_ecr = 1.0
            else:
                scr_ecr = 0.0

            return scr_ecr

        # Total number of loaded molecules
        n = self.nums()

        # Looping over all the elements of the selected matrix chunk
        for k in range(a, b + 1):

            # The linear index k is converted into the row and column indexes of
            # an hypothetical bidimensional symmetric matrix
            i = int(n - 2 - math.floor(math.sqrt(-8 * k + 4 * n * (n - 1) - 7) / 2.0 - 0.5))
            j = int(k + i + 1 - n * (n - 1) / 2 + (n - i) * ((n - i) - 1) / 2)
            # print 'k = %d , i = %d , j = %d' % (k,i,j)

            # The Rdkit molecules moli and molj are extracted form the molecule database
            moli = self[i].getMolecule()
            molj = self[j].getMolecule()

            logging.info('Processing molecules: %s-%s' % (self[i].getName(),self[j].getName()))

            # The Electrostatic score rule is calculated
            ecr_score = ecr(moli, molj)

            # If the prespecified links map has this link, and the value is >-1, then
            # we don't need to compute the score
            if (i,j) in self.prespecified_links and self.prespecified_links[(i,j)]>=-1:
                strict_scr = abs(self.prespecified_links[(i,j)])
                loose_scr = strict_scr
                logging.info('MCS molecules: %s - %s final score %s set in links file' %
                      (self[i].getName(), self[j].getName(), strict_scr))
            else:
                # The MCS is computed only if the passed molecules have the same charges 
                if ecr_score or self.options.ecrscore:
                    if ecr_score == 0.0 and self.options.ecrscore:
                        logging.critical('WARNING: Mutation between different charge molecules is enabled')
                        ecr_score = self.options.ecrscore

                    try:
                        if self.options.verbose == 'pedantic':
                            logging.info(50 * '-')
                            logging.info('MCS molecules: %s - %s' % (self[i].getName(), self[j].getName()))

                        # Maximum Common Subgraph (MCS) calculation
                        MC = mcs.MCS(moli, molj, options=self.options)
                        ml=MC.all_atom_match_list()
                        self.set_MCSmap(i,j,ml)
                        MCS_map[(i,j)]=ml

                    except Exception as e:
                        logging.warning(
                            'Skipping MCS molecules (exception): %s - %s\t\n\n%s' % (self[i].getName(), self[j].getName(), e))
                        logging.info(50 * '-')
                        continue
                else:
                    continue

                # The scoring between the two molecules is performed by using different rules.
                # The total score will be the product of all the single rules
                tmp_scr = ecr_score * MC.mncar() * MC.mcsr() * MC.atomic_number_rule() * MC.hybridization_rule()
                tmp_scr *= MC.sulfonamides_rule() * MC.heterocycles_rule() * MC.transmuting_methyl_into_ring_rule()
                tmp_scr *= MC.transmuting_ring_sizes_rule()
                # Note - no longer using tmcsr rule!
                strict_scr = tmp_scr * 1 #  MC.tmcsr(strict_flag=True)
                loose_scr = tmp_scr * 1 #  MC.tmcsr(strict_flag=False)
                logging.info(
                    'MCS molecules: %s - %s final score %s from ecr %s mncar %s mcsr %s tmcsr %s anum %s sulf %s het %s RingMe %s' % 
                      (self[i].getName(), self[j].getName(), strict_scr, ecr_score, MC.mncar(),MC.mcsr(),MC.tmcsr(strict_flag=True),
                        MC.atomic_number_rule(),MC.sulfonamides_rule(),MC.heterocycles_rule(),MC.transmuting_methyl_into_ring_rule()))

            strict_mtx[k] = strict_scr
            loose_mtx[k] = loose_scr
            true_strict_mtx[k] = strict_scr

            # process prespecified links now and overwrite the existing info
            if (i,j) in self.prespecified_links and self.prespecified_links[(i,j)]<0:
                print("Molecule pair",i,j,"forced to be included in the graph - score set to 1")
                strict_mtx[k] = 1.0
                loose_mtx[k] = 1.0
                # Note that true_strict_mtx holds the original strict_scr value
                continue

        return

    def build_matrices(self):
        """
        This function coordinates the calculation of the similarity score matrices
        by distributing chunks of the matrices between the allocated processes

        """

        logging.info('\nMatrix scoring in progress....\n')

        # The similarity score matrices are defined instances of the class SMatrix
        # which implements a basic class for symmetric matrices
        self.strict_mtx = SMatrix(shape=(self.nums(),))
        self.loose_mtx = SMatrix(shape=(self.nums(),))
        self.true_strict_mtx = SMatrix(shape=(self.nums(),))
        # The total number of the effective elements present in the symmetric matrix
        l = int(self.nums() * (self.nums() - 1) / 2)

        if self.options.parallel == 1:  # Serial execution
            MCS_map = {}
            self.compute_mtx(0, l - 1, self.strict_mtx, self.loose_mtx, self.true_strict_mtx, MCS_map)
            for idx in MCS_map:
              self.set_MCSmap(idx[0],idx[1],MCS_map[idx])
        else:
            # Parallel execution
            logging.info('Parallel mode is on')

            # Number of selected processes
            num_proc = self.options.parallel

            delta = int(l / num_proc)
            rem = l % num_proc

            if delta < 1:
                kmax = l
            else:
                kmax = num_proc
            proc = []

            with multiprocessing.Manager() as manager:

              # Shared memory array used by the different allocated processes
              # At the moment we're using a combination of Array and Manager, which is nasty
              strict_mtx = multiprocessing.Array('d', self.strict_mtx)
              loose_mtx = multiprocessing.Array('d', self.loose_mtx)
              true_strict_mtx = multiprocessing.Array('d', self.true_strict_mtx)
              MCS_map = manager.dict()
  
              # Chopping the indexes redistributing the remainder
              for k in range(0, kmax):
  
                  spc = delta + int(int(rem / (k + 1)) > 0)
  
                  if k == 0:
                      i = 0
                  else:
                      i = j + 1
  
                  if k != kmax - 1:
                      j = i + spc - 1
                  else:
                      j = l - 1
  
                  # Python multiprocessing allocation
                  p = multiprocessing.Process(target=self.compute_mtx,
                                              args=(i, j, strict_mtx, loose_mtx, true_strict_mtx, MCS_map, ))
                  p.start()
                  proc.append(p)
              # End parallel execution
              for p in proc:
                  p.join()

              # Copying back the results
              self.strict_mtx[:] = strict_mtx[:]
              self.loose_mtx[:] = loose_mtx[:]
              self.true_strict_mtx[:] = true_strict_mtx[:]
              for idx in MCS_map.keys():
                self.set_MCSmap(idx[0],idx[1],MCS_map[idx])

        return self.strict_mtx, self.loose_mtx

    def build_graph(self):
        """
        This function coordinates the Graph generation

        """
        logging.info('\nGenerating graph in progress....')

        # The Graph is build from an instance of the Class GraphGen by passing
        # the selected user options
        Gr = graphgen.GraphGen(self)

        # Writing the results is files
        if self.options.output:
            try:
                Gr.write_graph(self.options.output_no_images, self.options.output_no_graph)
                with open(self.options.name + ".pickle", "wb") as pickle_f:
                    pickle.dump(Gr, pickle_f)
            except Exception as e:
                logging.error(str(e))

        # Handle to the the NetworkX generated graph
        self.Graph = Gr.get_graph()

        # print self.Graph.nodes(data=True)

        # Display the graph by using Matplotlib
        if self.options.display:
            Gr.draw()

        return self.Graph

    def write_dic(self):
        """
        This function writes out a text file with the mapping between the
        generated molecule indexes and the corresponding molecule file names

        """

        try:
            file_txt = open(self.options.name + '.txt', 'w')
        except Exception:
            raise IOError('It was not possible to write out the mapping file')
        file_txt.write('#ID\tFileName\n')
        for key in self.dic_mapping:
            file_txt.write('%d\t%s\n' % (key, self.dic_mapping[key]))

        file_txt.close()

    # *************************


# Symmetric  Class
# *************************


class SMatrix(np.ndarray):
    """
    This class implements a "basic" interface for symmetric matrices 
    subclassing ndarray. The class internally stores a bi-dimensional
    numpy array as a linear array A[k], however the user can still 
    access to the matrix elements by using a two indeces notation A[i,j]
 
    """

    def __new__(subtype, shape, dtype=float, buffer=None, offset=0, strides=None, order=None):
        if len(shape) > 2:
            raise ValueError('The matrix shape is greater than two')

        elif len(shape) == 2:
            if shape[0] != shape[1]:
                raise ValueError('The matrix must be a squre matrix')

        l = int(shape[0] * (shape[0] - 1) / 2)

        shape = (l,)

        obj = np.ndarray.__new__(subtype, shape, dtype, buffer, offset, strides, order)

        # Array initialization
        obj = obj * 0.0

        return obj

    def __getitem__(self, *kargs):
        """
        This function retrieves the selected elements i,j from the symmetric
        matrix A[i,j]
        
        Parameters
        ----------
        *kargs : python tuples
           the passed elements i,j

        Returns
        -------
            : float
            the selected element extracted from the allocated linear array
            
        """

        if isinstance(kargs[0], int):
            k = kargs[0]
            return super(SMatrix, self).__getitem__(k)

        if isinstance(kargs[0], slice):
            k = kargs[0]
            return super(SMatrix, self).__getitem__(k)

        elif len(kargs[0]) > 2:
            raise ValueError('Two indices can be addressed')

        i = kargs[0][0]
        j = kargs[0][1]

        if i == j:
            return 0.0

        # Length of the linear array 
        l = self.size

        # Total number of elements in the corresponding bi-dimensional symmetric matrix
        n = int((1 + math.sqrt(1 + 8 * l)) / 2)

        if i > n - 1:
            raise ValueError('First index out of bound')

        if j > n - 1:
            raise ValueError('Second index out of bound')

        if i < j:
            k = int((n * (n - 1) / 2) - (n - i) * ((n - i) - 1) / 2 + j - i - 1)
        else:
            k = int((n * (n - 1) / 2) - (n - j) * ((n - j) - 1) / 2 + i - j - 1)

        return super(SMatrix, self).__getitem__(k)

    def __setitem__(self, *kargs):
        """
        This function set the matrix elements i,j to the passed value
        
        Parameters
        ----------
        *kargs : python tuples
           the passed elements i,j, value to set  

        """

        if isinstance(kargs[0], int):
            k = kargs[0]
            value = kargs[1]
            return super(SMatrix, self).__setitem__(k, value)

        elif isinstance(kargs[0], slice):
            start, stop, step = kargs[0].indices(len(self))
            value = kargs[1]
            return super(SMatrix, self).__setitem__(kargs[0], value)

        elif len(kargs[0]) > 2:
            raise ValueError('Two indices can be addressed')

        # Passed indexes and value to set
        i = kargs[0][0]
        j = kargs[0][1]
        value = kargs[1]

        # Length of the linear array 
        l = self.size

        # Total number of elements in the corresponding bi-dimensional symmetric matrix
        n = int((1 + math.sqrt(1 + 8 * l)) / 2)

        if i > n - 1:
            raise ValueError('First index out of bound')
        if j > n - 1:
            raise ValueError('Second index out of bound')

        if i < j:
            k = int((n * (n - 1) / 2) - (n - i) * ((n - i) - 1) / 2 + j - i - 1)
        else:
            k = int((n * (n - 1) / 2) - (n - j) * ((n - j) - 1) / 2 + i - j - 1)
        super(SMatrix, self).__setitem__(k, value)

    def to_numpy_2D_array(self):
        """
        This function returns the symmetric similarity score numpy matrix 
        generated from the linear array
        
        Returns
        -------
        np_mat : numpy matrix
           the symmetric similarity score numpy matrix built by using the linear
           array 
        
        """

        # Length of the linear array 
        l = self.size

        # Total number of elements in the corresponding bi-dimensional symmetric matrix
        n = int((1 + math.sqrt(1 + 8 * l)) / 2)

        np_mat = np.zeros((n, n))

        for i in range(0, n):
            for j in range(0, n):
                np_mat[i, j] = self[i, j]

        return np_mat

    def mat_size(self):
        """
        This function returns the size of the square similarity score matrix 
        
        Returns
        -------
        n : int
           the size of the similarity score matrix
        
        """

        # Length of the linear array 
        l = self.size

        # Total number of elements in the corresponding bi-dimensional symmetric matrix
        n = int((1 + math.sqrt(1 + 8 * l)) / 2)

        return n

    # *************************


# Molecule Class
# *************************


class Molecule(object):
    """
    This Class stores the Rdkit molecule objects, their identification number 
    and the total number of instantiated molecules 

    """

    # This variable is used to count the current total number of molecules
    # The variable is defined as private
    __total_molecules = 0

    def __init__(self, molecule, mol_id, molname):
        """
        Initialization class function 
        
        Parameters
        ----------
        molecule : Rdkit molecule object
           the molecule

        mol_id : int
           the molecule identification number

        molname : str
           the molecule file name
        
        """

        # Check Inputs
        if not isinstance(molecule, Chem.rdchem.Mol):
            raise ValueError('The passed molecule object is not a RdKit molecule')

        if not isinstance(molname, str):
            raise ValueError('The passed molecule name must be a string')

        # The variable __molecule saves the current RDkit molecule object
        # The variable is defined as private
        self.__molecule = molecule

        # The variable __ID saves the molecule identification number 
        # The variable is defined as private
        self.__ID = mol_id

        # The variable __name saves the molecule identification name 
        # The variable is defined as private
        self.__name = molname

        # The variable __active saves whether the molecule is a known active
        # The variable is defined as private
        self.__active = False


    def getID(self):
        """
        Get the molecule ID number

        Returns
        -------
           : int
           the molecule ID number

        """
        return self.__ID

    def getMolecule(self):
        """
        Get the Rdkit molecule object

        Returns
        -------
        mol_copy : Rdkit molecule object
           The copy of the RDkit molecule

        """
        mol_copy = Chem.Mol(self.__molecule)
        return mol_copy

    def getName(self):
        """
        Get the molecule file name

        Returns
        -------
           : str 
           the molecule string file name

        """

        return self.__name

    def isActive(self):
        """
        Get whether the molecule is active

        Returns
        -------
           : bool
           the molecule active status

        """

        return self.__active

    def setActive(self, active):
        """
        Set whether the molecule is active

        """

        self.__active=active

class CheckDir(argparse.Action):
    # Classes used to check some of the passed user options in the main function
    # Class used to check the input directory
    def __call__(self, parser, namespace, directory, option_string=None):
        if not os.path.isdir(directory):
            raise argparse.ArgumentTypeError('The directory name is not a valid path: %s' % directory)
        if os.access(directory, os.R_OK):
            setattr(namespace, self.dest, directory)
        else:
            raise argparse.ArgumentTypeError('The directory name is not readable: %s' % directory)


class CheckPos(argparse.Action):
    # Class used to check the parallel, time and max user options
    def __call__(self, parser, namespace, value, option_string=None):
        if value < 1:
            raise argparse.ArgumentTypeError('%s is not a positive integer number' % value)
        setattr(namespace, self.dest, value)


class CheckCutoff(argparse.Action):
    # Class used to check the cutoff user option
    def __call__(self, parser, namespace, value, option_string=None):
        if not isinstance(value, float) or value < 0.0:
            raise argparse.ArgumentTypeError('%s is not a positive real number' % value)
        setattr(namespace, self.dest, value)


class CheckEcrscore(argparse.Action):
    # Class used to check the handicap user option
    def __call__(self, parser, namespace, value, option_string=None):
        if not isinstance(value, float) or value < 0.0 or value > 1.0:
            raise argparse.ArgumentTypeError('%s is not a real number in the range [0.0, 1.0]' % value)
        setattr(namespace, self.dest, value)


def startup():
    # Options and arguments passed by the user
    ops = parser.parse_args()

    # Molecule DataBase initialized with the passed user options
    db_mol = DBMolecules(ops.directory, ops.parallel, ops.verbose, ops.time, ops.ecrscore, ops.threed, ops.max3d, 
                         ops.output, ops.name, ops.output_no_images, ops.output_no_graph, ops.display, 
                         ops.allow_tree, ops.max, ops.cutoff, ops.radial, ops.hub, ops.fast, ops.links_file, 
                         ops.known_actives_file, ops.max_dist_from_actives)
    # Similarity score linear array generation
    strict, loose = db_mol.build_matrices()

    # Get the 2D numpy matrices
    # strict.to_numpy_2D_array()
    # loose.to_numpy_2D_array()

    # Graph generation based on the similarity score matrix
    nx_graph = db_mol.build_graph()

    # print nx_graph.nodes(data=True)
    # print nx_graph.edges(data=True)


# Command line user interface
# ----------------------------------------------------------------
parser = argparse.ArgumentParser(description='Lead Optimization Mapper 2. A program to plan alchemical relative '
                                             'binding affinity calculations',
                                 prog='LOMAP v. %s' % get_versions()['version'])
parser.add_argument('directory', action=CheckDir, \
                    help='The mol2/sdf file directory')
parser.add_argument('-p', '--parallel', default=1, action=CheckPos, type=int, \
                    help='Set the parallel mode. If an integer number N is specified, N processes will be executed to '
                         'build the similarity matrices')
parser.add_argument('-v', '--verbose', default='info', type=str, \
                    choices=['off', 'info', 'pedantic'], help='verbose mode selection')

mcs_group = parser.add_argument_group('MCS setting')
mcs_group.add_argument('-t', '--time', default=20, action=CheckPos, type=int, \
                       help='Set the maximum time in seconds to perform the mcs search between pair of molecules')
mcs_group.add_argument('-e', '--ecrscore', default=0.0, action=CheckEcrscore, type=float, \
                       help='If different from 0.0 the value is use to set the electrostatic score between two molecules with different charges')
mcs_group.add_argument('-3', '--threed', default=False, action='store_true', \
                       help='Use the input 3D coordinates to guide the preferred MCS mappings')
mcs_group.add_argument('-x', '--max3d', default=1000, type=float, \
                       help='The MCS is trimmed to remove atoms which are further apart than this distance')

out_group = parser.add_argument_group('Output setting')
out_group.add_argument('-o', '--output', default=True, action='store_true', \
                       help='Generates output files')
out_group.add_argument('-n', '--name', type=str, default='out', \
                       help='File name prefix used to generate the output files')
out_group.add_argument('--output-no-images', default=False, action='store_true', \
                       help='Disable the generation on the image files, removed the dependency on Pillow')
out_group.add_argument('--output-no-graph', default=False, action='store_true', \
                       help='Disable the generation on the graph (.dot) file, removed the dependency on pygraphviz')

parser.add_argument('-d', '--display', default=False, action='store_true', \
                    help='Display the generated graph by using Matplotlib')

graph_group = parser.add_argument_group('Graph setting')
graph_group.add_argument('-T', '--allow-tree', default=False, action='store_true', \
                        help='Remove the requirement that all molecules be in a cycle, so that the returned '
                             'graph will be a tree instead.');
graph_group.add_argument('-m', '--max', default=6, action=CheckPos, type=int, \
                         help='The maximum diameter of the graph')
graph_group.add_argument('-A', '--max-dist-from-actives', default=2, action=CheckPos, type=int, \
                         help='The maximum distance of any molecule from an active (requires -k)')
graph_group.add_argument('-c', '--cutoff', default=0.4, action=CheckCutoff, type=float, \
                         help='The Minimum Similarity Score (MSS) used to build the graph')
graph_group.add_argument('-r', '--radial', default=False, action='store_true', \
                         help='Using the radial option to build the graph')
graph_group.add_argument('-b', '--hub', default=None, type=str, \
                         help='Using a radial graph approach with a manually specified hub compound')
graph_group.add_argument('-a', '--fast', default=False, action='store_true', \
                         help='Using the fast graphing when the lead compound is specified')
graph_group.add_argument('-l', '--links-file', type=str, default='', \
                          help='Specify a filename listing the pairs of molecule files that should be initialised as linked.'
                          'Each line can be "mol1 mol2", which indicates that Lomap should compute the score and mapping '
                          'but that this link must be used in the final graph, or "mol1 mol2 score", which indicates that '
                          'Lomap should use the provided score, or "mol1 mol2 score force", which indicates that Lomap '
                          'should use the provided score and force this link to be used in the final graph.')
graph_group.add_argument('-k', '--known-actives-file', type=str, default='', \
                          help='Specify a filename listing the molecule files that should be initialised as "known actives", one per line')

# ------------------------------------------------------------------


# Main function
if "__main__" == __name__:
    startup()
