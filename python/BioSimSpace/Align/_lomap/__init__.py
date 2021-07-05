"""
Lomap2
======

Alchemical free energy calculations hold increasing promise as an aid to drug 
discovery efforts. However, applications of these techniques in discovery 
projects have been relatively few, partly because of the difficulty of planning 
and setting up calculations. The Lead Optimization Mapper (LOMAP) is an 
automated algorithm to plan efficient relative free energy calculations between 
potential ligands within a substantial of compounds.

Authors: Gaetano Calabro' <gcalabro@uci.edu> 
         David Mobley     <dmobley@uci.edu>


Licence: LGPL

URL: https://github.com/nividic/Lomap


Using
-----
      Just write in Python 

      # Import Lomap
      import lomap

      # Generate the molecule database starting from 
      # a directory containing .mol2 files

      db_mol = lomap.DBMolecules("lomap/test/basic", output=True)

      # Calculate the similarity matrix betweeen the database 
      # molecules. Two molecules are generated related to the 
      # scrict rule and loose rule 

      strict, loose = db_mol.build_matrices()

      # Generate the NetworkX graph and output the results
      nx_graph = db_mol.build_graph() 


      # Calculate the Maximum Common Subgraph (MCS) between 
      # the first two molecules in the molecule database 
      # ignoring hydrogens and depicting the mapping in a file
    
      MC = lomap.MCS.getMapping(db_mol[0].getMolecule(), db_mol[1].getMolecule(), hydrogens=False, fname='mcs.png')

 
      # Alchemical transformation are usually performed between molecules with
      # the same charges. However, it is possible to allow these transformations
      # manually setting the electrostatic score for the whole set of molecules 
      # producing a connected graph. The electrostatic scrore must be in the 
      # range [0,1]


      db_mol = lomap.DBMolecules("python string pointing to a directory with mol2 files", output=True, ecrscore=0.1)
      strict, loose = db_mol.build_matrices()
      nx_graph = db_mol.build_graph() 
"""

from .dbmol import DBMolecules
from .dbmol import SMatrix
from .dbmol import Molecule
from .mcs import MCS

del dbmol
del mcs

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
