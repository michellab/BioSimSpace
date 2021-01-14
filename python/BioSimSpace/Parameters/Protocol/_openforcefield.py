######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2021
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# BioSimSpace is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BioSimSpace. If not, see <http://www.gnu.org/licenses/>.
#####################################################################

"""
Functionality for handling parameterisation protocols
for force fields from the Open Force Field Initiative.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["OpenForceField"]

# To override any protocols, just implement a custom "run" method in any
# of the classes.

import os as _os
import parmed as _parmed
import queue as _queue
import subprocess as _subprocess

import warnings as _warnings
# Suppress numpy warnings from RDKit import.
_warnings.filterwarnings("ignore", message="numpy.dtype size changed")
_warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
# Suppress duplicate to-Python converted warnings.
# Both Sire and RDKit register the same converter.
with _warnings.catch_warnings():
    _warnings.filterwarnings("ignore")
    from rdkit import Chem as _Chem
    from rdkit import RDLogger as _RDLogger

    # Disable RDKit warnings.
    _RDLogger.DisableLog('rdApp.*')

import sys as _sys
# Temporarily redirect stderr to suppress import warnings.
_sys.stderr = open(_os.devnull, "w")
from simtk.openmm.app import PDBFile as _PDBFile
from openforcefield.topology import Molecule as _OpenFFMolecule
from openforcefield.topology import Topology as _OpenFFTopology
from openforcefield.typing.engines.smirnoff import ForceField as _Forcefield
# Reset stderr.
_sys.stderr = _sys.__stderr__
del _sys

from Sire import IO as _SireIO
from Sire import Mol as _SireMol
from Sire import System as _SireSystem

from BioSimSpace import _isVerbose
from BioSimSpace import IO as _IO
from BioSimSpace._SireWrappers import Molecule as _Molecule

from . import _protocol

class OpenForceField(_protocol.Protocol):
    """A class for handling protocols for Open Force Field models."""

    def __init__(self, forcefield, property_map={}):
        """Constructor.

           Parameters
           ----------

           forcefield : str
               The name of the force field.

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(forcefield=forcefield, property_map=property_map)

        # Set the compatibility flags.
        self._tleap = False
        self._pdb2gmx = False

    def run(self, molecule, work_dir=None, queue=None):
        """Run the parameterisation protocol.

           Parameters
           ----------

           molecule : BioSimSpace._SireWrappers.Molecule
               The molecule to apply the parameterisation protocol to.

           work_dir : str
               The working directory.

           queue : queue.Queue
               The thread queue is which this method has been run.

           Returns
           -------

           molecule : BioSimSpace._SireWrappers.Molecule
               The parameterised molecule.
        """

        if type(molecule) is not _Molecule:
            raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'")

        if type(work_dir) is not None and type(work_dir) is not str:
            raise TypeError("'work_dir' must be of type 'str'")

        if type(queue) is not None and type(queue) is not _queue.Queue:
            raise TypeError("'queue' must be of type 'queue.Queue'")

        # Set work_dir to the current directory.
        if work_dir is None:
            work_dir = _os.getcwd()

        # Create the file prefix.
        prefix = work_dir + "/"

        # Create a copy of the molecule.
        new_mol = molecule.copy()

        # Write the molecule to a PDB file.
        try:
            pdb = _SireIO.PDB2(new_mol.toSystem()._sire_object)
            pdb.writeToFile(prefix + "molecule.pdb")
        except Exception as e:
            msg = "Failed to write the molecule to 'PDB' format."
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # Create an RDKit molecule from the PDB file.
        try:
            rdmol = _Chem.MolFromPDBFile(prefix + "molecule.pdb", removeHs=False)
        except Exception as e:
            raise IOError("RDKit was unable to read the molecular PDB file!") from None

        # Use RDKit to write back to PDB format so that we generate a CONECT
        # record, which is required by OpenFF.
        try:
            _Chem.MolToPDBFile(rdmol, prefix + "molecule.pdb")
        except Exception as e:
            raise IOError("RDKit was unable to write the molecular PDB file!") from None

        # Obtain the OpenMM Topology object from the PDB file.
        try:
            pdbfile = _PDBFile(prefix + "molecule.pdb")
            omm_topology = pdbfile.topology
        except Exception as e:
            raise IOError("OpenMM was unable to read the molecular PDB file!") from None

        # Use RDKit to generate the smiles string for the molecule.
        try:
            # Split on a dot in case multiple molecules were present.
            smiles = _Chem.MolToSmiles(rdmol).split(".")
        except Exception as e:
            raise IOError("RDKit was unable to generate a SMILE string for the molecule!") from None

        # Convert the SMILES string to an OpenFF molecule.
        try:
            off_molecules = []
            for s in smiles:
                off_molecules.append(_OpenFFMolecule.from_smiles(s))
        except Exception as e:
            raise IOError("Unable to convert SMILES to an OpenFF Molecule!") from None

        # Create the Open Forcefield Topology.
        try:
            off_topology = _OpenFFTopology.from_openmm(omm_topology, unique_molecules=off_molecules)
        except Exception as e:
            raise IOError("Unable to create OpenFF Topology!") from None

        # Load the force field.
        try:
            ff = self._forcefield + ".offxml"
            forcefield = _Forcefield(ff)
        except Exception as e:
            raise IOError(f"Unable to load force field: {ff}") from None

        # Create an OpenMM system.
        try:
            omm_system = forcefield.create_openmm_system(off_topology)
        except Exception as e:
            raise IOError("Unable to create OpenMM System!") from None

        # Convert the OpenMM System to a ParmEd structure.
        try:
            parmed_structure = _parmed.openmm.load_topology(omm_topology, omm_system, pdbfile.positions)
        except Exception as e:
            raise IOError("Unable to convert OpenMM System to ParmEd structure!") from None

        # Export AMBER format files.
        try:
            parmed_structure.save(prefix + "parmed.prmtop", overwrite=True)
            parmed_structure.save(prefix + "parmed.inpcrd", overwrite=True)
        except Exception as e:
            raise IOError("Unable to write ParmEd structure to AMBER format!!") from None

        # Load the parameterised molecule. (This could be a system of molecules.)
        try:
            par_mol = _IO.readMolecules([prefix + "parmed.prmtop", prefix + "parmed.inpcrd"])
            # Extract single molecules.
            if par_mol.nMolecules() == 1:
                par_mol = par_mol.getMolecules()[0]
        except Exception as e:
            msg = "Failed to read molecule from: 'parmed.prmtop', 'parmed.inpcrd'"
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # Make the parameterised molecule compatible with the original topology.
        new_mol.makeCompatibleWith(par_mol, property_map=self._property_map, overwrite=True, verbose=False)

        # Record the forcefield used to parameterise the molecule.
        new_mol._forcefield = self._forcefield

        if queue is not None:
            queue.put(new_mol)
        return new_mol
