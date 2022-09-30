######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2022
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
__email__ = "lester.hedges@gmail.com"

__all__ = ["OpenForceField"]

# To override any protocols, just implement a custom "run" method in any
# of the classes.

from ..._Utils import _try_import, _have_imported

import os as _os
_parmed = _try_import("parmed")
import queue as _queue
import subprocess as _subprocess

import warnings as _warnings
# Suppress numpy warnings from RDKit import.
_warnings.filterwarnings("ignore", message="numpy.dtype size changed")
_warnings.filterwarnings("ignore", message="numpy.ndarray size changed")
_warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
# Suppress duplicate to-Python converted warnings.
# Both Sire and RDKit register the same converter.
with _warnings.catch_warnings():
    _warnings.filterwarnings("ignore")
    _rdkit = _try_import("rdkit")

    if _have_imported(_rdkit):
        from rdkit import Chem as _Chem
        from rdkit import RDLogger as _RDLogger

        # Disable RDKit warnings.
        _RDLogger.DisableLog('rdApp.*')
    else:
        _Chem = _rdkit
        _RDLogger = _rdkit

import sys as _sys
# Temporarily redirect stderr to suppress import warnings.
_sys.stderr = open(_os.devnull, "w")

_openmm = _try_import("openmm")

if _have_imported(_openmm):
    from openmm.app import PDBFile as _PDBFile
else:
    _PDBFile = _openmm

_openff = _try_import("openff")

if _have_imported(_openff):
    from openff.toolkit.topology import Molecule as _OpenFFMolecule
    from openff.toolkit.topology import Topology as _OpenFFTopology
    from openff.toolkit.typing.engines.smirnoff import ForceField \
        as _Forcefield
else:
    _OpenFFMolecule = _openff
    _OpenFFTopology = _openff
    _Forcefield = _openff

# Reset stderr.
_sys.stderr = _sys.__stderr__
del _sys

from sire.legacy import IO as _SireIO
from sire.legacy import Mol as _SireMol
from sire.legacy import System as _SireSystem

from ... import _isVerbose
from ... import IO as _IO
from ..._Exceptions import IncompatibleError as _IncompatibleError
from ..._Exceptions import ThirdPartyError as _ThirdPartyError
from ..._SireWrappers import Molecule as _Molecule

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

           molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
               The molecule to parameterise, either as a Molecule object or SMILES
               string.

           work_dir : str
               The working directory.

           queue : queue.Queue
               The thread queue is which this method has been run.

           Returns
           -------

           molecule : BioSimSpace._SireWrappers.Molecule
               The parameterised molecule.
        """

        if not isinstance(molecule, (_Molecule, str)):
            raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule' or 'str'")

        if work_dir is not None and not isinstance(work_dir, str):
            raise TypeError("'work_dir' must be of type 'str'")

        if queue is not None and not isinstance(queue, _queue.Queue):
            raise TypeError("'queue' must be of type 'queue.Queue'")

        # Set work_dir to the current directory.
        if work_dir is None:
            work_dir = _os.getcwd()

        # Create the file prefix.
        prefix = work_dir + "/"

        # Flag whether the molecule is a SMILES string.
        if isinstance(molecule, str):
            is_smiles = True
        else:
            is_smiles = False

        # The following is adapted from the Open Force Field examples, where an
        # OpenFF system is converted to AMBER format files using ParmEd:
        # https://github.com/openforcefield/openff-toolkit/blob/master/examples/using_smirnoff_in_amber_or_gromacs/convert_to_amber_gromacs.ipynb

        if is_smiles:
            # Convert SMILES string to an OpenFF molecule.
            try:
                off_molecule = _OpenFFMolecule.from_smiles(molecule)
            except Exception as e:
                msg = "Failed to convert SMILES to Open Force Field Molecule."
                if _isVerbose():
                    msg += ": " + getattr(e, "message", repr(e))
                    raise IOError(msg) from e
                else:
                    raise IOError(msg) from None

            # Generate a single conformer.
            try:
                off_molecule.generate_conformers(n_conformers=1)
            except Exception as e:
                msg = "Unable to generate conformer from Open Force Field molecule."
                if _isVerbose():
                    msg += ": " + getattr(e, "message", repr(e))
                    raise IOError(msg) from e
                else:
                    raise IOError(msg) from None
        else:
            # Write the molecule to a PDB file.
            try:
                _IO.saveMolecules(prefix + "molecule", molecule, "pdb")
            except Exception as e:
                msg = "Failed to write the molecule to 'PDB' format."
                if _isVerbose():
                    msg += ": " + getattr(e, "message", repr(e))
                    raise IOError(msg) from e
                else:
                    raise IOError(msg) from None

            # Create an RDKit molecule from the PDB file.
            try:
                rdmol = _Chem.MolFromPDBFile(prefix + "molecule.pdb", removeHs=False)
            except Exception as e:
                msg = "RDKit was unable to read the molecular PDB file!"
                if _isVerbose():
                    msg += ": " + getattr(e, "message", repr(e))
                    raise _ThirdPartyError(msg) from e
                else:
                    raise _ThirdPartyError(msg) from None

            # Use RDKit to write back to SDF format.
            try:
                writer = _Chem.SDWriter(prefix + "molecule.sdf")
                writer.write(rdmol)
                writer.close()
            except Exception as e:
                msg = "RDKit was unable to write the molecular SDF file!"
                if _isVerbose():
                    msg += ": " + getattr(e, "message", repr(e))
                    raise _ThirdPartyError(msg) from e
                else:
                    raise _ThirdPartyError(msg) from None

            # Create the Open Forcefield Molecule from the intermediate SDF file,
            # as recommended by @j-wags and @mattwthompson.
            try:
                off_molecule = _OpenFFMolecule.from_file(prefix + "molecule.sdf")
            except Exception as e:
                msg = "Unable to create OpenFF Molecule!"
                if _isVerbose():
                    msg += ": " + getattr(e, "message", repr(e))
                    raise _ThirdPartyError(msg) from e
                else:
                    raise _ThirdPartyError(msg) from None

        # Extract the molecular topology.
        try:
            off_topology = off_molecule.to_topology()
        except Exception as e:
            msg = "Unable to create OpenFF Topology!"
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise _ThirdPartyError(msg) from e
            else:
                raise _ThirdPartyError(msg) from None

        # Load the force field.
        try:
            ff = self._forcefield + ".offxml"
            forcefield = _Forcefield(ff)
        except Exception as e:
            msg = f"Unable to load force field: {ff}"
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise _ThirdPartyError(msg) from e
            else:
                raise _ThirdPartyError(msg) from None

        # Obtain the OpenMM Topology object from the OpenFF topology.
        try:
            omm_topology = off_topology.to_openmm()
        except Exception as e:
            msg = "Unable to convert Open Force Field topology to OpenMM topology!"
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise _ThirdPartyError(msg) from e
            else:
                raise _ThirdPartyError(msg) from None

        # Create an OpenMM system.
        try:
            omm_system = forcefield.create_openmm_system(off_topology)
        except Exception as e:
            msg = "Unable to create OpenMM System!"
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise _ThirdPartyError(msg) from e
            else:
                raise _ThirdPartyError(msg) from None

        # Convert the OpenMM System to a ParmEd structure.
        try:
            if is_smiles:
                parmed_structure = _parmed.openmm.load_topology(omm_topology, omm_system, off_molecule.conformers[0])
            else:
                pdbfile = _PDBFile(prefix + "molecule.pdb")
                parmed_structure = _parmed.openmm.load_topology(omm_topology, omm_system, pdbfile.positions)
        except Exception as e:
            msg = "Unable to convert OpenMM System to ParmEd structure!"
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise _ThirdPartyError(msg) from e
            else:
                raise _ThirdPartyError(msg) from None

        # Export AMBER format files.
        try:
            parmed_structure.save(prefix + "parmed.prmtop", overwrite=True)
            parmed_structure.save(prefix + "parmed.inpcrd", overwrite=True)
        except Exception as e:
            msg = "Unable to write ParmEd structure to AMBER format!"
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise _ThirdPartyError(msg) from e
            else:
                raise _ThirdPartyError(msg) from None

        # Load the parameterised molecule. (This could be a system of molecules.)
        try:
            par_mol = _IO.readMolecules([prefix + "parmed.prmtop", prefix + "parmed.inpcrd"])
            # Extract single molecules.
            if par_mol.nMolecules() == 1:
                par_mol = par_mol.getMolecules()[0]
        except Exception as e:
            msg = "Failed to read molecule from: 'parmed.prmtop', 'parmed.inpcrd'"
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # Make the parameterised molecule compatible with the original topology.
        if is_smiles:
            new_mol = par_mol

            # We'll now add MolName and ResName info to the molecule, since
            # this will be missing.

            # Rename the molecule with the original SMILES string.
            edit_mol = new_mol._sire_object.edit()
            edit_mol = edit_mol.rename(molecule).molecule()

            # Rename the residue LIG.
            resname = _SireMol.ResName("LIG")
            edit_mol = edit_mol.residue(_SireMol.ResIdx(0)).rename(resname).molecule()

            # Commit the changes.
            new_mol._sire_object = edit_mol.commit()

        else:
            new_mol = molecule.copy()
            new_mol.makeCompatibleWith(par_mol, property_map=self._property_map, overwrite=True, verbose=False)

        # Record the forcefield used to parameterise the molecule.
        new_mol._forcefield = self._forcefield

        if queue is not None:
            queue.put(new_mol)
        return new_mol
