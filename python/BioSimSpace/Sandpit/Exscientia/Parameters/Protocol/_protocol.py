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
Functionality for handling parameterisation protocols.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Protocol"]

import os as _os
import queue as _queue
import shlex as _shlex
import subprocess as _subprocess

from ..._Utils import _try_import, _have_imported

# Temporarily redirect stderr to suppress import warnings.
import sys as _sys
_sys.stderr = open(_os.devnull, "w")

_openff = _try_import("openff")

if _have_imported(_openff):
    from openff.toolkit.topology import Molecule as _OpenFFMolecule
else:
    _OpenFFMolecule = _openff

# Reset stderr.
_sys.stderr = _sys.__stderr__
del _sys

from sire.legacy import IO as _SireIO
from sire.legacy import Mol as _SireMol
from sire.legacy import System as _SireSystem

from ... import _amber_home, _gmx_exe, _isVerbose
from ... import IO as _IO
from ... import _Utils
from ..._Exceptions import IncompatibleError as _IncompatibleError
from ..._Exceptions import MissingSoftwareError as _MissingSoftwareError
from ..._Exceptions import ParameterisationError as _ParameterisationError
from ..._Exceptions import ThirdPartyError as _ThirdPartyError
from ..._SireWrappers import Molecule as _Molecule

# Set the tLEaP cmd directory.
if _amber_home is not None:
    _cmd_dir = _amber_home + "/dat/leap/cmd"
    if not _os.path.isdir(_cmd_dir):
        raise IOError("Missing tLEaP dat directory: '%s'" % _cmd_dir)

# Search for all of the required executables.

# Search for the tLEaP exe.

_tleap_exe = None

if _amber_home is not None:
    _tleap_exe = "%s/bin/tleap" % _amber_home
    if not _os.path.isfile(_tleap_exe):
        raise IOError("Missing tLEaP executable: '%s'" % _tleap_exe)

# Search for the Antechamber exe.

_antechamber_exe = None

if _amber_home is not None:
    _antechamber_exe = "%s/bin/antechamber" % _amber_home
    if not _os.path.isfile(_antechamber_exe):
        raise IOError("Missing Antechamber executable: '%s'" % _antechamber_exe)

# Search for the parmchk exe.

_parmchk_exe = None

if _amber_home is not None:
    _parmchk_exe = "%s/bin/parmchk2" % _amber_home
    if not _os.path.isfile(_parmchk_exe):
        raise IOError("Missing parmchk executable: '%s'" % _parmchk_exe)

class Protocol():
    """A base class for parameterisation protocols."""

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

        # Don't allow user to create an instance of this base class.
        if type(self) is Protocol:
            raise Exception("<Protocol> must be subclassed.")

        # Validate and set the force field name.
        if not isinstance(forcefield, str):
            raise TypeError("'forcefield' must be of type 'str'")
        else:
            self._forcefield = forcefield

        # Validate and set the property map.
        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")
        else:
            self._property_map = property_map.copy()

        # Set default compatibility flags for the utility programs.
        # These must be overridden in the constructor of any derived classes.
        self._tleap = False
        self._pdb2gmx = False

        # Default to no water model for ion parameterisation.
        self._water_model = None

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

        if not isinstance(molecule, (_Molecule, str)):
            raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule' or 'str'")

        if work_dir is not None and not isinstance(work_dir, str):
            raise TypeError("'work_dir' must be of type 'str'")

        if queue is not None and not isinstance(queue, _queue.Queue):
            raise TypeError("'queue' must be of type 'queue.Queue'")

        # Set work_dir to the current directory.
        if work_dir is None:
            work_dir = _os.getcwd()

        # Flag whether the molecule is a SMILES string.
        if isinstance(molecule, str):
            is_smiles = True
        else:
            is_smiles = False

        # Create the file prefix.
        prefix = work_dir + "/"

        if not is_smiles:
            # Create a copy of the molecule.
            new_mol = molecule.copy()

        # Choose the program to run with depending on the force field compatibility.
        # If tLEaP and pdb2gmx are supported, default to tLEaP, then use pdb2gmx if
        # tLEaP fails to produce output.

        # First, try parameterise using tLEaP.
        if self._tleap:
            if _tleap_exe is not None:
                output = self._run_tleap(molecule, work_dir)
                if not is_smiles:
                    new_mol._ion_water_model = self._water_model
            # Otherwise, try using pdb2gmx.
            elif self._pdb2gmx:
                if _gmx_exe is not None:
                    output = self._run_pdb2gmx(molecule, work_dir)
                else:
                    raise _MissingSoftwareError("Cannot parameterise. Missing AmberTools and GROMACS.")

        # Parameterise using pdb2gmx.
        elif self._pdb2gmx:
            if _gmx_exe is not None:
                output = self._run_pdb2gmx(molecule, work_dir)
            else:
                raise _MissingSoftwareError("Cannot use pdb2gmx since GROMACS is not installed!")

        # Prepend the working directory to the output file names.
        output = [prefix + output[0], prefix + output[1]]

        try:
            # Load the parameterised molecule. (This could be a system of molecules.)
            par_mol = _IO.readMolecules(output)
            # Extract single molecules.
            if par_mol.nMolecules() == 1:
                par_mol = par_mol.getMolecules()[0]
        except Exception as e:
            msg = "Failed to read molecule from: '%s', '%s'" % (output[0], output[1])
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # Make the molecule 'mol' compatible with 'par_mol'. This will create
        # a mapping between atom indices in the two molecules and add all of
        # the new properties from 'par_mol' to 'mol'.
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
            new_mol.makeCompatibleWith(par_mol, property_map=self._property_map, overwrite=True, verbose=False)

        # Record the forcefield used to parameterise the molecule.
        new_mol._forcefield = self._forcefield

        if queue is not None:
            queue.put(new_mol)
        return new_mol

    def _has_missing_atoms(self, tleap_file):
        """Check whether tLEaP has added missing atoms.

           Parameters
           ----------

           tleap_file : str
               The path to the output file generated by tLEaP.

           Returns
           -------

           has_missing_atoms : bool
               Whether missing atoms were added.
        """

        if not isinstance(tleap_file, str):
            raise TypeError("'tleap_file' must be of type 'str'.")

        if not (_os.path.exists(tleap_file)):
            raise IOError(f"'tleap_file' doesn't exist: '{tleap_file}'")

        with open(tleap_file, "r") as file:
            for line in file:
                if "missing atoms" in line:
                    return True

        return False

    def _smiles_to_molecule(self, smiles, work_dir):
        """Convert a SMILES string to a Molecule.

           Parameters
           ----------

           smiles : str
               A SMILES representation of a molecule.

           work_dir : str
               The working directory.

           Returns
           -------

           molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
               The BioSimSpace representation of the molecule.
        """

        # Create an OpenFF molecule from the smiles string.
        off_molecule = _OpenFFMolecule.from_smiles(smiles)

        # Generate a single conformer for the molecule.
        off_molecule.generate_conformers(n_conformers=1)

        # Write to file.
        file = work_dir + "/off.pdb"
        off_molecule.to_file(file, "pdb")

        # Read PDB back into internal molecular representation.
        molecule = _IO.readMolecules(file)

        # Assume this is a single molecule.
        return molecule[0]

    def _run_tleap(self, molecule, work_dir):
        """Run using tLEaP.

           Parameters
           ----------

           molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
               The molecule to parameterise, either as a Molecule object or SMILES
               string.

           work_dir : str
               The working directory.
        """

        # Convert SMILES to a molecule.
        if isinstance(molecule, str):
            try:
                _molecule = self._smiles_to_molecule(molecule, work_dir)
            except Exception as e:
                msg = "Unable to convert SMILES to Molecule using Open Force Field."
                if _isVerbose():
                    msg += ": " + getattr(e, "message", repr(e))
                    raise _ThirdPartyError(msg) from e
                else:
                    raise _ThirdPartyError(msg) from None
        else:
            _molecule = molecule

        # Create a new system and molecule group.
        s = _SireSystem.System("BioSimSpace System")
        m = _SireMol.MoleculeGroup("all")

        # Add the molecule.
        m.add(_molecule._getSireObject())
        s.add(m)

        # Create the file prefix.
        prefix = work_dir + "/"

        # Write the system to a PDB file.
        try:
            pdb = _SireIO.PDB2(s, self._property_map)
            pdb.writeToFile(prefix + "leap.pdb")
        except Exception as e:
            msg = "Failed to write system to 'PDB' format."
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # Try to find a force field file.
        ff = _find_force_field(self._forcefield)

        # Write the LEaP input file.
        with open(prefix + "leap.txt", "w") as file:
            file.write("source %s\n" % ff)
            if self._water_model is not None:
                if self._water_model in ["tip4p", "tip5p"]:
                    file.write("source leaprc.water.tip4pew\n")
                else:
                    file.write("source leaprc.water.%s\n" % self._water_model)
            # Write extra user commands.
            if self._leap_commands is not None:
                for command in self._leap_commands:
                    file.write("%s\n" % command)
            file.write("mol = loadPdb leap.pdb\n")
            file.write("saveAmberParm mol leap.top leap.crd\n")
            file.write("quit")

        # Generate the tLEaP command.
        command = "%s -f leap.txt" % _tleap_exe

        with open(prefix + "README.txt", "w") as file:
            # Write the command to file.
            file.write("# tLEaP was run with the following command:\n")
            file.write("%s\n" % command)

        # Create files for stdout/stderr.
        stdout = open(prefix + "leap.out", "w")
        stderr = open(prefix + "leap.err", "w")

        # Run tLEaP as a subprocess.
        proc = _subprocess.run(_Utils.command_split(command), cwd=work_dir,
            shell=False, stdout=stdout, stderr=stderr)
        stdout.close()
        stderr.close()

        # tLEaP doesn't return sensible error codes, so we need to check that
        # the expected output was generated.
        if _os.path.isfile(prefix + "leap.top") and _os.path.isfile(prefix + "leap.crd"):
            # Check the output of tLEaP for missing atoms.
            if self._has_missing_atoms(prefix + "leap.out"):
                raise _ParameterisationError("tLEaP added missing atoms. The topology is now "
                                             "inconsistent with the original molecule. Please "
                                             "make sure that your initial molecule has a "
                                             "complete topology.")
            return ["leap.top", "leap.crd"]
        else:
            raise _ParameterisationError("tLEaP failed!")

    def _run_pdb2gmx(self, molecule, work_dir):
        """Run using pdb2gmx.

           Parameters
           ----------

           molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
               The molecule to parameterise, either as a Molecule object or SMILES
               string.

           work_dir : str
               The working directory.
        """

        # A list of supported force fields, mapping to their GROMACS ID string.
        # GROMACS supports a sub-set of the AMBER force fields.
        supported_ff = { "ff99"   : "amber99",
                         "ff99SB" : "amber99sb",
                         "ff03"   : "amber03"
                       }

        if self._forcefield not in supported_ff:
            raise _IncompatibleError("'pdb2gmx' does not support the '%s' force field." % self._forcefield)

        # Convert SMILES to a molecule.
        if isinstance(molecule, str):
            try:
                _molecule = self._smiles_to_molecule(molecule, work_dir)
            except Exception as e:
                msg = "Unable to convert SMILES to Molecule using Open Force Field."
                if _isVerbose():
                    msg += ": " + getattr(e, "message", repr(e))
                    raise _ThirdPartyError(msg) from e
                else:
                    raise _ThirdPartyError(msg) from None
        else:
            _molecule = molecule

        # Create a new system and molecule group.
        s = _SireSystem.System("BioSimSpace System")
        m = _SireMol.MoleculeGroup("all")

        # Add the molecule.
        m.add(_molecule._getSireObject())
        s.add(m)

        # Create the file prefix.
        prefix = work_dir + "/"

        # Write the system to a PDB file.
        try:
            pdb = _SireIO.PDB2(s, self._property_map)
            pdb.writeToFile(prefix + "input.pdb")
        except Exception as e:
            msg = "Failed to write system to 'PDB' format."
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # Generate the pdb2gmx command.
        command = "%s pdb2gmx -f input.pdb -o output.gro -p output.top -ignh -ff %s -water none" \
            % (_gmx_exe, supported_ff[self._forcefield])

        with open(prefix + "README.txt", "w") as file:
            # Write the command to file.
            file.write("# pdb2gmx was run with the following command:\n")
            file.write("%s\n" % command)

        # Create files for stdout/stderr.
        stdout = open(prefix + "pdb2gmx.out", "w")
        stderr = open(prefix + "pdb2gmx.err", "w")

        # Run pdb2gmx as a subprocess.
        proc = _subprocess.run(_Utils.command_split(command), cwd=work_dir,
            shell=False, stdout=stdout, stderr=stderr)
        stdout.close()
        stderr.close()

        # Check for the expected output.
        if _os.path.isfile(prefix + "output.gro") and _os.path.isfile(prefix + "output.top"):
            return ["output.gro", "output.top"]
        else:
            raise _ParameterisationError("pdb2gmx failed!")

def _find_force_field(forcefield):
    """Internal function to search LEaP compatible force field files.

       Parameters
       ----------

       forcefield : str
           The name of the force field.

       Returns
       -------

       file : str
           The full path of the matching force field file.
    """

    # Whether the force field is old.
    is_old = False

    # Search for a compatible force field file.
    ff = _IO.glob("%s/*.%s" % (_cmd_dir, forcefield))

    # Search the old force fields. First try a specific match.
    if len(ff) == 0:
        ff = _IO.glob("%s/oldff/leaprc.%s" % (_cmd_dir, forcefield))
        is_old = True

        # No matches, try globbing all files with matching extension.
        if len(ff) == 0:
            ff = _IO.glob("%s/oldff/*.%s" % (_cmd_dir, forcefield))

    # No force field found!
    if len(ff) == 0:
        raise ValueError("No force field file found for '%s'" % forcefield)

    # Multiple force fields found.
    elif len(ff) > 1:
        raise ValueError("Multiple force fields found for '%s': %s" % (forcefield, ff))

    # Create the force field name.
    ff = _os.path.basename(ff[0])
    if is_old:
        ff = "oldff/" + ff

    # Return the force field.
    return ff
