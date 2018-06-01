######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2018
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

# TODO:
# This should ultimately be a BioSimSpace.Process object since we need state
# that persists beyond the lifetime of a function call. This allows the user
# to query stdout/stderr when things go wrong, and examine intermediate files.

"""
Functionality for parameterising molecules using the Antechamber package.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

import Sire as _Sire

from BioSimSpace import _amber_home
from ._tleap import Tleap as _Tleap

import BioSimSpace.IO as _IO

import os as _os
import subprocess as _subprocess
import tempfile as _tempfile

__all__ = ["Antechamber"]

# Search for the Antechamber exe.

_antechamber_exe = None

if _amber_home is not None:
    if _os.path.isfile("%s/bin/antechamber" % _amber_home):
        _antechamber_exe = "%s/bin/antechamber" % _amber_home

if _antechamber_exe is None:
    # Search Sire bin directory.
    _bin_dir = _Sire.Base.getBinDir()
    _antechamber_exe = "%s/antechamber" % _bin_dir

    # Search system PATH.
    if not _os.path.isfile(_antechamber_exe):
        _antechamber_exe = _Sire.Base.findExe("antechamber").absoluteFilePath()

# Search for the parmchk exe.

_parmchk_exe = None

if _amber_home is not None:
    if _os.path.isfile("%s/bin/parmchk2" % _amber_home):
        _parmchk_exe = "%s/bin/parmchk2" % _amber_home

if _parmchk_exe is None:
    # Search Sire bin directory.
    _bin_dir = _Sire.Base.getBinDir()
    _parmchk_exe = "%s/parmchk2" % _bin_dir

    # Search system PATH.
    if not _os.path.isfile(_parmchk_exe):
        _parmchk_exe = _Sire.Base.findExe("parmchk_exe").absoluteFilePath()

class Antechamber():
    """A simple wrapper around the Antechamber package from AmberTools."""

    @staticmethod
    def parameterise(molecule, forcefield, frcmod=None, work_dir=None, verbose=False):
        """Parameterise a molecule with the general AMBER force field (GAFF).

           Positional arguments:

           molecule   -- The molecule to parameterise.
           forcefield -- The forcefield to use (GAFF or GAFF2).

           Keyword arguments:

           work_dir   -- The working directory for external processes.
           verbose    -- Whether to report stdout/stderr of external processes.
        """

        if type(molecule) is not _Sire.Mol._Mol.Molecule:
            raise TypeError("'molecule' must be of type 'Sire.Mol._Mol.Molecule'")

        if type(forcefield) is not str:
            raise TypeError("'forcefield' must be of type 'str'")
        else:
            if forcefield != "gaff" and forcefield != "gaff2":
                raise ValueError("Supported force fields are 'gaff' and 'gaff2'")

        # Create a new system and molecule group.
        s = _Sire.System.System("BioSimSpace molecule")
        m = _Sire.Mol.MoleculeGroup("all")

        # Add the molecule.
        m.add(molecule)
        s.add(m)

        # Create a temporary working directory and store the directory name.
        if work_dir is None:
            tmp_dir = _tempfile.TemporaryDirectory()
            work_dir = self._tmp_dir.name

        # User specified working directory.
        else:
            # Create the directory if it doesn't already exist.
            if not _os.path.isdir(work_dir):
                _os.makedirs(work_dir)

        # Store the current working directory.
        dir = _os.getcwd()

        # Change to the working directory.
        _os.chdir(work_dir)

        # Write the system to a PDB file.
        pdb = _Sire.IO.PDB2(s)
        pdb.writeToFile("antechamber.pdb")

        # Generate the Antechamber command.
        command = "%s -i antechamber.pdb -fi pdb -o antechamber.mol2 -fo mol2 -c bcc -s 2" % _antechamber_exe

        # Run Antechamber as a subprocess.
        if verbose:
            proc = _subprocess.run(command, shell=True)
        else:
            proc = _subprocess.run(command, shell=True,
                stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)

        # Antechamber doesn't return sensible error codes, so we need to check that
        # the expected output was generated.
        if _os.path.isfile("antechamber.mol2"):

            # Run parmchk to check for missing parameters.
            command = "%s -i antechamber.mol2 -f mol2 -o antechamber.frcmod" % _parmchk_exe

            # Run parmchk as a subprocess.
            if verbose:
                proc = _subprocess.run(command, shell=True)
            else:
                proc = _subprocess.run(command, shell=True,
                    stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)

            # The frcmod file was created.
            if _os.path.isfile("antechamber.frcmod"):

                # Change back to the current directory.
                _os.chdir(dir)

                # Now call tLEaP using the partially parameterised molecule and the frcmod file.
                # tLEap will run in the same working directory, using the Mol2 file generated by
                # Antechamber. As such we can pass None for the molecule argument.
                return _Tleap.parameterise(None, forcefield, work_dir=work_dir,
                    verbose=verbose, is_ante=True)
            else:
                raise RuntimeError("parmchk failed to generate the expected output!")

        else:
            # Change back to the current directory.
            _os.chdir(dir)

            raise RuntimeError("Antechamber failed to generate the expected output!")
