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

"""
Functionality for solvating molecular systems.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

import Sire as _Sire

from BioSimSpace import _bin_dir

from .._SireWrappers import System as _System
from .._SireWrappers import Molecule as _Molecule

import BioSimSpace.IO as _IO

import os as _os
import shutil as _shutil
import subprocess as _subprocess
import sys as _sys
import tempfile as _tempfile

# Search for the gmx exe.

_gmx_exe = "%s/gmx" % _bin_dir
if not _os.path.isfile(_gmx_exe):
    _gmx_exe = _Sire.Base.findExe("gmx").absoluteFilePath()

def spc(molecule=None, box=None):
    """Add solvent compatible with all three-point water models.

       Keyword arguments:

       molecule -- A molecule, or system of molecules.
       box      -- A list containing the box size in each dimension (in nm).
    """

    # Validate arguments.
    molecule, box = _validate_input(molecule, box)

    # Create the solvated system.
    return _solvate(molecule, box, "spc216")

def tip3p(molecule=None, box=None):
    """Add solvent compatible with all three-point water models.

       Keyword arguments:

       molecule -- A molecule, or system of molecules.
       box      -- A list containing the box size in each dimension (in nm).
    """

    # Validate arguments.
    molecule, box = _validate_input(molecule, box)

    # Create the solvated system.
    return _solvate(molecule, box, "spc216")

def tip4p(molecule=None, box=None):
    """Add TIP4P solvent.

       Keyword arguments:

       molecule -- A molecule, or system of molecules.
       box      -- A list containing the box size in each dimenion (in nm).
    """

    # Validate arguments.
    molecule, box = _validate_input(molecule, box)

    # Return the solvated system.
    return _solvate(molecule, box, "tip4p")

def _validate_input(molecule, box):
    """Internal function to validate function arguments.

       Positional arguments:

       molecule -- A molecule, or system of molecules.
       box      -- A list containing the box size in each dimension (in nm).
    """

    if molecule is not None:
        if type(molecule) is not _Molecule and type(molecule) is not _System:
            raise TypeError("'molecule' must be of type 'BioSimSpace.Molecule' "
                + "or 'BioSimSpace.System'")

        # Try to extract the box dimensions from the system.
        if type(molecule) is _System and box is None:
            try:
                box = system.property("space").dimensions()
            except:
                raise ValueError("The system has no box information. Please use "
                    + "the 'box' keyword argument.")
        else:
            if box is None:
                raise ValueError("Missing 'box' keyword argument!")
    else:
        if box is None:
            raise ValueError("Missing 'box' keyword argument!")

    if box is not None:
        if len(box) != 3:
            raise ValueError("The 'box' must have x, y, and z size information.")
        else:
            try:
                box = [float(x) for x in box]
            except:
                pass
            if not all(isinstance(x, float) for x in box):
                raise ValueError("The box dimensions must be of type 'float'")
            for size in box:
                if size < 0:
                    raise ValueError("The box size cannot be negative!.")

    return (molecule, box)

def _solvate(molecule, box, model, work_dir=None):
    """Internal function to add solvent using 'gmx solvate'.

       Positional arguments:

       molecule -- A molecule, or system of molecules.
       box      -- A list containing the box size in each dimension (in nm).
       model    -- The water model.

       Keyword arguments:

       work_dir -- The working directory for the process.
    """

    # Create a hash for the solvation run.
    #hash = hash((molecule, box)) % ((_sys.maxsize + 1) * 2)

    # Store the current working directory.
    dir = _os.getcwd()

    # Create a temporary working directory and store the directory name.
    if work_dir is None:
        tmp_dir = _tempfile.TemporaryDirectory()
        work_dir = tmp_dir.name

    # User specified working directory.
    else:
        # Create the directory if it doesn't already exist.
        if not _os.path.isdir(work_dir):
            _os.makedirs(work_dir)

    if work_dir is not None:
        # Change to the working directory for the process.
        # This avoid problems with relative paths.
        _os.chdir(work_dir)

    # Create the gmx command.
    command = "%s solvate -cs %s" % (_gmx_exe, model)

    if molecule is not None:
        # Write the molecule/system to a GRO and TOP files.
        _IO.saveMolecules("input", molecule, "gro87")
        _shutil.copyfile("input.gro87", "input.gro")

        # Update the command.
        command += " -cp input.gro"

        # Add the box information.
        if box is not None:
            command += " -box %f %f %f" % (box[0], box[1], box[2])

    # Just add box information.
    else:
        command += " -box %f %f %f" % (box[0], box[1], box[2])

    # Add the output file.
    command += " -o output.gro"

    # Run the command.
    proc = _subprocess.run(command, shell=True,
        stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)

    # gmx doesn't return sensible error codes, so we need to check that
    # the expected output was generated.
    if _os.path.isfile("output.gro"):
        system = _IO.readMolecules("output.gro")

        # Change back to the original directory.
        if work_dir is not None:
            _os.chdir(dir)

        return system
