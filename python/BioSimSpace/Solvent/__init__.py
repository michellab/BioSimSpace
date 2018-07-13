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

from BioSimSpace import _gmx_exe, _gromacs_path

from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from .._SireWrappers import System as _System
from .._SireWrappers import Molecule as _Molecule
from ..Types import Length as _Length

import BioSimSpace.IO as _IO

import os as _os
import re as _re
import shutil as _shutil
import subprocess as _subprocess
import sys as _sys
import tempfile as _tempfile
import warnings as _warnings

__all__ = ["solvate", "spc", "spce", "tip3p", "tip4p", "tip5p", "waterModels"]

def solvate(model, molecule=None, box=None, shell=None, ion_conc=0, is_neutral=True, map={}):
    """Add SPC solvent.

       Positional arguments:

       model      -- The name of the water model.

       Keyword arguments:

       molecule   -- A molecule, or system of molecules.
       box        -- A list containing the box size in each dimension (in nm).
       shell      -- Thickness of the water shell around the solute.
       ion_conc   -- The ion concentration in (mol per litre).
       is_neutral -- Whether to neutralise the system.
       map        -- A dictionary that maps system "properties" to their user defined
                     values. This allows the user to refer to properties with their
                     own naming scheme, e.g. { "charge" : "my-charge" }
    """

    if type(model) is not str:
        raise TypeError("'model' must be of type 'str'")
    else:
        if model not in waterModels():
            raise ValueError("Supported water models are: %s" % waterModels())

    return _model_dict[model](molecule, box, shell, ion_conc, is_neutral, map)

def spc(molecule=None, box=None, shell=None, ion_conc=0, is_neutral=True, map={}):
    """Add SPC solvent.

       Keyword arguments:

       molecule   -- A molecule, or system of molecules.
       box        -- A list containing the box size in each dimension (in nm).
       shell      -- Thickness of the water shell around the solute.
       ion_conc   -- The ion concentration in (mol per litre).
       is_neutral -- Whether to neutralise the system.
       map        -- A dictionary that maps system "properties" to their user defined
                     values. This allows the user to refer to properties with their
                     own naming scheme, e.g. { "charge" : "my-charge" }
    """

    if _gmx_exe is None or _gromacs_path is None:
        raise _MissingSoftwareError("'BioSimSpace.Solvent.spc' is not supported. "
            + "Please install GROMACS (http://www.gromacs.org).")

    # Validate arguments.
    molecule, box, shell = _validate_input(molecule, box, shell, ion_conc, is_neutral, map)

    # Create the solvated system.
    return _solvate(molecule, box, shell, "spc", 3, ion_conc, is_neutral)

def spce(molecule=None, box=None, shell=None, ion_conc=0, is_neutral=True, map={}):
    """Add SPC/E solvent.

       Keyword arguments:

       molecule   -- A molecule, or system of molecules.
       box        -- A list containing the box size in each dimension (in nm).
       shell      -- Thickness of the water shell around the solute.
       ion_conc   -- The ion concentration in (mol per litre).
       is_neutral -- Whether to neutralise the system.
       map        -- A dictionary that maps system "properties" to their user defined
                     values. This allows the user to refer to properties with their
                     own naming scheme, e.g. { "charge" : "my-charge" }
    """

    if _gmx_exe is None:
        raise _MissingSoftwareError("'BioSimSpace.Solvent.spce' is not supported. "
            + "Please install GROMACS (http://www.gromacs.org).")

    # Validate arguments.
    molecule, box, shell = _validate_input(molecule, box, shell, ion_conc, is_neutral, map)

    # Create the solvated system.
    return _solvate(molecule, box, shell, "spce", 3, ion_conc, is_neutral)

def tip3p(molecule=None, box=None, shell=None, ion_conc=0, is_neutral=True, map={}):
    """Add TIP3P solvent.

       Keyword arguments:

       molecule   -- A molecule, or system of molecules.
       box        -- A list containing the box size in each dimension (in nm).
       shell      -- Thickness of the water shell around the solute.
       ion_conc   -- The ion concentration in (mol per litre).
       is_neutral -- Whether to neutralise the system.
       map        -- A dictionary that maps system "properties" to their user defined
                     values. This allows the user to refer to properties with their
                     own naming scheme, e.g. { "charge" : "my-charge" }
    """

    if _gmx_exe is None:
        raise _MissingSoftwareError("'BioSimSpace.Solvent.tip3p' is not supported. "
            + "Please install GROMACS (http://www.gromacs.org).")

    # Validate arguments.
    molecule, box, shell = _validate_input(molecule, box, shell, ion_conc, is_neutral, map)

    # Create the solvated system.
    return _solvate(molecule, box, shell, "tip3p", 3, ion_conc, is_neutral)

def tip4p(molecule=None, box=None, shell=None, ion_conc=0, is_neutral=True, map={}):
    """Add TIP4P solvent.

       Keyword arguments:

       molecule   -- A molecule, or system of molecules.
       box        -- A list containing the box size in each dimension (in nm).
       shell      -- Thickness of the water shell around the solute.
       ion_conc   -- The ion concentration in (mol per litre).
       is_neutral -- Whether to neutralise the system.
       map        -- A dictionary that maps system "properties" to their user defined
                     values. This allows the user to refer to properties with their
                     own naming scheme, e.g. { "charge" : "my-charge" }
    """

    if _gmx_exe is None:
        raise _MissingSoftwareError("'BioSimSpace.Solvent.tip4p' is not supported. "
            + "Please install GROMACS (http://www.gromacs.org).")

    # Validate arguments.
    molecule, box, shell = _validate_input(molecule, box, shell, ion_conc, is_neutral, map)

    # Return the solvated system.
    return _solvate(molecule, box, shell, "tip4p", 4, ion_conc, is_neutral)

def tip5p(molecule=None, box=None, shell=None, ion_conc=0, is_neutral=True, map={}):
    """Add TIP5P solvent.

       Keyword arguments:

       molecule   -- A molecule, or system of molecules.
       box        -- A list containing the box size in each dimension (in nm).
       shell      -- Thickness of the water shell around the solute.
       ion_conc   -- The ion concentration in (mol per litre).
       is_neutral -- Whether to neutralise the system.
       map        -- A dictionary that maps system "properties" to their user defined
                     values. This allows the user to refer to properties with their
                     own naming scheme, e.g. { "charge" : "my-charge" }
    """

    if _gmx_exe is None:
        raise _MissingSoftwareError("'BioSimSpace.Solvent.tip5p' is not supported. "
            + "Please install GROMACS (http://www.gromacs.org).")

    # Validate arguments.
    molecule, box, shell = _validate_input(molecule, box, shell, ion_conc, is_neutral, map)

    # Return the solvated system.
    return _solvate(molecule, box, shell, "tip5p", 5, ion_conc, is_neutral)

def _validate_input(molecule, box, shell, ion_conc, is_neutral, map):
    """Internal function to validate function arguments.

       Positional arguments:

       molecule   -- A molecule, or system of molecules.
       box        -- A list containing the box size in each dimension (in nm).
       shell      -- Thickness of the water shell around the solute.
       ion_conc   -- The ion concentration in (mol per litre).
       is_neutral -- Whether to neutralise the system.
       map        -- A dictionary that maps system "properties" to their user defined
                     values. This allows the user to refer to properties with their
                     own naming scheme, e.g. { "charge" : "my-charge" }
    """

    if molecule is not None:
        if type(molecule) is not _Molecule and type(molecule) is not _System:
            raise TypeError("'molecule' must be of type 'BioSimSpace.Molecule' "
                + "or 'BioSimSpace.System'")

        # Try to extract the box dimensions from the system.
        if type(molecule) is _System and box is None:
            try:
                if "space" in map:
                    prop = map["space"]
                else:
                    prop = "space"
                box = system.property(prop).dimensions()
                # Convert to a list of Length objects.
                box = [_Length(box[0], "A"), _Length(box[1], "A"), _Length(box[2], "A")]
            except:
                raise ValueError("The system has no box information. Please use "
                    + "the 'box' keyword argument.")
        else:
            if box is None:
                raise ValueError("Missing 'box' keyword argument!")
    else:
        if box is None:
            raise ValueError("Missing 'box' keyword argument!")

        if shell is not None:
            _warnings.warn("Ignoring 'shell' keyword argument as solute is missing.")
            shell = None

    if box is not None:
        if len(box) != 3:
            raise ValueError("The 'box' must have x, y, and z size information.")
        else:
            if not all(isinstance(x, _Length) for x in box):
                raise ValueError("The box dimensions must be of type 'BioSimSpace.Types.Length'")

    if shell is not None:
        if type(shell) is not _Length:
            raise ValueError("'shell' must must be of type 'BioSimSpace.Types.Length'")

    if type(map) is not dict:
        raise TypeError("'map' must be of type 'dict'")

    if type(ion_conc) is not float and type(ion_conc) is not int:
        raise TypeError("'ion_conc' must be of type 'int' or 'float'.")

    if type(is_neutral) is not bool:
        raise TypeError("'is_neutral' must be of type 'bool'.")

    # Check that the box is large enough to hold the molecule.
    if molecule is not None and not _check_box_size(molecule, box):
        raise ValueError("The 'box' is not large enough to hold the 'molecule'")

    return (molecule, box, shell)

def _solvate(molecule, box, shell, model, num_point,
        ion_conc, is_neutral, work_dir=None, map={}):
    """Internal function to add solvent using 'gmx solvate'.

       Positional arguments:

       molecule   -- A molecule, or system of molecules.
       box        -- A list containing the box size in each dimension (in nm).
       shell      -- Thickness of the water shell around the solute.
       model      -- The name of the water model.
       num_point  -- The number of points in the model.
       ion_conc   -- The ion concentration in (mol per litre).
       is_neutral -- Whether to neutralise the system.

       Keyword arguments:

       work_dir -- The working directory for the process.
       map      -- A dictionary that maps system "properties" to their user defined
                   values. This allows the user to refer to properties with their
                   own naming scheme, e.g. { "charge" : "my-charge" }
    """

    if molecule is not None:
        # Store the centre of the molecule.
        center = molecule._getAABox().center()

        # Work out the vector from the centre of the molecule to the centre of the
        # water box, converting the distance in each direction to Angstroms.
        vec = []
        for x, y in zip(box, center):
            vec.append(0.5*x.angstroms().magnitude() - y)

        # Translate the molecule. This allows us to create a water box
        # around the molecule.
        molecule.translate(vec)

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
    if num_point == 3:
        mod = "spc216"
    else:
        mod = model
    command = "%s solvate -cs %s" % (_gmx_exe, mod)

    if molecule is not None:
        # Write the molecule/system to a GRO files.
        _IO.saveMolecules("input", molecule, "gro87")
        _shutil.copyfile("input.gro87", "input.gro")

        # Update the command.
        command += " -cp input.gro"

        # Add the box information.
        if box is not None:
            command += " -box %f %f %f" % (box[0].nanometers().magnitude(),
                                           box[1].nanometers().magnitude(),
                                           box[2].nanometers().magnitude())

        # Add the shell information.
        if shell is not None:
            command += " -shell %f" % shell.nanometers().magnitude()

    # Just add box information.
    else:
        command += " -box %f %f %f" % (box[0].nanometers().magnitude(),
                                       box[1].nanometers().magnitude(),
                                       box[2].nanometers().magnitude())

    # Add the output file.
    command += " -o output.gro"

    with open("README.txt", "w") as f:
        # Write the command to file.
        f.write("# gmx solvate was run with the following command:\n")
        f.write("%s\n" % command)

    # Run the command.
    proc = _subprocess.run(command, shell=True,
        stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)

    # gmx doesn't return sensible error codes, so we need to check that
    # the expected output was generated.
    if not _os.path.isfile("output.gro"):
        _os.chdir(dir)
        raise RuntimeError("'gmx solvate failed to generate output!")

    # Extract the water lines from the GRO file.
    water_lines = []
    with open("output.gro", "r") as file:
        for line in file:
            if _re.search("SOL", line):
                water_lines.append(line)

        # Add any box information. This is the last line in the GRO file.
        water_lines.append(line)

    # Write a GRO file that contains only the water atoms.
    if len(water_lines) - 1 > 0:
        with open("water.gro", "w") as file:
            file.write("BioSimSpace %s water box\n" % model.upper())
            file.write("%d\n" % (len(water_lines)-1))

            for line in water_lines:
                file.write("%s" % line)
    else:
        if work_dir is not None:
            _os.chdir(dir)
        raise ValueError("No water molecules were generated. Try increasing "
            + "the 'box' size or 'shell' thickness.")

    # Create a TOP file for the water model. By default we use the Amber03
    # force field to generate a dummy topology for the water model.
    with open("water.top", "w") as file:
        file.write("; Include AmberO3 force field\n")
        file.write('#include "amber03.ff/forcefield.itp"\n\n')
        file.write("; Include %s water topology\n" % model.upper())
        file.write('#include "amber03.ff/%s.itp"\n\n' % model)
        file.write("; Include ions\n")
        file.write('#include "amber03.ff/ions.itp"\n\n')
        file.write("[ system ] \n")
        file.write("BioSimSpace %s water box\n\n" % model.upper())
        file.write("[ molecules ] \n")
        file.write(";molecule name    nr.\n")
        file.write("SOL               %d\n" % ((len(water_lines)-1) / num_point))

    # Now we add ions to the system.
    if ion_conc > 0:

        # First write an mdp file.
        with open("ions.mdp", "w") as file:
            file.write("; Neighbour searching\n")
            file.write("cutoff-scheme           = Verlet\n")
            file.write("rlist                   = 1.1\n")
            file.write("pbc                     = xyz\n")
            file.write("verlet-buffer-tolerance = -1\n")
            file.write("\n; Electrostatics\n")
            file.write("coulombtype             = PME\n")
            file.write("pme-order               = 4\n")
            file.write("fourierspacing          = 0.10\n")
            file.write("rcoulomb                = 1.0\n")
            file.write("\n; VdW\n")
            file.write("rvdw                    = 1.0\n")

        # Create the grompp command.
        command = "gmx grompp -f ions.mdp -po ions.out.mdp -c water.gro -p water.top -o ions.tpr"

        with open("README.txt", "a") as f:
            # Write the command to file.
            f.write("\n# gmx grompp was run with the following command:\n")
            f.write("%s\n" % command)

        # Run the command.
        proc = _subprocess.run(command, shell=True,
            stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)

        # Check for the tpr output file.
        if not _os.path.isfile("ions.tpr"):
            _os.chdir(dir)
            raise RuntimeError("'gmx grommp' failed to generate output! "
                + "Perhaps your box is too small?")

        # Create the genion command.
        command = "echo 2 | gmx genion -s ions.tpr -o water_ions.gro -p water.top -%s -conc %f" \
            % ("neutral" if is_neutral else "noneutral", ion_conc)

        # Now run genion using the ions.tpr file as input.
        proc = _subprocess.run(command, shell=True,
            stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)

        # Check for the tpr output file.
        if not _os.path.isfile("water_ions.gro"):
            _os.chdir(dir)
            raise RuntimeError("'gmx genion' failed to add ions! Perhaps your box is too small?")

        # Load the water plus ions box.
        water = _IO.readMolecules(["water_ions.gro", "water.top"])

    else:
        # Load the water box.
        water = _IO.readMolecules(["water.gro", "water.top"])

    if molecule is not None:
        # Translate the molecule and water back to the original position.
        vec = [-x for x in vec]
        molecule.translate(vec)
        water.translate(vec)

    # Create a new system by adding the water to the original molecule.
    if molecule is not None:
        if type(molecule) is _System:
            system = _System(molecule + water.getMolecules())
        else:
            system = molecule + water.getMolecules()

        if "space" in map:
            prop = map["space"]
        else:
            prop = "space"

        # Add the space property from the water system.
        system._sire_system.setProperty(prop, water._sire_system.property(prop))
    else:
        system = water

    # Change back to the original directory.
    if work_dir is not None:
        _os.chdir(dir)

    return system

def _check_box_size(molecule, box):
    """Internal function to check that box is big enough for the molecule.

       Positional arguments:

       molecule -- A molecule, or system of molecules.
       box      -- A list containing the box size in each dimension (in nm).
    """

    # Get the axis-aligned bounding box of the molecule/system.
    aabox = molecule._getAABox()

    # Calculate the box size in each dimension, storing each component as a
    # length in Angstroms.
    mol_box = [_Length(2*x," A") for x in aabox.halfExtents()]

    # Make sure the box is big enough in each dimension.
    for len1, len2 in zip(box, mol_box):
        if len1 < len2:
            return False

    # We made it this far, all dimensions are large enough.
    return True

# Create a list of the water models names.
# This needs to come after all of the solvation functions.
_models = []
_model_dict = {}
import sys as _sys
_namespace = _sys.modules[__name__]
for _var in dir():
    if _var[0] != "_" and _var != "solvate" and _var[0] != "M":
        _models.append(_var)
        _model_dict[_var] = getattr(_namespace, _var)
del(_namespace)
del(_sys)
del(_var)

def waterModels():
    "Return a list of the supported water models"
    return _models
