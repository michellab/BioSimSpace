######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2019
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
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["solvate", "spc", "spce", "tip3p", "tip4p", "tip5p", "waterModels"]

import os as _os
import re as _re
import subprocess as _subprocess
import sys as _sys
import tempfile as _tempfile
import warnings as _warnings

from Sire import Base as _SireBase
from Sire import IO as _SireIO
from Sire import Mol as _SireMol

from BioSimSpace import _gmx_exe, _gromacs_path

from BioSimSpace._Exceptions import MissingSoftwareError as _MissingSoftwareError
from BioSimSpace._SireWrappers import System as _System
from BioSimSpace._SireWrappers import Molecule as _Molecule
from BioSimSpace._SireWrappers import Molecules as _Molecules
from BioSimSpace.Types import Coordinate as _Coordinate
from BioSimSpace.Types import Length as _Length

from BioSimSpace import IO as _IO
from BioSimSpace import _Utils as _Utils

def solvate(model, molecule=None, box=None, shell=None,
        ion_conc=0, is_neutral=True, work_dir=None, property_map={}):
    """Add SPC solvent.

       Parameters
       ----------

       model : str
           The name of the water model.

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                  :class:`Molecule <BioSimSpace._SireWrappers.Molecules>`, \
                  :class:`System <BioSimSpace._SireWrappers.System>`
           A molecule, or container/system of molecules.

       box : [:class:`Length <BioSimSpace.Types.Length>`]
           A list containing the box size in each dimension.

       shell : :class:`Length` <BioSimSpace.Types.Length>`
           Thickness of the water shell around the solute. Note that the
           base length of the resulting box must be at least twice as large
           as the cutoff used by the chosen molecular dynamics engine. As such,
           the shell option is often unsuitable for small molecules.

       ion_conc : float
           The ion concentration in (mol per litre).

       is_neutral : bool
           Whether to neutralise the system.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       system : :class:`System <BioSimSpace._SireWrappers.System>`
           The solvated molecular system.
    """

    if type(model) is not str:
        raise TypeError("'model' must be of type 'str'")
    else:
        # Strip whitespace and convert to lower case.
        model = model.replace(" ", "").lower()

        # Check that this water model is supported.
        if model not in _models_lower:
            raise ValueError("Supported water models are: %s" % waterModels())

    return _model_dict[model](molecule, box, shell, ion_conc, is_neutral, work_dir, property_map)

def spc(molecule=None, box=None, shell=None, ion_conc=0,
        is_neutral=True, work_dir=None, property_map={}):
    """Add SPC solvent.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                  :class:`Molecule <BioSimSpace._SireWrappers.Molecules>`, \
                  :class:`System <BioSimSpace._SireWrappers.System>`
           A molecule, or container/system of molecules.

       box : [:class:`Length <BioSimSpace.Types.Length>`]
           A list containing the box size in each dimension.

       shell : :class:`Length` <BioSimSpace.Types.Length>`
           Thickness of the water shell around the solute. Note that the
           base length of the resulting box must be at least twice as large
           as the cutoff used by the chosen molecular dynamics engine. As such,
           the shell option is often unsuitable for small molecules.

       ion_conc : float
           The ion concentration in (mol per litre).

       is_neutral : bool
           Whether to neutralise the system.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       system : :class:`System <BioSimSpace._SireWrappers.System>`
           The solvated molecular system.
    """

    if _gmx_exe is None or _gromacs_path is None:
        raise _MissingSoftwareError("'BioSimSpace.Solvent.spc' is not supported. "
                                    "Please install GROMACS (http://www.gromacs.org).")

    # Validate arguments.
    molecule, box, shell, work_dir, property_map = \
        _validate_input(molecule, box, shell, ion_conc, is_neutral, work_dir, property_map)

    # Create the solvated system.
    return _solvate(molecule, box, shell, "spc", 3, ion_conc,
            is_neutral, work_dir=work_dir, property_map=property_map)

def spce(molecule=None, box=None, shell=None, ion_conc=0, is_neutral=True,
        work_dir=None, property_map={}):
    """Add SPC/E solvent.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                  :class:`Molecule <BioSimSpace._SireWrappers.Molecules>`, \
                  :class:`System <BioSimSpace._SireWrappers.System>`
           A molecule, or container/system of molecules.

       box : [:class:`Length <BioSimSpace.Types.Length>`]
           A list containing the box size in each dimension.

       shell : :class:`Length` <BioSimSpace.Types.Length>`
           Thickness of the water shell around the solute. Note that the
           base length of the resulting box must be at least twice as large
           as the cutoff used by the chosen molecular dynamics engine. As such,
           the shell option is often unsuitable for small molecules.

       ion_conc : float
           The ion concentration in (mol per litre).

       is_neutral : bool
           Whether to neutralise the system.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       system : :class:`System <BioSimSpace._SireWrappers.System>`
           The solvated molecular system.
    """

    if _gmx_exe is None:
        raise _MissingSoftwareError("'BioSimSpace.Solvent.spce' is not supported. "
                                    "Please install GROMACS (http://www.gromacs.org).")

    # Validate arguments.
    molecule, box, shell, work_dir, property_map = \
        _validate_input(molecule, box, shell, ion_conc, is_neutral, work_dir, property_map)

    # Create the solvated system.
    return _solvate(molecule, box, shell, "spce", 3, ion_conc,
            is_neutral, work_dir=work_dir, property_map=property_map)

def tip3p(molecule=None, box=None, shell=None, ion_conc=0,
        is_neutral=True, work_dir=None, property_map={}):
    """Add TIP3P solvent.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                  :class:`Molecule <BioSimSpace._SireWrappers.Molecules>`, \
                  :class:`System <BioSimSpace._SireWrappers.System>`
           A molecule, or container/system of molecules.

       box : [:class:`Length <BioSimSpace.Types.Length>`]
           A list containing the box size in each dimension.

       shell : :class:`Length` <BioSimSpace.Types.Length>`
           Thickness of the water shell around the solute. Note that the
           base length of the resulting box must be at least twice as large
           as the cutoff used by the chosen molecular dynamics engine. As such,
           the shell option is often unsuitable for small molecules.

       ion_conc : float
           The ion concentration in (mol per litre).

       is_neutral : bool
           Whether to neutralise the system.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       system : :class:`System <BioSimSpace._SireWrappers.System>`
           The solvated molecular system.
    """

    if _gmx_exe is None:
        raise _MissingSoftwareError("'BioSimSpace.Solvent.tip3p' is not supported. "
                                    "Please install GROMACS (http://www.gromacs.org).")

    # Validate arguments.
    molecule, box, shell, work_dir, property_map = \
        _validate_input(molecule, box, shell, ion_conc, is_neutral, work_dir, property_map)

    # Create the solvated system.
    return _solvate(molecule, box, shell, "tip3p", 3, ion_conc,
            is_neutral, work_dir=work_dir, property_map=property_map)

def tip4p(molecule=None, box=None, shell=None, ion_conc=0,
        is_neutral=True, work_dir=None, property_map={}):
    """Add TIP4P solvent.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                  :class:`Molecule <BioSimSpace._SireWrappers.Molecules>`, \
                  :class:`System <BioSimSpace._SireWrappers.System>`
           A molecule, or container/system of molecules.

       box : [:class:`Length <BioSimSpace.Types.Length>`]
           A list containing the box size in each dimension.

       shell : :class:`Length` <BioSimSpace.Types.Length>`
           Thickness of the water shell around the solute. Note that the
           base length of the resulting box must be at least twice as large
           as the cutoff used by the chosen molecular dynamics engine. As such,
           the shell option is often unsuitable for small molecules.

       ion_conc : float
           The ion concentration in (mol per litre).

       is_neutral : bool
           Whether to neutralise the system.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       system : :class:`System <BioSimSpace._SireWrappers.System>`
           The solvated molecular system.
    """

    if _gmx_exe is None:
        raise _MissingSoftwareError("'BioSimSpace.Solvent.tip4p' is not supported. "
                                    "Please install GROMACS (http://www.gromacs.org).")

    # Validate arguments.
    molecule, box, shell, work_dir, property_map = \
        _validate_input(molecule, box, shell, ion_conc, is_neutral, work_dir, property_map)

    # Return the solvated system.
    return _solvate(molecule, box, shell, "tip4p", 4, ion_conc,
            is_neutral, work_dir=work_dir, property_map=property_map)

def tip5p(molecule=None, box=None, shell=None, ion_conc=0,
        is_neutral=True, work_dir=None, property_map={}):
    """Add TIP5P solvent.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                  :class:`Molecule <BioSimSpace._SireWrappers.Molecules>`, \
                  :class:`System <BioSimSpace._SireWrappers.System>`
           A molecule, or container/system of molecules.

       box : [:class:`Length <BioSimSpace.Types.Length>`]
           A list containing the box size in each dimension.

       shell : :class:`Length` <BioSimSpace.Types.Length>`
           Thickness of the water shell around the solute. Note that the
           base length of the resulting box must be at least twice as large
           as the cutoff used by the chosen molecular dynamics engine. As such,
           the shell option is often unsuitable for small molecules.

       ion_conc : float
           The ion concentration in (mol per litre).

       is_neutral : bool
           Whether to neutralise the system.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       system : :class:`System <BioSimSpace._SireWrappers.System>`
           The solvated molecular system.
    """

    if _gmx_exe is None:
        raise _MissingSoftwareError("'BioSimSpace.Solvent.tip5p' is not supported. "
                                    "Please install GROMACS (http://www.gromacs.org).")

    # Validate arguments.
    molecule, box, shell, work_dir, property_map = \
        _validate_input(molecule, box, shell, ion_conc, is_neutral, work_dir, property_map)

    # Return the solvated system.
    return _solvate(molecule, box, shell, "tip5p", 5, ion_conc,
            is_neutral, work_dir=work_dir, property_map=property_map)

def _validate_input(molecule, box, shell, ion_conc, is_neutral, work_dir, property_map):
    """Internal function to validate function arguments.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                  :class:`Molecule <BioSimSpace._SireWrappers.Molecules>`, \
                  :class:`System <BioSimSpace._SireWrappers.System>`
           A molecule, or container/system of molecules.

       box : [:class:`Length <BioSimSpace.Types.Length>`]
           A list containing the box size in each dimension.

       shell : :class:`Length` <BioSimSpace.Types.Length>`
           Thickness of the water shell around the solute. Note that the
           base length of the resulting box must be at least twice as large
           as the cutoff used by the chosen molecular dynamics engine. As such,
           the shell option is often unsuitable for small molecules.

       ion_conc : float
           The ion concentration in (mol per litre).

       is_neutral : bool
           Whether to neutralise the system.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       (molecule, box, shell, work_dir, property_map) : tuple
           The validated input arguments.
    """

    # Whether to check the box size.
    check_box = True

    # Validate the molecule and create a local copy called _molecule to ensure
    # that the passed molecule is preserved.
    if molecule is not None:
        if type(molecule) is _Molecule:
            _molecule = _Molecule(molecule)
        elif type(molecule) is _Molecules:
            _molecule = molecule.toSystem()
        elif type(molecule) is _System:
            _molecule = _System(molecule)
        else:
            raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule' "
                            "'BioSimSpace._SireWrappers.Molecules', or 'BioSimSpace._SireWrappers.System'")

        # Try to extract the box dimensions from the system.
        if type(molecule) is _System and box is None and shell is None:
            try:
                check_box = False
                prop = property_map.get("space", "space")
                box = molecule._sire_object.property(prop).dimensions()
                # Convert to a list of Length objects.
                box = [_Length(box[0], "A"), _Length(box[1], "A"), _Length(box[2], "A")]
            except:
                raise ValueError("The system has no box information. Please use "
                                 "the 'box' keyword argument.")
        else:
            if box is None and shell is None:
                raise ValueError("Missing 'box' keyword argument!")
    else:
        _molecule = None

        if box is None:
            raise ValueError("Missing 'box' keyword argument!")

        if shell is not None:
            _warnings.warn("Ignoring 'shell' keyword argument as solute is missing.")
            shell = None

    if box is not None:
        # Convert tuple to list.
        if type(box) is tuple:
            box = list(box)

        # Convert Coordinate to list.
        if type(box) is _Coordinate:
            box = [box.x(), box.y(), box.z()]

        # Validate.
        if len(box) != 3:
            raise ValueError("The 'box' must have x, y, and z size information.")
        else:
            if not all(isinstance(x, _Length) for x in box):
                raise ValueError("The box dimensions must be of type 'BioSimSpace.Types.Length'")

    if shell is not None:
        if type(shell) is not _Length:
            raise ValueError("'shell' must must be of type 'BioSimSpace.Types.Length'")

        if box is not None:
            _warnings.warn("Ignoring 'box' keyword argument as 'shell' takes precendence.")

        # Work out the box size based on axis-aligned bounding box.
        # We take the maximum dimension as the base length of our box.
        base_length = max(2 * molecule._getAABox().halfExtents())

        # Now add the shell thickness.
        base_length = _Length(base_length, "A") + shell

        # If we need to add ions, make sure the box is at least 2.54 nanometers
        # wide, i.e. twice the rlist cutoff used by GROMACS protocols.
        if base_length < _Length(2.54, "nm"):
            base_length = _Length(2.54, "nm")

        # Create the dimensions for a cubic box.
        box = 3*[base_length]

    # Check that the ion concentration is valid.
    if type(ion_conc) is not float and type(ion_conc) is not int:
        raise TypeError("'ion_conc' must be of type 'int' or 'float'.")
    elif ion_conc < 0:
        raise ValueError("'ion_conc' cannot be negative!")

    if type(is_neutral) is not bool:
        raise TypeError("'is_neutral' must be of type 'bool'.")

    # Check that the working directory is valid.
    if work_dir is not None and type(work_dir) is not str:
        raise TypeError("'work_dir' must be of type 'str'")

    # Check that the property map is valid.
    if type(property_map) is not dict:
        raise TypeError("'property_map' must be of type 'dict'")

    # Check that the box is large enough to hold the molecule.
    if check_box:
        if molecule is not None and shell is None and not _check_box_size(molecule, box, property_map):
            raise ValueError("The 'box' is not large enough to hold the 'molecule'")

    return (_molecule, box, shell, work_dir, property_map)

def _solvate(molecule, box, shell, model, num_point,
        ion_conc, is_neutral, work_dir=None, property_map={}):
    """Internal function to add solvent using 'gmx solvate'.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                  :class:`System <BioSimSpace._SireWrappers.System>`
           A molecule, or system of molecules.

       box : [:class:`Length <BioSimSpace.Types.Length>`]
           A list containing the box size in each dimension.

       shell : :class:`Length` <BioSimSpace.Types.Length>`
           Thickness of the water shell around the solute.

       model : str
           The name of the water model.

       num_point : int
           The number of atoms in the water model.

       ion_conc : float
           The ion concentration in (mol per litre).

       is_neutral : bool
           Whether to neutralise the system.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }


       Returns
       -------

       system : :class:`System <BioSimSpace._SireWrappers.System>`
           The solvated system.
    """

    if molecule is not None:
        # Store the centre of the molecule.
        center = molecule._getAABox(property_map).center()

        # Work out the vector from the centre of the molecule to the centre of the
        # water box, converting the distance in each direction to Angstroms.
        vec = []
        for x, y in zip(box, center):
            vec.append(0.5*x.angstroms().magnitude() - y)

        # Translate the molecule. This allows us to create a water box
        # around the molecule.
        molecule.translate(vec, property_map)

        if type(molecule) is _System:

            # Reformat all of the water molecules so that they match the
            # expected GROMACS topology template.
            waters = _SireIO.setGromacsWater(molecule._sire_object.search("water"), model)

            # Convert to a BioSimSpace molecules container.
            waters = _Molecules(waters.toMolecules())

            # Remove the old water molecules then add those with the updated topology.
            molecule.removeWaterMolecules()
            molecule.addMolecules(waters)

    # Create a temporary working directory and store the directory name.
    if work_dir is None:
        tmp_dir = _tempfile.TemporaryDirectory()
        work_dir = tmp_dir.name

    # Run the solvation in the working directory.
    with _Utils.cd(work_dir):

        # Create the gmx command.
        if num_point == 3:
            mod = "spc216"
        else:
            mod = model
        command = "%s solvate -cs %s" % (_gmx_exe, mod)

        if molecule is not None:
            # Write the molecule/system to a GRO files.
            _IO.saveMolecules("input", molecule, "gro87")
            _os.rename("input.gro87", "input.gro")

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

        with open("README.txt", "w") as file:
            # Write the command to file.
            file.write("# gmx solvate was run with the following command:\n")
            file.write("%s\n" % command)

        # Create files for stdout/stderr.
        stdout = open("solvate.out", "w")
        stderr = open("solvate.err", "w")

        # Run gmx solvate as a subprocess.
        proc = _subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)
        stdout.close()
        stderr.close()

        # gmx doesn't return sensible error codes, so we need to check that
        # the expected output was generated.
        if not _os.path.isfile("output.gro"):
            raise RuntimeError("'gmx solvate failed to generate output!")

        # Extract the water lines from the GRO file.
        water_lines = []
        with open("output.gro", "r") as file:
            for line in file:
                if _re.search("SOL", line):
                    # Store the SOL atom record.
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
            raise ValueError("No water molecules were generated. Try increasing "
                            "the 'box' size or 'shell' thickness.")

        # Create a TOP file for the water model. By default we use the Amber03
        # force field to generate a dummy topology for the water model.
        with open("water_ions.top", "w") as file:
            file.write("#define FLEXIBLE 1\n\n")
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

        # Load the water box.
        water = _IO.readMolecules(["water.gro", "water_ions.top"])

        # Create a new system by adding the water to the original molecule.
        if molecule is not None:
            # Translate the molecule and water back to the original position.
            vec = [-x for x in vec]
            molecule.translate(vec, property_map)
            water.translate(vec)

            if type(molecule) is _System:
                # Extract the non-water molecules from the original system.
                non_waters = _Molecules(molecule.search("not water")._sire_object.toMolecules())

                # Create a system by adding these to the water molecules from
                # gmx solvate, which will include the original waters.
                system = non_waters.toSystem() + water

            else:
                system = molecule.toSystem() + water

            # Add all of the water box properties to the new system.
            for prop in water._sire_object.propertyKeys():
                prop = property_map.get(prop, prop)

                # Add the space property from the water system.
                system._sire_object.setProperty(prop, water._sire_object.property(prop))
        else:
            system = water

        # Now we add ions to the system and neutralise the charge.
        if ion_conc > 0 or is_neutral:

            # Write the molecule + water system to file.
            _IO.saveMolecules("solvated", system, "gro87")
            _IO.saveMolecules("solvated", system, "grotop")
            _os.rename("solvated.gro87", "solvated.gro")
            _os.rename("solvated.grotop", "solvated.top")

            # First write an mdp file.
            with open("ions.mdp", "w") as file:
                file.write("; Neighbour searching\n")
                file.write("cutoff-scheme           = Verlet\n")
                file.write("rlist                   = 1.1\n")
                file.write("pbc                     = xyz\n")
                file.write("verlet-buffer-tolerance = -1\n")
                file.write("\n; Electrostatics\n")
                file.write("coulombtype             = cut-off\n")
                file.write("\n; VdW\n")
                file.write("rvdw                    = 1.0\n")

            # Create the grompp command.
            command = "%s grompp -f ions.mdp -po ions.out.mdp -c solvated.gro -p solvated.top -o ions.tpr" % _gmx_exe

            with open("README.txt", "a") as file:
                # Write the command to file.
                file.write("\n# gmx grompp was run with the following command:\n")
                file.write("%s\n" % command)

            # Create files for stdout/stderr.
            stdout = open("grommp.out", "w")
            stderr = open("grommp.err", "w")

            # Run grompp as a subprocess.
            proc = _subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)
            stdout.close()
            stderr.close()

            # Flag whether to break out of the ion adding stage.
            is_break = False

            # Check for the tpr output file.
            if not _os.path.isfile("ions.tpr"):
                if shell is None:
                    raise RuntimeError("'gmx grommp' failed to generate output! "
                                       "Perhaps your box is too small?")
                else:
                    is_break = True
                    _warnings.warn("Unable to achieve target ion concentration, try using "
                                   "'box' option instead of 'shell'.")

            # Only continue if grommp was successful. This allows us to skip the remainder
            # of the code if the ion addition failed when the 'shell' option was chosen, i.e.
            # because the estimated simulation box was too small.
            if not is_break:
                is_break = False

                # The ion concentration is unset.
                if ion_conc == 0:
                    # Get the current molecular charge.
                    charge = system.charge()

                    # Round to the nearest integer value.
                    charge = round(charge.magnitude())

                    # Create the genion command.
                    command = "echo SOL | %s genion -s ions.tpr -o solvated_ions.gro -p solvated.top -neutral" % _gmx_exe

                    # Add enough counter ions to neutralise the charge.
                    if charge > 0:
                        command += " -nn %d" % abs(charge)
                    else:
                        command += " -np %d" % abs(charge)
                else:
                    # Create the genion command.
                    command = "echo SOL | %s genion -s ions.tpr -o solvated_ions.gro -p solvated.top -%s -conc %f" \
                        % (_gmx_exe, "neutral" if is_neutral else "noneutral", ion_conc)

                with open("README.txt", "a") as file:
                    # Write the command to file.
                    file.write("\n# gmx genion was run with the following command:\n")
                    file.write("%s\n" % command)

                # Create files for stdout/stderr.
                stdout = open("genion.out", "w")
                stderr = open("genion.err", "w")

                # Run genion as a subprocess.
                proc = _subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)
                stdout.close()
                stderr.close()

                # Check for the output GRO file.
                if not _os.path.isfile("solvated_ions.gro"):
                    if shell is None:
                        raise RuntimeError("'gmx genion' failed to add ions! Perhaps your box is too small?")
                    else:
                        is_break = True
                        _warnings.warn("Unable to achieve target ion concentration, try using "
                                       "'box' option instead of 'shell'.")

                if not is_break:
                    # Counters for the number of SOL, NA, and CL atoms.
                    num_sol = 0
                    num_na = 0
                    num_cl = 0

                    # We now need to loop through the GRO file to extract the lines
                    # corresponding to water or ion atoms.
                    water_ion_lines = []

                    with open("solvated_ions.gro", "r") as file:
                        for line in file:
                            # This is a Sodium atom.
                            if _re.search("NA", line):
                                water_ion_lines.append(line)
                                num_na += 1

                            # This is a Chlorine atom.
                            if _re.search("CL", line):
                                water_ion_lines.append(line)
                                num_cl += 1

                            # This is a water atom.
                            elif _re.search("SOL", line):
                                water_ion_lines.append(line)
                                num_sol += 1

                    # Add any box information. This is the last line in the GRO file.
                    water_ion_lines.append(line)

                    # Write a GRO file that contains only the water and ion atoms.
                    if len(water_ion_lines) - 1 > 0:
                        with open("water_ions.gro", "w") as file:
                            file.write("BioSimSpace %s water box\n" % model.upper())
                            file.write("%d\n" % (len(water_ion_lines)-1))

                            for line in water_ion_lines:
                                file.write("%s" % line)

                    # Ions have been added. Update the TOP file fo the water model
                    # with the new atom counts.
                    if num_na > 0 or num_cl > 0:
                        with open("water_ions.top", "w") as file:
                            file.write("#define FLEXIBLE 1\n\n")
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
                            file.write("SOL               %d\n" % (num_sol / num_point))
                            if num_na > 0:
                                file.write("NA                %d\n" % num_na)
                            if num_cl > 0:
                                file.write("CL                %d\n" % num_cl)

                    # Load the water/ion box.
                    water_ions = _IO.readMolecules(["water_ions.gro", "water_ions.top"])

                    # Create a new system by adding the water to the original molecule.
                    if molecule is not None:

                        if type(molecule) is _System:
                            # Extract the non-water molecules from the original system.
                            non_waters = _Molecules(molecule.search("not water")._sire_object.toMolecules())

                            # Create a system by adding these to the water and ion
                            # molecules from gmx solvate, which will include the
                            # original waters.
                            system = non_waters.toSystem() + water_ions
                        else:
                            system = molecule.toSystem() + water_ions

                        # Add all of the water molecules' properties to the new system.
                        for prop in water_ions._sire_object.propertyKeys():
                            prop = property_map.get(prop, prop)

                            # Add the space property from the water system.
                            system._sire_object.setProperty(prop, water_ions._sire_object.property(prop))
                    else:
                        system = water_ions

        # Store the name of the water model as a system property.
        system._sire_object.setProperty("water_model", _SireBase.wrap(model))

    return system

def _check_box_size(molecule, box, property_map={}):
    """Internal function to check that box is big enough for the molecule.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                  :class:`System <BioSimSpace._SireWrappers.System>`
           A molecule, or system of molecules.

       box : [:class:`Length <BioSimSpace.Types.Length>`]
           A list containing the box size in each dimension.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       is_okay : True
           Whether the box is large enough.
    """

    # Get the axis-aligned bounding box of the molecule/system.
    aabox = molecule._getAABox(property_map)

    # Calculate the box size in each dimension, storing each component as a
    # length in Angstroms.
    mol_box = [_Length(2*x," A") for x in aabox.halfExtents()]

    # Make sure the box is big enough in each dimension.
    for len1, len2 in zip(box, mol_box):
        if len1 < len2:
            return False

    # We made it this far, all dimensions are large enough.
    return True

def _rename_water_molecule(molecule):
    """Internal function to rename residues/atoms in a water molecule to match
       the naming conventions used by GROMACS.

       Parameters
       ----------

       molecule : Sire.Mol.Molecule
           A Sire Molecule object.

       Returns
       -------

       molecule : Sire.Mol.Molecule
           The updated Sire Molecule object.
    """

    # Make the molecule editable.
    molecule = molecule.edit()

    # In GROMACS, all water molecules must be given the residue label "SOL".
    # We extract all of the waters from the system and relabel the
    # residues as appropriate.
    #
    # We need to work out what to do if existing water molecules don't contain
    # all of the required atoms, e.g. if we have crystal water oxygen sites
    # from a PDB file, or if the system is to be solvated with a 4- or 5-point
    # water model (where we need to add virtual sites).

    # Update the molecule with the new residue name.
    molecule = molecule.residue(_SireMol.ResIdx(0))     \
                       .rename(_SireMol.ResName("SOL")) \
                       .molecule()

    # Index for the hydrogen atoms.
    hydrogen_idx = 1

    # Gromacs water models use HW1/HW2 for hydrogen atoms at OW for water.
    for atom in molecule.atoms():
        try:
            # Hydrogen.
            if atom.property("element") == _SireMol.Element("H"):
                molecule = molecule.atom(atom.number())                              \
                                   .rename(_SireMol.AtomName("HW%d" % hydrogen_idx)) \
                                   .molecule()
                hydrogen_idx += 1
            # Oxygen.
            elif atom.property("element") == _SireMol.Element("O"):
                molecule = molecule.atom(atom.number())             \
                                   .rename(_SireMol.AtomName("OW")) \
                                   .molecule()

        # Otherwise, try to infer the element from the atom name.
        except:
            # Strip all digits from the name.
            name = "".join([x for x in atom.name().value() if not x.isdigit()])

            # Remove any whitespace.
            name = name.replace(" ", "")

            # Try to infer the element.
            element = _SireMol.Element.biologicalElement(name)

            # Hydrogen.
            if element == _SireMol.Element("H"):
                molecule = molecule.atom(atom.number())                              \
                                   .rename(_SireMol.AtomName("HW%d" % hydrogen_idx)) \
                                   .molecule()
                hydrogen_idx += 1
            # Oxygen.
            elif element == _SireMol.Element("O"):
                molecule = molecule.atom(atom.number())             \
                                   .rename(_SireMol.AtomName("OW")) \
                                   .molecule()

    # Commit and return the updated molecule.
    return molecule.commit()

# Create a list of the water models names.
# This needs to come after all of the solvation functions.
_models = []        # List of water models (actual names).
_models_lower = []  # List of lower case names.
_model_dict = {}    # Mapping between lower case names and functions.
import sys as _sys
_namespace = _sys.modules[__name__]
for _var in dir():
    if _var[0] != "_" and _var != "solvate" and _var[0] != "M":
        _models.append(_var)
        _models_lower.append(_var.lower())
        _model_dict[_var.lower()] = getattr(_namespace, _var)
del _namespace
del _sys
del _var

def waterModels():
    """Return a list of the supported water models.

       Returns
       -------

       models : [str]
          A list of the supported water models.
    """
    return _models
