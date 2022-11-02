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
Functionality for parameterising molecules.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["parameterise",
           "ff03",
           "ff99",
           "ff99SB",
           "ff99SBildn",
           "ff14SB",
           "gaff",
           "gaff2",
           "forceFields",
           "amberForceFields",
           "openForceFields"]

from .. import _amber_home, _gmx_exe, _gmx_path, _isVerbose

from .._Exceptions import IncompatibleError as _IncompatibleError
from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from .._SireWrappers import Molecule as _Molecule
from ..Solvent import waterModels as _waterModels
from ..Types import Charge as _Charge
from .._Utils import _try_import, _have_imported
from .. import _Utils

from ._process import Process as _Process
from . import Protocol as _Protocol

def parameterise(molecule, forcefield, water_model=None, work_dir=None, property_map={}):
    """Parameterise a molecule using a specified force field.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
           The molecule to parameterise, either as a Molecule object or SMILES
           string.

       water_model : str
           The water model used to parameterise any structural ions. This
           will be ignored when it is not supported by the chosen force field,
           or when ions aren't present. Run 'BioSimSpace.Solvent.waterModels()'
           to see the supported water models.

       forcefield : str
           The force field. Run BioSimSpace.Parameters.forceFields() to get a
           list of the supported force fields.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The parameterised molecule.
    """

    if not isinstance(forcefield, str):
        raise TypeError("'forcefield' must be of type 'str'")
    else:
        # Strip whitespace and convert to lower case.
        forcefield = forcefield.replace(" ", "").lower()

        if forcefield not in _forcefields_lower:
            raise ValueError("Supported force fields are: %s" % forceFields())

    # Check whether the forcefield supports using a water model to parameterise
    # structural ions.
    if _requires_water_model[forcefield]:
        return _forcefield_dict[forcefield](molecule, water_model=water_model, work_dir=work_dir, property_map=property_map)
    else:
        return _forcefield_dict[forcefield](molecule, work_dir=work_dir, property_map=property_map)

def ff99(molecule, water_model=None, leap_commands=None, work_dir=None, property_map={}):
    """Parameterise using the ff99 force field.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
           The molecule to parameterise, either as a Molecule object or SMILES
           string.

       water_model : str
           The water model used to parameterise any structural ions.
           Run 'BioSimSpace.Solvent.waterModels()' to see the supported
           water models. This is ignored if ions are not present.

       leap_commands : [str]
           An optional list of extra commands for the LEaP program. These
           will be added after any default commands and can be used to, e.g.,
           load additional parameter files. When this option is set, we can no
           longer fall back on GROMACS's pdb2gmx.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The parameterised molecule.
    """

    if _amber_home is None and (_gmx_exe is None or _gmx_path is None):
        raise _MissingSoftwareError("'BioSimSpace.Parameters.ff99' is not supported. "
                                    "Please install AmberTools (http://ambermd.org) or "
                                    "GROMACS (http://www.gromacs.org).")

    # Validate arguments.
    _validate(molecule=molecule, water_model=water_model,
        check_ions=True, leap_commands=leap_commands, property_map=property_map)

    # Create a default protocol.
    protocol = _Protocol.FF99(water_model=water_model,
                              leap_commands=leap_commands,
                              property_map=property_map)

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, work_dir=work_dir, auto_start=True)

def ff99SB(molecule, water_model=None, leap_commands=None, work_dir=None, property_map={}):
    """Parameterise using the ff99SB force field.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
           The molecule to parameterise, either as a Molecule object or SMILES
           string.

       water_model : str
           The water model used to parameterise any structural ions.
           Run 'BioSimSpace.Solvent.waterModels()' to see the supported
           water models. This is ignored if ions are not present.

       leap_commands : [str]
           An optional list of extra commands for the LEaP program. These
           will be added after any default commands and can be used to, e.g.,
           load additional parameter files. When this option is set, we can no
           longer fall back on GROMACS's pdb2gmx.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The parameterised molecule.
    """

    if _amber_home is None and (_gmx_exe is None or _gmx_path is None):
        raise _MissingSoftwareError("'BioSimSpace.Parameters.ff99SB' is not supported. "
                                    "Please install AmberTools (http://ambermd.org) "
                                    "or GROMACS (http://www.gromacs.org).")

    # Validate arguments.
    _validate(molecule=molecule, water_model=water_model,
        check_ions=True, leap_commands=leap_commands, property_map=property_map)

    # Create a default protocol.
    protocol = _Protocol.FF99SB(water_model=water_model,
                                leap_commands=leap_commands,
                                property_map=property_map)

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, work_dir=work_dir, auto_start=True)

def ff99SBildn(molecule, water_model=None, leap_commands=None, work_dir=None, property_map={}):
    """Parameterise using the ff99SBildn force field.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
           The molecule to parameterise, either as a Molecule object or SMILES
           string.

       water_model : str
           The water model used to parameterise any structural ions.
           Run 'BioSimSpace.Solvent.waterModels()' to see the supported
           water models. This is ignored if ions are not present.

       leap_commands : [str]
           An optional list of extra commands for the LEaP program. These
           will be added after any default commands and can be used to, e.g.,
           load additional parameter files. When this option is set, we can no
           longer fall back on GROMACS's pdb2gmx.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The parameterised molecule.
    """

    if _amber_home is None and (_gmx_exe is None or _gmx_path is None):
        raise _MissingSoftwareError("'BioSimSpace.Parameters.ff99SBildn' is not supported. "
                                    "Please install AmberTools (http://ambermd.org) "
                                    "or GROMACS (http://www.gromacs.org).")

    # Validate arguments.
    _validate(molecule=molecule, water_model=water_model,
        check_ions=True, leap_commands=leap_commands, property_map=property_map)

    # Create a default protocol.
    protocol = _Protocol.FF99SBILDN(water_model=water_model,
                                    leap_commands=leap_commands,
                                    property_map=property_map)

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, work_dir=work_dir, auto_start=True)

def ff03(molecule, water_model=None, leap_commands=None, work_dir=None, property_map={}):
    """Parameterise using the ff03 force field.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
           The molecule to parameterise, either as a Molecule object or SMILES
           string.

       water_model : str
           The water model used to parameterise any structural ions.
           Run 'BioSimSpace.Solvent.waterModels()' to see the supported
           water models. This is ignored if ions are not present.

       leap_commands : [str]
           An optional list of extra commands for the LEaP program. These
           will be added after any default commands and can be used to, e.g.,
           load additional parameter files. When this option is set, we can no
           longer fall back on GROMACS's pdb2gmx.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The parameterised molecule.
    """

    if _amber_home is None and (_gmx_exe is None or _gmx_path is None):
        raise _MissingSoftwareError("'BioSimSpace.Parameters.ff03' is not supported. "
                                    "Please install AmberTools (http://ambermd.org) "
                                    "or GROMACS (http://www.gromacs.org).")

    # Validate arguments.
    _validate(molecule=molecule, water_model=water_model,
        check_ions=True, leap_commands=leap_commands, property_map=property_map)

    # Create a default protocol.
    protocol = _Protocol.FF03(water_model=water_model,
                              leap_commands=leap_commands,
                              property_map=property_map)

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, work_dir=work_dir, auto_start=True)

def ff14SB(molecule, water_model=None, leap_commands=None, work_dir=None, property_map={}):
    """Parameterise using the ff14SB force field.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
           The molecule to parameterise, either as a Molecule object or SMILES
           string.

       water_model : str
           The water model used to parameterise any structural ions.
           Run 'BioSimSpace.Solvent.waterModels()' to see the supported
           water models. This is ignored if ions are not present.

       leap_commands : [str]
           An optional list of extra commands for the LEaP program. These
           will be added after any default commands and can be used to, e.g.,
           load additional parameter files. When this option is set, we can no
           longer fall back on GROMACS's pdb2gmx.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The parameterised molecule.
    """

    if _amber_home is None:
        raise _MissingSoftwareError("'BioSimSpace.Parameters.ff14SB' is not supported. "
                                    "Please install AmberTools (http://ambermd.org).")

    # Validate arguments.
    _validate(molecule=molecule, water_model=water_model,
        check_ions=True, leap_commands=leap_commands, property_map=property_map)

    # Create a default protocol.
    protocol = _Protocol.FF14SB(water_model=water_model,
                                leap_commands=leap_commands,
                                property_map=property_map)

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, work_dir=work_dir, auto_start=True)

def gaff(molecule, work_dir=None, net_charge=None, charge_method="BCC", property_map={}):
    """Parameterise using the gaff force field.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
           The molecule to parameterise, either as a Molecule object or SMILES
           string.

       net_charge : int, :class:`Charge <BioSimSpace.Types.Charge>`
           The net charge on the molecule.

       charge_method : str
               The method to use when calculating atomic charges:
               "RESP", "CM2", "MUL", "BCC", "ESP", "GAS"

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The parameterised molecule.
    """

    if _amber_home is None:
        raise _MissingSoftwareError("'BioSimSpace.Parameters.gaff' is not supported. "
                                    "Please install AmberTools (http://ambermd.org).")

    # Validate arguments.
    _validate(molecule=molecule, check_ions=False, property_map=property_map)

    if net_charge is not None:
        # Get the value of the charge.
        if isinstance(net_charge, _Charge):
            net_charge = net_charge.value()

        if isinstance(net_charge, float):
            if net_charge % 1 != 0:
                raise ValueError("'net_charge' must be integer valued.")

        # Try to convert to int.
        try:
            net_charge = int(net_charge)
        except:
            raise TypeError("'net_charge' must be of type 'int', or `BioSimSpace.Types.Charge'")

    # Create a default protocol.
    protocol = _Protocol.GAFF(net_charge=net_charge,
                              charge_method=charge_method,
                              property_map=property_map)

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, work_dir=work_dir, auto_start=True)

def gaff2(molecule, work_dir=None, net_charge=None, charge_method="BCC", property_map={}):
    """Parameterise using the gaff force field.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
           The molecule to parameterise, either as a Molecule object or SMILES
           string.

       net_charge : int, :class:`Charge <BioSimSpace.Types.Charge>`
           The net charge on the molecule.

       charge_method : str
               The method to use when calculating atomic charges:
               "RESP", "CM2", "MUL", "BCC", "ESP", "GAS"

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The parameterised molecule.
    """

    if _amber_home is None:
        raise _MissingSoftwareError("'BioSimSpace.Parameters.gaff2' is not supported. "
                                    "Please install AmberTools (http://ambermd.org).")

    # Validate arguments.
    _validate(molecule=molecule, check_ions=False, property_map=property_map)

    if net_charge is not None:
        # Get the value of the charge.
        if isinstance(net_charge, _Charge):
            net_charge = net_charge.value()

        if isinstance(net_charge, float):
            if net_charge % 1 != 0:
                raise ValueError("'net_charge' must be integer valued.")

        # Try to convert to int.
        try:
            net_charge = int(net_charge)
        except:
            raise TypeError("'net_charge' must be of type 'int', or `BioSimSpace.Types.Charge'")

        if net_charge % 1 != 0:
            raise ValueError("'net_charge' must be integer valued.")

    # Create a default protocol.
    protocol = _Protocol.GAFF2(net_charge=net_charge,
                               charge_method=charge_method,
                               property_map=property_map)

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, work_dir=work_dir, auto_start=True)

def _parameterise_openff(molecule, forcefield, work_dir=None, property_map={}):
    """Parameterise a molecule using a force field from the Open Force Field
       Initiative.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
           The molecule to parameterise, either as a Molecule object or SMILES
           string.

       forcefield : str
           The force field. Run BioSimSpace.Parameters.openForceFields() to get a
           list of the supported force fields.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The parameterised molecule.
    """

    from sire.legacy.Base import findExe as _findExe

    try:
        _findExe("antechamber")
    except:
        raise _MissingSoftwareError(f"'{forcefield}' is not supported. AmberTools "
                                     "(http://ambermd.org) is needed for charge "
                                     "calculation and 'antechamber' executable "
                                     "must be in your PATH.") from None

    # Check the Antechamber version. Open Force Field requires Antechamber >= 20.0.
    try:
        # Antechamber returns an exit code of 1 when requesting version information.
        # As such, we wrap the call within a try-except block in case it fails.

        import shlex as _shlex
        import subprocess as _subprocess

        # Generate the command-line string. (Antechamber must be in the PATH,
        # so no need to use AMBERHOME.
        command = "antechamber -v"

        # Run the command as a subprocess.
        proc = _subprocess.run(_Utils.command_split(command), shell=False, text=True,
            stdout=_subprocess.PIPE, stderr=_subprocess.STDOUT)

        # Get stdout and split into lines.
        lines = proc.stdout.split("\n")

        # If present, version information is on line 1.
        string = lines[1]

        # Delete the welcome message.
        string = string.replace("Welcome to antechamber", "")

        # Extract the version and convert to float.
        version = float(string.split(":")[0])

        # The version is okay, enable Open Force Field support.
        if version >= 20:
            is_compatible = True
        # Disable Open Force Field support.
        else:
            is_compatible = False

        del _shlex
        del _subprocess

    # Something went wrong, disable Open Force Field support.
    except:
        is_compatible = False
        raise

    if not is_compatible:
        raise _IncompatibleError(f"'{forcefield}' requires Antechamber >= 20.0")

    # Validate arguments.

    if not isinstance(molecule, (_Molecule, str)):
        raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule' or 'str'")

    if not isinstance(forcefield, str):
        raise TypeError("'forcefield' must be of type 'str'")
    else:
        # Strip whitespace and convert to lower case.
        forcefield = forcefield.replace(" ", "").lower()

        if forcefield not in _forcefields_lower:
            raise ValueError("Supported force fields are: %s" % openForceFields())

    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'")

    # Create a default protocol.
    protocol = _Protocol.OpenForceField(forcefield, property_map=property_map)

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, work_dir=work_dir, auto_start=True)

# Create a list of the force field names.
# This needs to come after all of the force field functions.
_forcefields = []           # List of force fields (actual names).
_amber_forcefields = []     # List of the supported AMBER force fields.
_forcefields_lower = []     # List of lower case names.
_forcefield_dict = {}       # Mapping between lower case names and functions.
_requires_water_model = {}  # Whether a given forcefield requires a water model
                            # to parameterise structural ions.
import sys as _sys
_namespace = _sys.modules[__name__]
for _var in dir():
    if _var[0] != "_" and _var[0].upper() != "P":
        _forcefields.append(_var)
        _amber_forcefields.append(_var)
        _forcefields_lower.append(_var.lower())
        if _var not in ["gaff", "gaff2"]:
            _requires_water_model[_var.lower()] = True
        else:
            _requires_water_model[_var.lower()] = False
        _forcefield_dict[_var.lower()] = getattr(_namespace, _var)

# Wrapper function to dynamically generate functions with a given name.
# Here "name" refers to the name of a supported force field from the Open
# Force Field Initiative. The force field name has been "tidied" so that
# it conforms to sensible function naming standards, i.e. "-" and "."
# characters replaced by underscores.
def _make_function(name):
    def _function(molecule, work_dir=None, property_map={}):
        """Parameterise a molecule using the named force field from the
           Open Force Field initiative.

           Parameters
           ----------

           molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
               The molecule to parameterise, either as a Molecule object or SMILES
               string.

           work_dir : str
               The working directory for the process.

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }

           Returns
           -------

           molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
               The parameterised molecule.
        """
        return _parameterise_openff(molecule, name, work_dir, property_map)
    return _function

# Dynamically create functions for all available force fields from the Open
# Force Field Initiative.
_openforcefields = _try_import("openforcefields")

if _have_imported(_openforcefields):
    from glob import glob as _glob
    import os as _os
    _openff_dirs = _openforcefields.get_forcefield_dirs_paths()
    _open_forcefields = []
    # Loop over all force field directories.
    for _dir in _openff_dirs:
        # Glob all offxml files in the directory.
        _ff_list = _glob(f"{_dir}" + "/*.offxml")
        for _ff in _ff_list:
            # Get the force field name (base name minus extension).
            _base = _os.path.basename(_ff)
            _ff = _os.path.splitext(_base)[0]

            # Only include unconstrained force-fields since we need to go via
            # an intermediate ParmEd conversion. This means ParmEd must receive
            # bond parameters. We can then choose to constrain later, if required.
            # See, e.g: https://github.com/openforcefield/openff-toolkit/issues/603

            if "unconstrained" in _ff:
                # Append to the list of available force fields.
                _forcefields.append(_ff)
                _open_forcefields.append(_ff)

                # Create a sane function name, i.e. replace "-" and "."
                # characters with "_".
                _func_name = _ff.replace("-", "_")
                _func_name = _func_name.replace(".", "_")

                # Generate the function and bind it to the namespace.
                _function = _make_function(_ff)
                setattr(_namespace, _func_name, _function)

                # Expose the function to the user.
                __all__.append(_func_name)

                # Convert force field name to lower case and map to its function.
                _forcefields_lower.append(_ff.lower())
                _forcefield_dict[_ff.lower()] = getattr(_namespace, _func_name)
                _requires_water_model[_ff.lower()] = False

    # Clean up redundant attributes.
    del _base
    del _dir
    del _ff
    del _ff_list
    del _function
    del _func_name
    del _glob
    del _make_function
    del _openff_dirs
    del _os
    del _namespace
    del _sys
    del _var
elif _isVerbose():
    print("openforcefields not available as this module cannot be loaded.")

del _openforcefields


def forceFields():
    """Return a list of the supported force fields.

       Returns
       -------

       force_fields : [str]
           A list of the supported force fields.
    """
    return _forcefields

def amberForceFields():
    """Return a list of the supported AMBER force fields.

       Returns
       -------

       force_fields : [str]
           A list of the supported force fields.
    """
    return _amber_forcefields

def openForceFields():
    """Return a list of the supported force fields from the Open Force Field
       Initiative.

       Returns
       -------

       force_fields : [str]
           A list of the supported force fields.
    """
    return _open_forcefields

def _validate_water_model(water_model):
    """Internal helper function used to validate the water model chosen
       for parameterising structural ions.

       Parameters
       ----------

       water_model : str
           The chosen water model.

       Returns
       -------

       is_valid : bool
           Whether the water model is supported.
    """
    # Srip whitespace and convert to lower case.
    water_model = water_model.replace(" ", "").lower()

    # Check that we support this water model.
    if water_model in _waterModels():
        return True
    else:
        return False

def _has_ions(molecule, property_map={}):
    """Internal helper function to check whether a molecule contains ions.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           A molecule object.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       has_ions : bool
           Whether the molecule contains ions.

       ions: [str]
           A list of the ions that were found.
    """

    # Store a list of ionic elements.
    # (Taken from an AMBER leaprc.water.* file.)

    elements = ["F",
                "Cl",
                "Br",
                "I",
                "Li",
                "Na",
                "K",
                "Rb",
                "Cs",
                "Mg",
                "Tl",
                "Cu",
                "Ag",
                "Be",
                "Cu",
                "Ni",
                "Pt",
                "Zn",
                "Co",
                "Pd",
                "Ag",
                "Cr",
                "Fe",
                "Mg",
                "V",
                "Mn",
                "Hg",
                "Cd",
                "Yb",
                "Ca",
                "Sn",
                "Pb",
                "Eu",
                "Sr",
                "Sm",
                "Ba",
                "Ra",
                "Al",
                "Fe",
                "Cr",
                "In",
                "Tl",
                "Y",
                "La",
                "Ce",
                "Pr",
                "Nd",
                "Sm",
                "Eu",
                "Gd",
                "Tb",
                "Dy",
                "Er",
                "Tm",
                "Lu",
                "Hf",
                "Zr",
                "Ce",
                "U",
                "Pu",
                "Th"]

    # We need to search for ions individually since Sire can't
    # handle or'ed search strings beyond a certain size.

    # A list of ions that we've found.
    ions = []

    # Whether the molecule is a string.
    if isinstance(molecule, str):
        is_string = True
        molecule = molecule.upper()
    else:
        is_string = False

    for element in elements:
        if is_string:
            if element.upper() in molecule:
                ions.append(element)
        else:
            try:
                molecule.search(f"element {element}", property_map=property_map)
                ions.append(element)
            except:
                pass

    # Check whether we found any ions.
    if len(ions) > 0:
        return True, ions
    else:
        return False, ions

def _validate(molecule=None, water_model=None, check_ions=False,
        leap_commands=None, work_dir=None, property_map=None):
    """
    Internal function to validate arguments.

    Parameters
    ----------

    molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
        The molecule to parameterise, either as a Molecule object or SMILES
        string.

    water_model : str
        The water model used to parameterise any structural ions.
        Run 'BioSimSpace.Solvent.waterModels()' to see the supported
        water models. This is ignored if ions are not present.

    check_ions: bool
        Whether to check for the presence of structural ions. This is only
        required when parameterising with protein force fields.

    leap_commands : [str]
        An optional list of extra commands for the LEaP program. These
        will be added after any default commands and can be used to, e.g.,
        load additional parameter files. When this option is set, we can no
        longer fall back on GROMACS's pdb2gmx.

    work_dir : str
        The working directory for the process.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }
    """

    if molecule is not None:
        if not isinstance(molecule, (_Molecule, str)):
            raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule' or 'str'")

    if water_model is not None:
        if not isinstance(water_model, str):
            raise TypeError("'water_model' must be of type 'str'.")
        if not _validate_water_model(water_model):
            water_models = ", ".join(_waterModels())
            raise ValueError(f"'{water_model}' is unsupported. Supported models are: {water_models}")
    elif check_ions:
        has_ions, ions = _has_ions(molecule, property_map=property_map)
        if has_ions:
            ion_string = ", ".join(ions)
            raise ValueError(f"The molecule contains the following ions: {ion_string}. "
                              "Please choose a 'water_model' for the ion parameters.")

    if leap_commands is not None:
        if not isinstance(leap_commands, (list, tuple)):
            raise TypeError("'leap_commands' must be a 'list' of 'str' types.")
        else:
            if not all(isinstance(x, str) for x in leap_commands):
                raise TypeError("'leap_commands' must be a 'list' of 'str' types.")

    if property_map is not None:
        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")
