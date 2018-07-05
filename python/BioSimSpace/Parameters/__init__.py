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
Functionality for parameterising molecules.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

from BioSimSpace import _amber_home, _gmx_exe, _gromacs_path

from .._SireWrappers import System as _System
from .._SireWrappers import Molecule as _Molecule

from ._process import Process as _Process
from . import Protocol as _Protocol

from warnings import warn as _warn

def parameterise(molecule, forcefield, options={}, map={}):
    """Parameterise using the a specified force field.

       Positional arguments:

       molecule   -- The molecule to parameterise.
       forcefield -- The force field to parameterise with.

       Keyword arguments:

       options  -- A dictionary of keyword options to override the protocol defaults.
       map      -- A dictionary that maps system "properties" to their user defined
                   values. This allows the user to refer to properties with their
                   own naming scheme, e.g. { "charge" : "my-charge" }
    """

    if type(forcefield) is not str:
        raise TypeError("'forcefield' must be of type 'str'")
    else:
        if forcefield not in forceFields():
            raise ValueError("Supported force fields are: %s" % forceFields())

    return _forcefield_dict[forcefield](molecule, options=options, map=map)

def ff99(molecule, options={}, map={}):
    """Parameterise using the ff99 force field.

       Positional arguments:

       molecule -- The molecule to parameterise.

       Keyword arguments:

       options  -- A dictionary of keyword options to override the protocol defaults.
       map      -- A dictionary that maps system "properties" to their user defined
                   values. This allows the user to refer to properties with their
                   own naming scheme, e.g. { "charge" : "my-charge" }
    """

    if _amber_home is None and (_gmx_exe is None or _gromacs_path is None):
        _warn("'BioSimSpace.Parameters.ff99' is not supported. Please install "
            + "AMBER (http://ambermd.org) or GROMACS (http://www.gromacs.org).")
        return None

    # Validate arguments.

    if type(molecule) is not _Molecule:
        raise TypeError("'molecule' must be of type 'BioSimSpace.Molecule'")

    if type(options) is not dict:
        raise TypeError("'options' must be of type 'dict'")

    if type(map) is not dict:
        raise TypeError("'map' must be of type 'dict'")

    # Create a default protocol.
    protocol = _Protocol.FF99(map=map)

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, autostart=True)

def ff99SB(molecule, options={}, map={}):
    """Parameterise using the ff99SB force field.

       Positional arguments:

       molecule -- The molecule to parameterise.

       Keyword arguments:

       options  -- A dictionary of keyword options to override the protocol defaults.
       map      -- A dictionary that maps system "properties" to their user defined
                   values. This allows the user to refer to properties with their
                   own naming scheme, e.g. { "charge" : "my-charge" }
    """

    if _amber_home is None and (_gmx_exe is None or _gromacs_path is None):
        _warn("'BioSimSpace.Parameters.ff99' is not supported. Please install "
            + "AMBER (http://ambermd.org) or GROMACS (http://www.gromacs.org).")
        return None

    # Validate arguments.

    if type(molecule) is not _Molecule:
        raise TypeError("'molecule' must be of type 'BioSimSpace.Molecule'")

    if type(options) is not dict:
        raise TypeError("'options' must be of type 'dict'")

    if type(map) is not dict:
        raise TypeError("'map' must be of type 'dict'")

    # Create a default protocol.
    protocol = _Protocol.FF99SB(map=map)

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, autostart=True)

def ff03(molecule, options={}, map={}):
    """Parameterise using the ff03 force field.

       Positional arguments:

       molecule -- The molecule to parameterise.

       Keyword arguments:

       options  -- A dictionary of keyword options to override the protocol defaults.
       map      -- A dictionary that maps system "properties" to their user defined
                   values. This allows the user to refer to properties with their
                   own naming scheme, e.g. { "charge" : "my-charge" }
    """

    if _amber_home is None and (_gmx_exe is None or _gromacs_path is None):
        _warn("'BioSimSpace.Parameters.ff99' is not supported. Please install "
            + "AMBER (http://ambermd.org) or GROMACS (http://www.gromacs.org).")
        return None

    # Validate arguments.

    if type(molecule) is not _Molecule:
        raise TypeError("'molecule' must be of type 'BioSimSpace.Molecule'")

    if type(options) is not dict:
        raise TypeError("'options' must be of type 'dict'")

    if type(map) is not dict:
        raise TypeError("'map' must be of type 'dict'")

    # Create a default protocol.
    protocol = _Protocol.FF03(map=map)

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, autostart=True)

def ff14SB(molecule, options={}, map={}):
    """Parameterise using the ff14SB force field.

       Positional arguments:

       molecule -- The molecule to parameterise.

       Keyword arguments:

       options  -- A dictionary of keyword options to override the protocol defaults.
       map      -- A dictionary that maps system "properties" to their user defined
                   values. This allows the user to refer to properties with their
                   own naming scheme, e.g. { "charge" : "my-charge" }
    """

    if _amber_home is None:
        _warn("'BioSimSpace.Parameters.ff14SB' is not supported. Please install "
            + "AMBER (http://ambermd.org).")
        return None

    # Validate arguments.

    if type(molecule) is not _Molecule:
        raise TypeError("'molecule' must be of type 'BioSimSpace.Molecule'")

    if type(options) is not dict:
        raise TypeError("'options' must be of type 'dict'")

    if type(map) is not dict:
        raise TypeError("'map' must be of type 'dict'")

    # Create a default protocol.
    protocol = _Protocol.FF14SB(map=map)

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, autostart=True)

def gaff(molecule, options={}, map={}):
    """Parameterise using the gaff force field.

       Positional arguments:

       molecule -- The molecule to parameterise.

       Keyword arguments:

       options  -- A dictionary of keyword options to override the protocol defaults.
       map      -- A dictionary that maps system "properties" to their user defined
                   values. This allows the user to refer to properties with their
                   own naming scheme, e.g. { "charge" : "my-charge" }
    """

    if _amber_home is None:
        _warn("'BioSimSpace.Parameters.gaff' is not supported. Please install "
            + "AMBER (http://ambermd.org).")
        return None

    # Validate arguments.

    if type(molecule) is not _Molecule:
        raise TypeError("'molecule' must be of type 'BioSimSpace.Molecule'")

    if type(options) is not dict:
        raise TypeError("'options' must be of type 'dict'")

    if type(map) is not dict:
        raise TypeError("'map' must be of type 'dict'")

    # Create a default protocol.
    protocol = _Protocol.GAFF(map=map)

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, autostart=True)

def gaff2(molecule, options={}, map={}):
    """Parameterise using the gaff force field.

       Positional arguments:

       molecule -- The molecule to parameterise.

       Keyword arguments:

       options  -- A dictionary of keyword options to override the protocol defaults.
       map      -- A dictionary that maps system "properties" to their user defined
                   values. This allows the user to refer to properties with their
                   own naming scheme, e.g. { "charge" : "my-charge" }
    """

    if _amber_home is None:
        _warn("'BioSimSpace.Parameters.gaff2' is not supported. Please install "
            + "AMBER (http://ambermd.org).")
        return None

    # Validate arguments.

    if type(molecule) is not _Molecule:
        raise TypeError("'molecule' must be of type 'BioSimSpace.Molecule'")

    if type(options) is not dict:
        raise TypeError("'options' must be of type 'dict'")

    if type(map) is not dict:
        raise TypeError("'map' must be of type 'dict'")

    # Create a default protocol.
    protocol = _Protocol.GAFF2(map=map)

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, autostart=True)

# Create a list of the force field names.
# This needs to come after all of the force field functions.
_forcefields = []
_forcefield_dict = {}
import sys as _sys
_namespace = _sys.modules[__name__]
for _var in dir():
    if _var[0] != "_" and _var[0].upper() != "P":
        _forcefields.append(_var)
        _forcefield_dict[_var] = getattr(_namespace, _var)
del(_namespace)
del(_sys)

def forceFields():
    "Return a list of the supported force fields"
    return _forcefields
