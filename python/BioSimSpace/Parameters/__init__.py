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

from .._SireWrappers import System as _System
from .._SireWrappers import Molecule as _Molecule

from ._process import Process as _Process
from . import Protocol as _Protocol

def ff99(molecule, options={}, background=True):
    """Parameterise using the ff99 force field.

       Positional arguments:

       molecule -- The molecule to parameterise.

       Keyword arguments:

       options  -- A dictionary of keyword options to override the protocol defaults.
    """

    # Validate arguments.

    if type(molecule) is not _Molecule:
        raise TypeError("'molecule' must be of type 'BioSimSpace.Molecule'")

    if type(options) is not dict:
        raise TypeError("'options' must be of type 'dict'")

    # Create a default protocol.
    protocol = _Protocol.FF99()

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, autostart=True)

def ff99SB(molecule, options={}, background=True):
    """Parameterise using the ff99SB force field.

       Positional arguments:

       molecule -- The molecule to parameterise.

       Keyword arguments:

       options  -- A dictionary of keyword options to override the protocol defaults.
    """

    # Validate arguments.

    if type(molecule) is not _Molecule:
        raise TypeError("'molecule' must be of type 'BioSimSpace.Molecule'")

    if type(options) is not dict:
        raise TypeError("'options' must be of type 'dict'")

    # Create a default protocol.
    protocol = _Protocol.FF99SB()

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, autostart=True)

def ff03(molecule, options={}, background=True):
    """Parameterise using the ff03 force field.

       Positional arguments:

       molecule -- The molecule to parameterise.

       Keyword arguments:

       options  -- A dictionary of keyword options to override the protocol defaults.
    """

    # Validate arguments.

    if type(molecule) is not _Molecule:
        raise TypeError("'molecule' must be of type 'BioSimSpace.Molecule'")

    if type(options) is not dict:
        raise TypeError("'options' must be of type 'dict'")

    # Create a default protocol.
    protocol = _Protocol.FF03()

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, autostart=True)

def ff14SB(molecule, options={}, background=True):
    """Parameterise using the ff14SB force field.

       Positional arguments:

       molecule -- The molecule to parameterise.

       Keyword arguments:

       options  -- A dictionary of keyword options to override the protocol defaults.
    """

    # Validate arguments.

    if type(molecule) is not _Molecule:
        raise TypeError("'molecule' must be of type 'BioSimSpace.Molecule'")

    if type(options) is not dict:
        raise TypeError("'options' must be of type 'dict'")

    # Create a default protocol.
    protocol = _Protocol.FF14SB()

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, autostart=True)

def gaff(molecule, options={}, background=True):
    """Parameterise using the gaff force field.

       Positional arguments:

       molecule -- The molecule to parameterise.

       Keyword arguments:

       options  -- A dictionary of keyword options to override the protocol defaults.
    """

    # Validate arguments.

    if type(molecule) is not _Molecule:
        raise TypeError("'molecule' must be of type 'BioSimSpace.Molecule'")

    if type(options) is not dict:
        raise TypeError("'options' must be of type 'dict'")

    # Create a default protocol.
    protocol = _Protocol.GAFF()

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, autostart=True)

def gaff2(molecule, options={}, background=True):
    """Parameterise using the gaff force field.

       Positional arguments:

       molecule -- The molecule to parameterise.

       Keyword arguments:

       options  -- A dictionary of keyword options to override the protocol defaults.
    """

    # Validate arguments.

    if type(molecule) is not _Molecule:
        raise TypeError("'molecule' must be of type 'BioSimSpace.Molecule'")

    if type(options) is not dict:
        raise TypeError("'options' must be of type 'dict'")

    # Create a default protocol.
    protocol = _Protocol.GAFF2()

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, autostart=True)

# Create a list of the force field names.
# This needs to come after all of the force field functions.
_forcefields = []
for _var in dir():
    if _var[0] != "_":
        _forcefields.append(_var)

def forceFields():
    "Print a list of the supported force fields"
    print(", ".join(_forcefields))
