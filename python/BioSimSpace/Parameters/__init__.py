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
Functionality for parameterising molecular systems.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

from .._SireWrappers import System as _System
from .._SireWrappers import Molecule as _Molecule

from .._Utils import Tleap as _Tleap

import multiprocessing as _multiprocessing

def ff99(molecules):
    """Parameterise using the ff99 force field.

       Positional arguments:

       molecules -- The molecules to parameterise.
    """
    return _parameterise(molecules, "ff99")

def ff99SB(molecules):
    """Parameterise using the ff99SB force field.

       Positional arguments:

       molecules -- The molecules to parameterise.
    """
    return _parameterise(molecules, "ff99SB")

def ff03(molecules):
    """Parameterise using the ff03 force field.

       Positional arguments:

       molecules -- The molecules to parameterise.
    """
    return _parameterise(molecules, "ff03")

def ff14(molecules):
    """Parameterise using the ff14 force field.

       Positional arguments:

       molecules -- The molecules to parameterise.
    """
    return _parameterise(molecules, "ff14")

def ff14SB(molecules):
    """Parameterise using the ff14SB force field.

       Positional arguments:

       molecules -- The molecules to parameterise.
    """
    return _parameterise(molecules, "ff14SB")

def _parameterise(molecules, forcefield):
    """Internal function to parameterise a set of molecules using a given force field.

       Positional arguments:

       molecules  -- The molecules to parameterise.
       forcefield -- The force field to use.
    """

    if type(molecules) is _Molecule:
        return _Molecule(_Tleap.run(molecules._getSireMolecule(), forcefield))
