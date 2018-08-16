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
A charge type.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

import Sire.Units as _Units

from ._type import Type as _Type

__all__ = ["Charge"]

class Charge(_Type):
    # Dictionary of allowed units.
    _supported_units = { "ELECTRON CHARGE" : _Units.e_charge,
                         "COULOMB"         : _Units.coulomb }

    # Map unit abbreviations to the full name.
    _abbreviations = { "E" : "ELECTRON CHARGE",
                       "C" : "COULOMB" }

    # Print formatting.
    _print_format = { "ELECTRON CHARGE" : "|e|",
                      "COULOMB"         : "C" }

    def __init__(self, *args):
        """Constructor.

           Positional arguments
           --------------------

           magnitude : float
               The magnitude.

           unit : str
               The unit.

           or

           string : str
               A string representation of the charge.
        """

        # Call the base class constructor.
        super().__init__(*args)

    def electron_charge(self):
        """Return the energy in electron charge."""
        return Charge((self._magnitude * self._supported_units[self._unit]).to(_Units.e_charge), "ELECTRON CHARGE")

    def coulomb(self):
        """Return the energy in Coulomb."""
        return Charge((self._magnitude * self._supported_units[self._unit]).to(_Units.coulomb), "COULOMB")

    def _default_unit(self, mag=None):
        """Internal method to return an object of the same type in the default unit.

           Positional argument:

           mag -- The magnitude (optional).
        """
        if mag is None:
            return self.electron_charge()
        else:
            return Charge(mag, "ELECTRON CHARGE")

    def _convert_to(self, unit):
        """Return the charge in a different unit.

           Positional arguments:

           unit -- The unit to convert to.
        """
        if unit == "ELECTRON CHARGE":
            return self.electron_charge()
        elif unit == "COULOMB":
            return self.coulomb()
        else:
            raise ValueError("Supported units are: '%s'" % list(self._supported_units.keys()))

    def _validate_unit(self, unit):
        """Validate that the unit are supported."""

        # Strip whitespace and convert to upper case.
        unit = unit.replace(" ", "").upper()

        # Replace all instances of ELECTRON with E.
        unit = unit.replace("ELECTRON", "E")

        # Strip all instances of CHARGE.
        unit = unit.replace("CHARGE", "");

        # Strip all instances of -.
        unit = unit.replace("-", "")

        # Strip all instances of S.
        unit = unit.replace("S", "")

        # Replace all instance of COULOMB with C.
        unit = unit.replace("COULOMB", "C")

        # Check that the unit is supported.
        if unit in self._abbreviations:
            return self._abbreviations[unit]
        else:
            raise ValueError("Supported units are: '%s'" % list(self._supported_units.keys()))
