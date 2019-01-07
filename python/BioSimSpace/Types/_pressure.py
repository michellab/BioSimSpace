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
An energy type.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

import Sire.Units as _Units

from ._type import Type as _Type

__all__ = ["Pressure"]

class Pressure(_Type):
    # Dictionary of allowed units.
    _supported_units = { "ATMOSPHERE" : _Units.atm,
                         "BAR"        : _Units.bar }

    # Map unit abbreviations to the full name.
    _abbreviations = { "ATM" : "ATMOSPHERE",
                       "BAR" : "BAR" }

    # Print formatting.
    _print_format = { "ATMOSPHERE" : "atm",
                      "BAR"        : "bar" }

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
               A string representation of the pressure.
        """

        # Call the base class constructor.
        super().__init__(*args)

    def atm(self):
        """Return the atmospheric pressure."""
        return Pressure((self._magnitude * self._supported_units[self._unit]).to(_Units.atm), "ATMOSPHERE")

    def bar(self):
        """Return the pressure in bar."""
        return Pressure((self._magnitude * self._supported_units[self._unit]).to(_Units.bar), "BAR")

    def _default_unit(self, mag=None):
        """Internal method to return an object of the same type in the default unit.

           Positional argument:

           mag -- The magnitude (optional).
        """
        if mag is None:
            return self.atm()
        else:
            return Pressure(mag, "ATMOSPHERE")

    def _convert_to(self, unit):
        """Return the pressure in a different unit.


           Positional arguments:

           unit -- The unit to convert to.
        """
        if unit == "ATMOSPHERE":
            return self.atm()
        elif unit == "BAR":
            return self.bar()
        else:
            raise ValueError("Supported units are: '%s'" % list(self._supported_units.keys()))

    def _validate_unit(self, unit):
        """Validate that the unit are supported."""

        # Strip whitespace and convert to upper case.
        unit = unit.replace(" ", "").upper()

        # Strip all instances of PRESSURE.
        unit = unit.replace("PRESSURE", "")
        unit = unit.replace("PRESS", "")
        unit = unit.replace("PRES", "")

        # Replace all instance of ATMOSPHERIC/ATMOSPHERE with ATM.
        unit = unit.replace("ATMOSPHERIC", "ATM")
        unit = unit.replace("ATMOSPHERE", "ATM")

        # Strip all "S" characters.
        unit = unit.replace("S", "")

        # Check that the unit is supported.
        if unit in self._abbreviations:
            return self._abbreviations[unit]
        else:
            raise ValueError("Supported units are: '%s'" % list(self._supported_units.keys()))
