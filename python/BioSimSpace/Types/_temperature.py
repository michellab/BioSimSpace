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
A temperature type.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

import Sire.Units as _Units

from ._type import Type as _Type

import re as _re

__all__ = ["Temperature"]

class Temperature(_Type):
    # Dictionary of allowed units.
    _supported_units = { "KELVIN"     : _Units.kelvin,
                         "CELSIUS"    : _Units.celsius,
                         "FAHRENHEIT" : _Units.fahrenheit }

    # Map unit abbreviations to the full name.
    _abbreviations = { "K" : "KELVIN",
                       "C" : "CELSIUS",
                       "F" : "FAHRENHEIT" }

    # Print formatting.
    _print_format = { "KELVIN"     : "K",
                      "CELSIUS"    : "C",
                      "FAHRENHEIT" : "F" }

    def __init__(self, *args):
        """Constructor.

           Positional arguments:

           magnitude -- The magnitude.
           unit      -- The unit.

           or

           string    -- A string representation of the temperature.
        """

        # Call the base class constructor.
        super().__init__(*args)

        # Check that the temperature is above absolute zero.
        if self._kelvin() < 0:
            raise ValueError("The temperature cannot be less than absolute zero (0 Kelvin).")

    def _kelvin(self):
        """Return the magnitude of the temperature in Kelvin."""
        return (self._magnitude * self._supported_units[self._unit]).value()

    def kelvin(self):
        """Return the temperature in Kelvin."""
        return Temperature((self._magnitude * self._supported_units[self._unit]).value(), "KELVIN")

    def celsius(self):
        """Return the temperature in Celsius."""
        return Temperature((self._magnitude * self._supported_units[self._unit]).to(_Units.celsius), "CELSIUS")

    def fahrenheit(self):
        """Return the temperature in Fahrenheit."""
        return Temperature((self._magnitude * self._supported_units[self._unit]).to(_Units.fahrenheit), "FAHRENHEIT")

    def _default_unit(self, mag=None):
        """Internal method to return an object of the same type in the default unit.

           Positional argument:

           mag -- The magnitude (optional).
        """
        if mag is None:
            return self.kelvin()
        else:
            return Temperature(mag, "KELVIN")

    def _convert_to(self, unit):
        """Return the temperature in a different unit.

           Positional arguments:

           unit -- The unit to convert to.
        """
        if unit == "KELVIN":
            return self.kelvin()
        elif unit == "CELSIUS":
            return self.celsius()
        elif unit == "FAHRENHEIT":
            return self.fahrenheit()
        else:
            raise ValueError("Supported units are: '%s'" % list(self._supported_units.keys()))

    def _validate_unit(self, unit):
        """Validate that the unit are supported."""

        # Strip whitespace and convert to upper case.
        unit = unit.replace(" ", "").upper()

        # Check that the unit is supported.
        if unit in self._supported_units:
            return unit
        elif unit in self._abbreviations:
            return self._abbreviations[unit]
        else:
            raise ValueError("Supported units are: '%s'" % list(self._supported_units.keys()))
