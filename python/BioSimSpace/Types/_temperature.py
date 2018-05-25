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

import Sire.Units as _Units

import re as _re

__all__ = ["Temperature"]

class Temperature:
    # Dictionary of allowed units.
    _supported_units = { "KELVIN"     : _Units.kelvin,
                         "CELSIUS"    : _Units.celsius,
                         "FAHRENHEIT" : _Units.fahrenheit }

    # Map unit abbreviations to the full name.
    _abbreviations  = { "K" : "KELVIN",
                        "C" : "CELSIUS",
                        "F" : "FAHRENHEIT" }

    def __init__(self, magnitude, unit):
        """Constructor.

           Positional arguments:

           magnitude -- The magnitude.
           unit      -- The unit.
        """

        # Check that the magnitude is valid.
        if type(magnitude) is int:
            self._magnitude = float(magnitude)
        elif type(magnitude) is float:
            self._magnitude = magnitude
        else:
            raise TypeError("'magnitude' must be of type 'int' or 'float'")

        # Don't support negative temperatures.
        if magnitude < 0:
            raise ValueError("The temperature cannot be negative!")

        # Check that the unit is supported.
        self._unit = self._validate_unit(unit)

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "%.2f %s" % (self._magnitude, self._unit[0].upper())

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return "BioSimSpace.Types.Temperature(%f, '%s')" % (self._magnitude, self._unit)

    def __add__(self, other):
        """Addition operator."""

        # Add the magnitudes in a common unit.
        mag = self.kelvin() + other.kelvin()

        # Return a new temperature object.
        return Temperature(mag, "KELVIN")

    def __sub__(self, other):
        """Subtraction operator."""

        # Subtract the magnitudes in a common unit.
        mag = self.kelvin() - other.kelvin()

        # Return a new temperature object.
        return Temperature(mag, "KELVIN")

    def __mul__(self, other):
        """Multiplication operator."""

        if type(other) is Temperature:
            # Multiply the magnitudes in a common unit.
            mag = self.kelvin() * other.kelvin()

            # Return a new temperature object.
            return Temperature(mag, "KELVIN")

        else:
            return self.__rmul__(other)

    def __rmul__(self, other):
        """Multiplication operator."""

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Only support multiplication by float.
        if type(other) is float:
            mag = self._magnitude * other
            return Temperature(mag, "KELVIN")

        else:
            raise NotImplementedError

    def __truediv__(self, other):
        """Division operator."""

        if type(other) is Temperature:
            # Divide the magnitudes in a common unit.
            mag = self.kelvin() / other.kelvin()

            # Return a new temperature object.
            return Temperature(mag, "KELVIN")

        else:
            # Convert int to float.
            if type(other) is int:
                other = float(other)

            # Only support division by float.
            if type(other) is float:
                mag = self._magnitude / other
                return Temperature(mag, "KELVIN")

            else:
                raise NotImplementedError

    def magnitude(self):
        """Return the magnitude."""
        return self._magnitude

    def unit(self):
        """Return the unit."""
        return self._unit

    def kelvin(self):
        """Return the temperature in Kelvin."""
        return (self._magnitude * self._supported_units[self._unit]).value()

    def celsius(self):
        """Return the temperature in Celsius."""
        return (self._magnitude * self._supported_units[self._unit]).to(_Units.celsius)

    def fahrenheit(self):
        """Return the temperature in Fahrenheit."""
        return (self._magnitude * self._supported_units[self._unit]).to(_Units.fahrenheit)

    def _validate_unit(self, unit):
        """Validate that the unit are supported."""

        # Strip whitespace and convert to upper case.
        unit = unit.replace(" ", "").upper()

        # Check that unit is supported.
        if not unit in self._supported_units:
            if not unit in self._abbreviations:
                raise ValueError("Supported units are: '%s'" % list(self._supported_units.keys()))
            else:
                unit = self._abbreviations[unit]

        return unit
