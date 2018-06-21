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

__all__ = ["Temperature"]

class Temperature:
    # Dictionary of allowed units.
    _supported_units = { "KELVIN"     : _Units.kelvin,
                         "CELSIUS"    : _Units.celsius,
                         "FAHRENHEIT" : _Units.fahrenheit }

    # Map unit abbreviations to the full name.
    _abbreviations = { "K" : "KELVIN",
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

        # Check that the unit is supported.
        self._unit = self._validate_unit(unit)

        # Check that the temperature is above absolute zero.
        if self._kelvin() < 0:
            raise ValueError("The temperature cannot be less than absolute zero (0 Kelvin).")

    def __str__(self):
        """Return a human readable string representation of the object."""
        if self._magnitude > 1e6 or abs(self._magnitude) < 1e-6:
            return "%.4e %s" % (self._magnitude, self._unit[0].upper())
        else:
            return "%.2f %s" % (self._magnitude, self._unit[0].upper())

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._magnitude > 1e6 or abs(self._magnitude) < 1e-6:
            return "BioSimSpace.Types.Temperature(%.4e, '%s')" % (self._magnitude, self._unit)
        else:
            return "BioSimSpace.Types.Temperature(%f, '%s')" % (self._magnitude, self._unit)

    def __add__(self, other):
        """Addition operator."""

        # Add the magnitudes in a common unit.
        mag = self.kelvin().magnitude() + other.kelvin().magnitude()

        # Get new magnitude in the original unit.
        # Left-hand operand takes precedence.
        mag = Temperature(mag, "KELVIN")._convert_to(self._unit).magnitude()

        # Return a new temperature object.
        return Temperature(mag, self._unit)

    def __sub__(self, other):
        """Subtraction operator."""

        # Subtract the magnitudes in a common unit.
        mag = self.kelvin().magnitude() - other.kelvin().magnitude()

        # Get new magnitude in the original unit.
        # Left-hand operand takes precedence.
        mag = Temperature(mag, "KELVIN")._convert_to(self._unit).magnitude()

        # Return a new temperature object.
        return Temperature(mag, self._unit)

    def __mul__(self, other):
        """Multiplication operator."""

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Only support multiplication by float.
        if type(other) is float:
            # Convert to Kelvin and multiply.
            mag = self.kelvin().magnitude() * other

            # Get new magnitude in the original unit.
            mag = Temperature(mag, "KELVIN")._convert_to(self._unit).magnitude()

            # Return the new temperature.
            return Temperature(mag, self._unit)

        else:
            raise NotImplementedError

    def __rmul__(self, other):
        """Multiplication operator."""

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Only support multiplication by float.
        if type(other) is float:
            # Convert to Kelvin and multiply.
            mag = self.kelvin().magnitude() * other

            # Get new magnitude in the original unit.
            mag = Temperature(mag, "KELVIN")._convert_to(self._unit).magnitude()

            # Return the new temperature.
            return Temperature(mag, self._unit)

        else:
            raise NotImplementedError

    def __truediv__(self, other):
        """Division operator."""

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Only support division by float.
        if type(other) is float:
            # Convert to Kelvin and divide.
            mag = self.kelvin().magnitude() / other

            # Get new magnitude in the original unit.
            mag = Temperature(mag, "KELVIN")._convert_to(self._unit).magnitude()

            # Return the new temperature.
            return Temperature(mag, self._unit)

        else:
            raise NotImplementedError

    def __lt__(self, other):
        """Less than operator."""
        return self.kelvin().magnitude() < other.kelvin().magnitude()

    def __le__(self, other):
        """Less than or equal to operator."""
        return self.kelvin().magnitude() <= other.kelvin().magnitude()

    def __eq__(self, other):
        """Equals to operator."""
        return self.kelvin().magnitude() == other.kelvin().magnitude()

    def __ne__(self, other):
        """Not equals to operator."""
        return self.kelvin().magnitude() != other.kelvin().magnitude()

    def __ge__(self, other):
        """Greater than or equal to operator."""
        return self.kelvin().magnitude() >= other.kelvin().magnitude()

    def __gt__(self, other):
        """Gretear than operator."""
        return self.kelvin().magnitude() > other.kelvin().magnitude()

    def magnitude(self):
        """Return the magnitude."""
        return self._magnitude

    def unit(self):
        """Return the unit."""
        return self._unit

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
