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

__all__ = ["Pressure"]

class Pressure:
    # Dictionary of allowed units.
    _supported_units = { "ATMOSPHERE" : _Units.atm,
                         "BAR"        : _Units.bar }

    # Map unit abbreviations to the full name.
    _abbreviations = { "ATM" : "ATMOSPHERE",
                       "BAR" : "BAR" }

    def __init__(self, *args):
        """Constructor.

           Positional arguments:

           magnitude -- The magnitude.
           unit      -- The unit.

           or

           string    -- A string representation of the pressure.
        """

        # The user has passed a magnitude and a unit.
        if len(args) > 1:
            magnitude = args[0]
            unit = args[1]

            # Check that the magnitude is valid.
            if type(magnitude) is int:
                self._magnitude = float(magnitude)
            elif type(magnitude) is float:
                self._magnitude = magnitude
            else:
                raise TypeError("'magnitude' must be of type 'int' or 'float'")

            # Check that the unit is supported.
            self._unit = self._validate_unit(unit)

        # The user has passed a string representation of the temperature.
        elif len(args) == 1:
            if type(args[0]) != str:
                raise TypeError("'string' must be of type 'str'")

            # Convert the string to a Pressure object.
            press = self._from_string(args[0])

            # Store the magnitude and unit.
            self._magnitude = press._magnitude
            self._unit = press._unit

        # No arguments.
        else:
            raise TypeError("__init__() missing positional argument(s): 'magnitude' and 'unit', or 'string'")

        # Store the abbreviated unit.
        self._abbrev = list(self._abbreviations.keys())[list(self._abbreviations.values()).index(self._unit)].lower()

    def __str__(self):
        """Return a human readable string representation of the object."""
        if self._magnitude > 1e4 or self._magnitude < 1e-4:
            return "%.4e %s" % (self._magnitude, self._abbrev)
        else:
            return "%5.4f %s" % (self._magnitude, self._abbrev)

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._magnitude > 1e4 or abs(self._magnitude) < 1e-4:
            return "BioSimSpace.Types.Pressure(%.4e, '%s')" % (self._magnitude, self._unit)
        else:
            return "BioSimSpace.Types.Pressure(%5.4f, '%s')" % (self._magnitude, self._unit)

    def __add__(self, other):
        """Addition operator."""

        # Addition of another Pressure object.
        if type(other) is Pressure:
            # Add the magnitudes in a common unit.
            mag = self.atm().magnitude() + other.atm().magnitude()

            # Get new magnitude in the original unit.
            # Left-hand operand takes precedence.
            mag = Pressure(mag, "atm")._convert_to(self._unit).magnitude()

            # Return a new temperature object.
            return Pressure(mag, self._unit)

        # Addition of a string.
        elif type(other) is str:
            temp = self._from_string(other)
            return self + temp

        else:
            raise NotImplementedError

    def __sub__(self, other):
        """Subtraction operator."""

        # Subtraction of another Pressure object.
        if type(other) is Pressure:
            # Subtract the magnitudes in a common unit.
            mag = self.atm().magnitude() - other.atm().magnitude()

            # Get new magnitude in the original unit.
            # Left-hand operand takes precedence.
            mag = Pressure(mag, "atm")._convert_to(self._unit).magnitude()

            # Return a new temperature object.
            return Pressure(mag, self._unit)

        # Addition of a string.
        elif type(other) is str:
            temp = self._from_string(other)
            return self - temp

        else:
            raise NotImplementedError

    def __mul__(self, other):
        """Multiplication operator."""

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Only support multiplication by float.
        if type(other) is float:
            # Convert to atm and multiply.
            mag = self.atm().magnitude() * other

            # Get new magnitude in the original unit.
            mag = Pressure(mag, "atm")._convert_to(self._unit).magnitude()

            # Return the new temperature.
            return Pressure(mag, self._unit)

        else:
            raise NotImplementedError

    def __rmul__(self, other):
        """Multiplication operator."""

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Only support multiplication by float.
        if type(other) is float:
            # Convert to atm and multiply.
            mag = self.atm().magnitude() * other

            # Get new magnitude in the original unit.
            mag = Pressure(mag, "atm")._convert_to(self._unit).magnitude()

            # Return the new temperature.
            return Pressure(mag, self._unit)

        else:
            raise NotImplementedError

    def __truediv__(self, other):
        """Division operator."""

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Float division.
        if type(other) is float:
            # Convert to atm and divide.
            mag = self.atm().magnitude() / other

            # Get new magnitude in the original unit.
            mag = Pressure(mag, "atm")._convert_to(self._unit).magnitude()

            # Return the new temperature.
            return Pressure(mag, self._unit)

        # Division by another temperature.
        elif type(other) is Pressure:
            return self.atm().magnitude() / other.atm().magnitude()

        # Division by a string.
        elif type(other) is str:
            press = self._from_string(other)
            return self / press

        else:
            raise NotImplementedError

    def __lt__(self, other):
        """Less than operator."""

        # Compare to another Pressure object.
        if type(other) is Pressure:
            return self.atm().magnitude() < other.atm().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.atm().magnitude() < self._from_string(other).atm().magnitude()

        else:
            raise NotImplementedError

    def __le__(self, other):
        """Less than or equal to operator."""

        # Compare to another Pressure object.
        if type(other) is Pressure:
            return self.atm().magnitude() <= other.atm().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.atm().magnitude() <= self._from_string(other).atm().magnitude()

        else:
            raise NotImplementedError

    def __eq__(self, other):
        """Equals to operator."""

        # Compare to another Pressure object.
        if type(other) is Pressure:
            return self.atm().magnitude() == other.atm().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.atm().magnitude() == self._from_string(other).atm().magnitude()

        else:
            raise NotImplementedError

    def __ne__(self, other):
        """Not equals to operator."""

        # Compare to another Pressure object.
        if type(other) is Pressure:
            return self.atm().magnitude() != other.atm().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.atm().magnitude() != self._from_string(other).atm().magnitude()

        else:
            raise NotImplementedError

    def __ge__(self, other):
        """Greater than or equal to operator."""

        # Compare to another Pressure object.
        if type(other) is Pressure:
            return self.atm().magnitude() >= other.atm().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.atm().magnitude() >= self._from_string(other).atm().magnitude()

        else:
            raise NotImplementedError

    def __gt__(self, other):
        """Gretear than operator."""

        # Compare to another Pressure object.
        if type(other) is Pressure:
            return self.atm().magnitude() > other.atm().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.atm().magnitude() > self._from_string(other).atm().magnitude()

        else:
            raise NotImplementedError

    def magnitude(self):
        """Return the magnitude."""
        return self._magnitude

    def unit(self):
        """Return the unit."""
        return self._unit

    def atm(self):
        """Return the atmospheric pressure."""
        return Pressure((self._magnitude * self._supported_units[self._unit]).to(_Units.atm), "ATMOSPHERE")

    def bar(self):
        """Return the pressure in bar."""
        return Pressure((self._magnitude * self._supported_units[self._unit]).to(_Units.bar), "BAR")

    def _from_string(self, string):
        """Convert a string to a Pressure object.

           Positional arguments:

           string -- The string to interpret.
        """

        if type(string) is str:
            # Strip white space from the string.
            string = string.replace(" ", "")

            # Try to match scientific format.
            match = _re.search("(\-?\d+\.?\d*e\-?\d+)(.*)", string, _re.IGNORECASE)

            # Try to match decimal format.
            if match is None:
                match = _re.search("(\-?\d+\.?\d*)(.*)", string, _re.IGNORECASE)

                # No matches, raise an error.
                if match is None:
                    raise ValueError("Could not interpret %s: '%s'" % (unit_type, value))

            # Extract the value and unit.
            value, unit = match.groups()

            # Convert the value to a float.
            value = float(value)

            # Create and return a new Pressure object.
            return Pressure(value, unit)

        else:
            raise TypeError("'string' must be of type 'str'")

    def _convert_to(self, unit):
        """Return the temperature in a different unit.

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
