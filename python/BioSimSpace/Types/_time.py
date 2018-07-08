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

__all__ = ["Time"]

class Time:
    # Dictionary of allowed units.
    _supported_units = { "DAY"         : _Units.day,
                         "HOUR"        : _Units.hour,
                         "MINUTE"      : _Units.minute,
                         "SECOND"      : _Units.second,
                         "MILLISECOND" : _Units.millisecond,
                         "NANOSECOND"  : _Units.nanosecond,
                         "PICOSECOND"  : _Units.picosecond,
                         "FEMTOSECOND" : _Units.femtosecond }

    # Map unit abbreviations to the full name.
    _abbreviations = { "HR"  : "HOUR",
                       "MIN" : "MINUTE",
                       "SEC" : "SECOND",
                       "MS"  : "MILLISECOND",
                       "NS"  : "NANOSECOND",
                       "PS"  : "PICOSECOND",
                       "FS"  : "FEMTOSECOND" }

    def __init__(self, *args):
        """Constructor.

           Positional arguments:

           magnitude -- The magnitude.
           unit      -- The unit.

           or

           string    -- A string representation of the time.
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

            # Don't support negative times.
            if magnitude < 0:
                raise ValueError("The time cannot be negative!")

            # Check that the unit is supported.
            self._unit = self._validate_unit(unit)

        # The user has passed a string representation of the time.
        elif len(args) == 1:
            if type(args[0]) != str:
                raise TypeError("'string' must be of type 'str'")

            # Convert the string to a Time object.
            time = self._from_string(args[0])

            # Store the magnitude and unit.
            self._magnitude = time._magnitude
            self._unit = time._unit

        # No arguments.
        else:
            raise TypeError("__init__() missing positional argument(s): 'magnitude' and 'unit', or 'string'")

        # Store the abbreviated unit.
        try:
            self._abbrev = list(self._abbreviations.keys())[list(self._abbreviations.values()).index(self._unit)].lower()
        except:
            self._abbrev = self._unit.lower()
        if self._magnitude != 1:
            if self._abbrev[-1] != "s":
                self._abbrev += "s"

    def __str__(self):
        """Return a human readable string representation of the object."""
        if self._magnitude > 1e4 or self._magnitude < 1e-4:
            return "%.4e %s" % (self._magnitude, self._abbrev)
        else:
            return "%5.4f %s" % (self._magnitude, self._abbrev)

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._magnitude > 1e4 or self._magnitude < 1e-4:
            return "BioSimSpace.Types.Time(%.4e, '%s')" % (self._magnitude, self._abbrev)
        else:
            return "BioSimSpace.Types.Time(%5.4f, '%s')" % (self._magnitude, self._abbrev)

    def __add__(self, other):
        """Addition operator."""

        # Addition of another Time object.
        if type(other) is Time:
            # Add the magnitudes in a common unit.
            mag = self.picoseconds().magnitude() + other.picoseconds().magnitude()

            # Get new magnitude in the original unit.
            # Left-hand operand takes precedence.
            mag = Time(mag, "PICOSECOND")._convert_to(self._unit).magnitude()

            # Return a new time object.
            return Time(mag, self._unit)

        # Addition of a string.
        elif type(other) is str:
            time = self._from_string(other)
            return self + time

        else:
            raise NotImplementedError

    def __sub__(self, other):
        """Subtraction operator."""

        # Subtraction of another Time object.
        if type(other) is Time:
            # Subtract the magnitudes in a common unit.
            mag = self.picoseconds().magnitude() - other.picoseconds().magnitude()

            # Get new magnitude in the original unit.
            # Left-hand operand takes precedence.
            mag = Time(mag, "PICOSECOND")._convert_to(self._unit).magnitude()

            # Return a new time object.
            return Time(mag, self._unit)

        # Subtraction of a string.
        elif type(other) is str:
            time = self._from_string(other)
            return self - time

        else:
            raise NotImplementedError

    def __mul__(self, other):
        """Multiplication operator."""

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Only support multiplication by float.
        if type(other) is float:
            mag = self._magnitude * other
            return Time(mag, self._unit)

        else:
            raise NotImplementedError

    def __rmul__(self, other):
        """Multiplication operator."""

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Only support multiplication by float.
        if type(other) is float:
            mag = self._magnitude * other
            return Time(mag, self._unit)

        else:
            raise NotImplementedError

    def __truediv__(self, other):
        """Division operator."""

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Float division.
        if type(other) is float:
            mag = self._magnitude / other
            return Time(mag, self._unit)

        # Divide by time.
        elif type(other) is Time:
            return self.picoseconds().magnitude() / other.picoseconds().magnitude()

        # Divide by string.
        elif type(other) is str:
            time = self._from_string(other)
            return self / time

        else:
            raise NotImplementedError

    def __lt__(self, other):
        """Less than operator."""

        # Compare to another Time object.
        if type(other) is Time:
            return self.picoseconds().magnitude() < other.picoseconds().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.picoseconds().magnitude() < self._from_string(other).picoseconds().magnitude()

        else:
            raise NotImplementedError

    def __le__(self, other):
        """Less than or equal to operator."""

        # Compare to another Time object.
        if type(other) is Time:
            return self.picoseconds().magnitude() <= other.picoseconds().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.picoseconds().magnitude() <= self._from_string(other).picoseconds().magnitude()

        else:
            raise NotImplementedError

    def __eq__(self, other):
        """Equals to operator."""

        # Compare to another Time object.
        if type(other) is Time:
            return self.picoseconds().magnitude() == other.picoseconds().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.picoseconds().magnitude() == self._from_string(other).picoseconds().magnitude()

        else:
            raise NotImplementedError

    def __ne__(self, other):
        """Not equals to operator."""

        # Compare to another Time object.
        if type(other) is Time:
            return self.picoseconds().magnitude() != other.picoseconds().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.picoseconds().magnitude() != self._from_string(other).picoseconds().magnitude()

        else:
            raise NotImplementedError

    def __ge__(self, other):
        """Greater than or equal to operator."""

        # Compare to another Time object.
        if type(other) is Time:
            return self.picoseconds().magnitude() >= other.picoseconds().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.picoseconds().magnitude() >= self._from_string(other).picoseconds().magnitude()

        else:
            raise NotImplementedError

    def __gt__(self, other):
        """Gretear than operator."""

        # Compare to another Time object.
        if type(other) is Time:
            return self.picoseconds().magnitude() > other.picoseconds().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.picoseconds().magnitude() > self._from_string(other).picoseconds().magnitude()

        else:
            raise NotImplementedError

    def magnitude(self):
        """Return the magnitude."""
        return self._magnitude

    def unit(self):
        """Return the unit."""
        return self._unit

    def weeks(self):
        """Return the time in weeks."""
        return Time((self._magnitude * self._supported_units[self._unit]).to(_Units.week), "WEEK")

    def days(self):
        """Return the time in days."""
        return Time((self._magnitude * self._supported_units[self._unit]).to(_Units.day), "DAY")

    def hours(self):
        """Return the time in hours."""
        return Time((self._magnitude * self._supported_units[self._unit]).to(_Units.hour), "HOUR")

    def minutes(self):
        """Return the time in minutes."""
        return Time((self._magnitude * self._supported_units[self._unit]).to(_Units.minute), "MINUTE")

    def seconds(self):
        """Return the time in seconds."""
        return Time((self._magnitude * self._supported_units[self._unit]).to(_Units.second), "SECOND")

    def milliseconds(self):
        """Return the time in milliseconds."""
        return Time((self._magnitude * self._supported_units[self._unit]).to(_Units.millisecond), "MILLISECOND")

    def nanoseconds(self):
        """Return the time in nanoseconds."""
        return Time((self._magnitude * self._supported_units[self._unit]).to(_Units.nanosecond), "NANOSECOND")

    def picoseconds(self):
        """Return the time in picoseconds."""
        return Time((self._magnitude * self._supported_units[self._unit]).to(_Units.picosecond), "PICOSECOND")

    def femtoseconds(self):
        """Return the time in femtoseconds."""
        return Time((self._magnitude * self._supported_units[self._unit]).to(_Units.femtosecond), "FEMTOSECOND")

    def _from_string(self, string):
        """Convert a string to a Time object.

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

            # Create and return a new Time object.
            return Time(value, unit)

        else:
            raise TypeError("'string' must be of type 'str'")

    def _convert_to(self, unit):
        """Return the time in a different unit.

           Positional arguments:

           unit -- The unit to convert to.
        """
        if unit == "WEEK":
            return self.weeks()
        elif unit == "DAY":
            return self.days()
        elif unit == "HOUR":
            return self.hours()
        elif unit == "MINUTE":
            return self.minutes()
        elif unit == "SECOND":
            return self.seconds()
        elif unit == "MILLISECOND":
            return self.milliseconds()
        elif unit == "NANOSECOND":
            return self.nanoseconds()
        elif unit == "PICOSECOND":
            return self.picoseconds()
        elif unit == "FEMTOSECOND":
            return self.femtoseconds()
        else:
            raise ValueError("Supported units are: '%s'" % list(self._supported_units.keys()))

    def _validate_unit(self, unit):
        """Validate that the unit are supported."""

        # Strip whitespace and convert to upper case.
        unit = unit.replace(" ", "").upper()

        # Check that the unit is supported.
        if unit in self._supported_units:
            return unit
        elif unit[:-1] in self._supported_units:
            return unit[:-1]
        elif unit in self._abbreviations:
            return self._abbreviations[unit]
        elif unit[:-1] in self._abbreviations:
            return self._abbreviations[unit[:-1]]
        else:
            raise ValueError("Supported units are: '%s'" % list(self._supported_units.keys()))
