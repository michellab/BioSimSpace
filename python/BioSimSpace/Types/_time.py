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
    _abbreviations  = { "HR"  : "HOUR",
                        "MIN" : "MINUTE",
                        "SEC" : "SECOND",
                        "MS"  : "MILLISECOND",
                        "NS"  : "NANOSECOND",
                        "PS"  : "PICOSECOND",
                        "FS"  : "FEMTOSECOND" }

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

        # Don't support negative times.
        if magnitude < 0:
            raise ValueError("The time cannot be negative!")

        # Check that the unit is supported.
        self._unit = self._validate_unit(unit)

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
        if self._magnitude > 1e6:
            return "%.4e %s" % (self._magnitude, self._abbrev)
        else:
            return "%.2f %s" % (self._magnitude, self._abbrev)

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._magnitude > 1e6:
            return "BioSimSpace.Types.Time(%.4e, '%s')" % (self._magnitude, self._abbrev)
        else:
            return "BioSimSpace.Types.Time(%f, '%s')" % (self._magnitude, self._abbrev)

    def __add__(self, other):
        """Addition operator."""

        # Add the magnitudes in a common unit.
        mag = self.picoseconds().magnitude() + other.picoseconds().magnitude()

        # Get new magnitude in the original unit.
        # Left-hand operand takes precedence.
        mag = Time(mag, "PICOSECOND")._convert_to(self._unit).magnitude()

        # Return a new time object.
        return Time(mag, self._unit)

    def __sub__(self, other):
        """Subtraction operator."""

        # Subtract the magnitudes in a common unit.
        mag = self.picoseconds().magnitude() - other.picoseconds().magnitude()

        # Get new magnitude in the original unit.
        # Left-hand operand takes precedence.
        mag = Time(mag, "PICOSECOND")._convert_to(self._unit).magnitude()

        # Return a new time object.
        return Time(mag, self._unit)

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

        # Only support division by float.
        if type(other) is float:
            mag = self._magnitude / other
            return Time(mag, self._unit)

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

    def _convert_to(self, unit):
        """Return the temperature in a different unit.

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
