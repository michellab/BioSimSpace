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

__all__ = ["Length"]

class Length:
    # Dictionary of allowed units.
    _supported_units = { "METER"      : _Units.meter,
                         "CENTIMETER" : _Units.centimeter,
                         "MILLIMETER" : _Units.millimeter,
                         "NANOMETER"  : _Units.nanometer,
                         "ANGSTROM"   : _Units.angstrom,
                         "PICOMETER"  : _Units.picometer }

    # Map unit abbreviations to the full name.
    _abbreviations = { "M"  : "METER",
                       "CM" : "CENTIMETER",
                       "MM" : "MILLIMETER",
                       "NM" : "NANOMETER",
                       "A"  : "ANGSTROM",
                       "PM" : "PICOMETER" }

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
                raise ValueError("The length cannot be negative!")

            # Check that the unit is supported.
            self._unit = self._validate_unit(unit)

        # The user has passed a string representation of the length.
        elif len(args) == 1:
            if type(args[0]) != str:
                raise TypeError("'string' must be of type 'str'")

            # Convert the string to a Length object.
            length = self._from_string(args[0])

            # Store the magnitude and unit.
            self._magnitude = length._magnitude
            self._unit = length._unit

        # No arguments.
        else:
            raise TypeError("__init__() missing positional argument(s): 'magnitude' and 'unit', or 'string'")

        # Store the abbreviated unit.
        try:
            self._abbrev = list(self._abbreviations.keys())[list(self._abbreviations.values()).index(self._unit)].lower()
        except:
            self._abbrev = self._unit.lower()

        # Handle Angstrom separately.
        if self._abbrev == "a":
            self._abbrev = "A"

    def __str__(self):
        """Return a human readable string representation of the object."""
        if self._magnitude > 1e4 or self._magnitude < 1e-4:
            return "%.4e %s" % (self._magnitude, self._abbrev)
        else:
            return "%5.4f %s" % (self._magnitude, self._abbrev)

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._magnitude > 1e4 or self._magnitude < 1e-4:
            return "BioSimSpace.Types.Length(%.4e, '%s')" % (self._magnitude, self._abbrev)
        else:
            return "BioSimSpace.Types.Length(%5.4f, '%s')" % (self._magnitude, self._abbrev)

    def __add__(self, other):
        """Addition operator."""

        # Addition of another Length object.
        if type(other) is Length:
            # Add the magnitudes in a common unit.
            mag = self.angstroms().magnitude() + other.angstroms().magnitude()

            # Get new magnitude in the original unit.
            # Left-hand operand takes precedence.
            mag = Length(mag, "ANGSTROM")._convert_to(self._unit).magnitude()

            # Return a new length object.
            return Length(mag, self._unit)

        # Addition of a string.
        elif type(other) is str:
            length = self._from_string(other)
            return self + length

        else:
            raise NotImplementedError

    def __sub__(self, other):
        """Subtraction operator."""

        # Subtraction of another Length object.
        if type(other) is Length:
            # Subtract the magnitudes in a common unit.
            mag = self.angstroms().magnitude() - other.angstroms().magnitude()

            # Get new magnitude in the original unit.
            # Left-hand operand takes precedence.
            mag = Length(mag, "ANGSTROM")._convert_to(self._unit).magnitude()

            # Return a new length object.
            return Length(mag, self._unit)

        # Addition of a string.
        elif type(other) is str:
            length = self._from_string(other)
            return self - length

        else:
            raise NotImplementedError

    def __mul__(self, other):
        """Multiplication operator."""

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Multiplication by float.
        if type(other) is float:
            mag = self._magnitude * other
            return Length(mag, self._unit)

        # Multiplication by another Length.
        elif type(other) is Length:
            mag = self.angstroms().magnitude() * other.angstroms().magnitude()
            return _Area(mag, "A2")

        # Multiplication by an Area.
        elif type(other) is _Area:
            mag = self.angstroms().magnitude() * other.angstroms2().magnitude()
            return _Volume(mag, "A3")

        # Multiplication by a string.
        elif type(other) is str:
            try:
                length = Length(other)
                return self * length
            except:
                try:
                    area = _Area(other)
                    return self * area
                except:
                    raise ValueError("Could not convert the string to a 'BioSimSpace.Length' "
                        + "or 'BioSimSpace.Area' type.")
        else:
            raise NotImplementedError

    def __rmul__(self, other):
        """Multiplication operator."""

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Multiplication by float.
        if type(other) is float:
            mag = self._magnitude * other
            return Length(mag, self._unit)

        # Multiplication by another Length.
        elif type(other) is Length:
            mag = self.angstroms().magnitude() * other.angstroms().magnitude()
            return _Area(mag, "A2")

        # Multiplication by an Area.
        elif type(other) is _Area:
            mag = self.angstroms().magnitude() * other.angstroms2().magnitude()
            return _Volume(mag, "A3")

        # Multiplication by a string.
        elif type(other) is str:
            try:
                length = Length(other)
                return self * length
            except:
                try:
                    area = _Area(other)
                    return self * area
                except:
                    raise ValueError("Could not convert the string to a "
                        + "'BioSimSpace.Types.Length' or a 'BioSimSpace.Types.Area'.")
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
            return Length(mag, self._unit)

        # Division by another Length.
        elif type(other) is Length:
            return self.angstroms().magnitude() / other.angstroms().magnitude()

        # Division by a string.
        elif type(other) is str:
            try:
                length = self._from_string(other)
                return self / length
            except:
                raise ValueError("Could not convert the string to a 'BioSimSpace.Types.Length'.")
        else:
            raise NotImplementedError

    def __lt__(self, other):
        """Less than operator."""

        # Compare with another Length object.
        if type(other) is Length:
            return self.angstroms().magnitude() < other.angstroms().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.angstroms().magnitude() < self._from_string(other).angstroms().magnitude()

        else:
            raise NotImplementedError

    def __le__(self, other):
        """Less than or equal to operator."""

        # Compare with another Length object.
        if type(other) is Length:
            return self.angstroms().magnitude() <= other.angstroms().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.angstroms().magnitude() <= self._from_string(other).angstroms().magnitude()

        else:
            raise NotImplementedError

    def __eq__(self, other):
        """Equals to operator."""

        # Compare with another Length object.
        if type(other) is Length:
            return self.angstroms().magnitude() == other.angstroms().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.angstroms().magnitude() == self._from_string(other).angstroms().magnitude()

        else:
            raise NotImplementedError

    def __ne__(self, other):
        """Not equals to operator."""

        # Compare with another Length object.
        if type(other) is Length:
            return self.angstroms().magnitude() != other.angstroms().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.angstroms().magnitude() != self._from_string(other).angstroms().magnitude()

        else:
            raise NotImplementedError

    def __ge__(self, other):
        """Greater than or equal to operator."""

        # Compare with another Length object.
        if type(other) is Length:
            return self.angstroms().magnitude() >= other.angstroms().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.angstroms().magnitude() >= self._from_string(other).angstroms().magnitude()

        else:
            raise NotImplementedError

    def __gt__(self, other):
        """Gretear than operator."""

        # Compare with another Length object.
        if type(other) is Length:
            return self.angstroms().magnitude() > other.angstroms().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.angstroms().magnitude() > self._from_string(other).angstroms().magnitude()

        else:
            raise NotImplementedError

    def magnitude(self):
        """Return the magnitude."""
        return self._magnitude

    def unit(self):
        """Return the unit."""
        return self._unit

    def meters(self):
        """Return the length in meters."""
        return Length((self._magnitude * self._supported_units[self._unit]).to(_Units.meter), "METER")

    def centimeters(self):
        """Return the length in centimeters."""
        return Length((self._magnitude * self._supported_units[self._unit]).to(_Units.centimeter), "CENTIMETER")

    def millimeters(self):
        """Return the length in millimeters."""
        return Length((self._magnitude * self._supported_units[self._unit]).to(_Units.millimeter), "MILLIMETER")

    def nanometers(self):
        """Return the length in nanometers."""
        return Length((self._magnitude * self._supported_units[self._unit]).to(_Units.nanometer), "NANOMETER")

    def angstroms(self):
        """Return the length in angstroms."""
        return Length((self._magnitude * self._supported_units[self._unit]).to(_Units.angstrom), "ANGSTROM")

    def picometers(self):
        """Return the length in picometers."""
        return Length((self._magnitude * self._supported_units[self._unit]).to(_Units.picometer), "PICOMETER")

    def _from_string(self, string):
        """Convert a string to a Length object.

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

            # Create and return a new Length object.
            return Length(value, unit)

        else:
            raise TypeError("'string' must be of type 'str'")

    def _convert_to(self, unit):
        """Return the length in a different unit.

           Positional arguments:

           unit -- The unit to convert to.
        """
        if unit == "METER":
            return self.meters()
        elif unit == "CENTIMETER":
            return self.centimeters()
        elif unit == "MILLIMETER":
            return self.millimeters()
        elif unit == "NANOMETER":
            return self.nanometers()
        elif unit == "ANGSTROM":
            return self.angstroms()
        elif unit == "PICOMETER":
            return self.picometers()
        else:
            raise ValueError("Supported units are: '%s'" % list(self._supported_units.keys()))

    def _validate_unit(self, unit):
        """Validate that the unit are supported."""

        # Strip whitespace and convert to upper case.
        unit = unit.replace(" ", "").upper()

        # Strip any "S" characters.
        unit = unit.replace("S", "").upper()

        # Fix for ANGSTROM (since it contains an "S").
        if unit[0:3] == "ANG":
            unit = "ANGS" + unit[3:]

        # Check that the unit is supported.
        if unit in self._supported_units:
            return unit
        elif unit in self._abbreviations:
            return self._abbreviations[unit]
        else:
            raise ValueError("Supported units are: '%s'" % list(self._supported_units.keys()))

# Import at bottom of module to avoid circular dependency.
from ._area import Area as _Area
from ._volume import Volume as _Volume
