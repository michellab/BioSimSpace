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

__all__ = ["Volume"]

class Volume:
    # Dictionary of allowed units.
    _supported_units = { "METER3"      : _Units.meter2,
                         "NANOMETER3"  : _Units.nanometer2,
                         "ANGSTROM3"   : _Units.angstrom2,
                         "PICOMETER3"  : _Units.picometer2 }

    # Map unit abbreviations to the full name.
    _abbreviations = { "M3"  : "METER3",
                       "NM3" : "NANOMETER3",
                       "A3"  : "ANGSTROM3",
                       "PM3" : "PICOMETER3" }

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
                raise ValueError("The area cannot be negative!")

            # Check that the unit is supported.
            self._unit = self._validate_unit(unit)

        # The user has passed a string representation of the area.
        elif len(args) == 1:
            if type(args[0]) != str:
                raise TypeError("'string' must be of type 'str'")

            # Convert the string to a Volume object.
            volume = self._from_string(args[0])

            # Store the magnitude and unit.
            self._magnitude = volume._magnitude
            self._unit = volume._unit

        # No arguments.
        else:
            raise TypeError("__init__() missing positional argument(s): 'magnitude' and 'unit', or 'string'")

        # Store the abbreviated unit.
        try:
            self._abbrev = list(self._abbreviations.keys())[list(self._abbreviations.values()).index(self._unit)].lower()
            self._abbrev = self._abbrev[0:-1] + "^" + self._abbrev[-1]
        except:
            self._abbrev = self._unit.lower()

        # Handle Angstrom separately.
        if self._abbrev == "a^3":
            self._abbrev = "A^3"

    def __str__(self):
        """Return a human readable string representation of the object."""
        if self._magnitude > 1e6 or self._magnitude < 1e-6:
            return "%.4e %s" % (self._magnitude, self._abbrev)
        else:
            return "%.2f %s" % (self._magnitude, self._abbrev)

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._magnitude > 1e6 or self._magnitude < 1e-6:
            return "BioSimSpace.Types.Volume(%.4e, '%s')" % (self._magnitude, self._abbrev)
        else:
            return "BioSimSpace.Types.Volume(%f, '%s')" % (self._magnitude, self._abbrev)

    def __add__(self, other):
        """Addition operator."""

        # Addition of another Volume object.
        if type(other) is Volume:
            # Add the magnitudes in a common unit.
            mag = self.angstroms3().magnitude() + other.angstroms3().magnitude()

            # Get new magnitude in the original unit.
            # Left-hand operand takes precedence.
            mag = Volume(mag, "ANGSTROM3")._convert_to(self._unit).magnitude()

            # Return a new length object.
            return Volume(mag, self._unit)

        # Addition of a string.
        elif type(other) is str:
            volume = self._from_string(other)
            return self + volume

        else:
            raise NotImplementedError

    def __sub__(self, other):
        """Subtraction operator."""

        # Addition of another Volume object.
        if type(other) is Volume:
            # Subtract the magnitudes in a common unit.
            mag = self.angstroms3().magnitude() - other.angstroms3().magnitude()

            # Get new magnitude in the original unit.
            # Left-hand operand takes precedence.
            mag = Volume(mag, "ANGSTROM3")._convert_to(self._unit).magnitude()

            # Return a new Volume object.
            return Volume(mag, self._unit)

        # Subtraction of a string.
        elif type(other) is str:
            volume = self._from_string(other)
            return self - volume

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
            return Volume(mag, self._unit)

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
            return Volume(mag, self._unit)

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
            return Volume(mag, self._unit)

        # Division by another Volume.
        elif type(other) is Volume:
            return self.angstroms3().magnitude() / other.angstroms3().magnitude()

        # Division by an Area.
        elif type(other) is _Area:
            mag = self.angstroms3().magnitude() / other.angstroms2().magnitude()
            return _Length(mag, "A")

        # Division by a Length.
        elif type(other) is _Length:
            mag = self.angstroms3().magnitude() / other.angstroms().magnitude()
            return _Area(mag, "A2")

        # Division by a string.
        elif type(other) is str:
            try:
                length = _Length(other)
                return self / length
            except:
                try:
                    area = _Area(other)
                    return self / area
                except:
                    try:
                        volume = Volume(other)
                        return self / volume
                    except:
                        raise ValueError("Could not convert the string to a "
                            + "'BioSimSpace.Types.Length', 'BioSimSpace.Types.Area', "
                            + "or 'BioSimSpace.Types.Volume'")
        else:
            raise NotImplementedError

    def __lt__(self, other):
        """Less than operator."""

        # Compare with another Volume object.
        if type(other) is Volume:
            return self.angstroms3().magnitude() < other.angstroms3().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.angstroms3().magnitude() < self._from_string(other).angstroms3().magnitude()

        else:
            raise NotImplementedError

    def __le__(self, other):
        """Less than or equal to operator."""

        # Compare with another Volume object.
        if type(other) is Volume:
            return self.angstroms3().magnitude() <= other.angstroms3().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.angstroms3().magnitude() <= self._from_string(other).angstroms3().magnitude()

        else:
            raise NotImplementedError

    def __eq__(self, other):
        """Equals to operator."""

        # Compare with another Volume object.
        if type(other) is Volume:
            return self.angstroms3().magnitude() == other.angstroms3().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.angstroms3().magnitude() == self._from_string(other).angstroms3().magnitude()

        else:
            raise NotImplementedError

    def __ne__(self, other):
        """Not equals to operator."""

        # Compare with another Volume object.
        if type(other) is Volume:
            return self.angstroms3().magnitude() != other.angstroms3().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.angstroms3().magnitude() != self._from_string(other).angstroms3().magnitude()

        else:
            raise NotImplementedError

    def __ge__(self, other):
        """Greater than or equal to operator."""

        # Compare with another Volume object.
        if type(other) is Volume:
            return self.angstroms3().magnitude() >= other.angstroms3().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.angstroms3().magnitude() >= self._from_string(other).angstroms3().magnitude()

        else:
            raise NotImplementedError

    def __gt__(self, other):
        """Gretear than operator."""

        # Compare with another Volume object.
        if type(other) is Volume:
            return self.angstroms3().magnitude() > other.angstroms3().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.angstroms3().magnitude() > self._from_string(other).angstroms3().magnitude()

        else:
            raise NotImplementedError

    def magnitude(self):
        """Return the magnitude."""
        return self._magnitude

    def unit(self):
        """Return the unit."""
        return self._unit

    def meters3(self):
        """Return the volume in cubic meters."""
        return Volume((self._magnitude * self._supported_units[self._unit]).to(_Units.meter2), "METER3")

    def nanometers3(self):
        """Return the volume in cubic nanometers."""
        return Volume((self._magnitude * self._supported_units[self._unit]).to(_Units.nanometer2), "NANOMETER3")

    def angstroms3(self):
        """Return the volume in cubic angstroms."""
        return Volume((self._magnitude * self._supported_units[self._unit]).to(_Units.angstrom2), "ANGSTROM3")

    def picometers3(self):
        """Return the volume in cubic picometers."""
        return Volume((self._magnitude * self._supported_units[self._unit]).to(_Units.picometer2), "PICOMETER3")

    def _from_string(self, string):
        """Convert a string to a Volume object.

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

            # Create and return a new Volume object.
            return Volume(value, unit)

        else:
            raise TypeError("'string' must be of type 'str'")

    def _convert_to(self, unit):
        """Return the volume in a different unit.

           Positional arguments:

           unit -- The unit to convert to.
        """
        if unit == "METER3":
            return self.meters3()
        elif unit == "NANOMETER3":
            return self.nanometers3()
        elif unit == "ANGSTROM3":
            return self.angstroms3()
        elif unit == "PICOMETER3":
            return self.picometers3()
        else:
            raise ValueError("Supported units are: '%s'" % list(self._supported_units.keys()))

    def _validate_unit(self, unit):
        """Validate that the unit are supported."""

        # Strip whitespace and convert to upper case.
        unit = unit.replace(" ", "").upper()

        # Replace any occurence of cubed with 3.
        unit = unit.replace("CUBED", "3").upper()
        unit = unit.replace("CUBE", "3").upper()

        # Strip "^" character.
        unit = unit.replace("^", "").upper()

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
from ._length import Length as _Length
