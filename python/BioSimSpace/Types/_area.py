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

__all__ = ["Area"]

class Area:
    # Dictionary of allowed units.
    _supported_units = { "METER2"      : _Units.meter2,
                         "NANOMETER2"  : _Units.nanometer2,
                         "ANGSTROM2"   : _Units.angstrom2,
                         "PICOMETER2"  : _Units.picometer2 }

    # Map unit abbreviations to the full name.
    _abbreviations = { "M^2"  : "METER2",
                       "NM^2" : "NANOMETER2",
                       "A^2"  : "ANGSTROM2",
                       "PM^2" : "PICOMETER2" }

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

        # Don't support negative areas.
        if magnitude < 0:
            raise ValueError("The area cannot be negative!")

        # Check that the unit is supported.
        self._unit = self._validate_unit(unit)

        # Store the abbreviated unit.
        try:
            self._abbrev = list(self._abbreviations.keys())[list(self._abbreviations.values()).index(self._unit)].lower()
        except:
            self._abbrev = self._unit.lower()

        # Handle Angstrom separately.
        if self._abbrev == "a^2":
            self._abbrev = "A^2"

    def __str__(self):
        """Return a human readable string representation of the object."""
        if self._magnitude > 1e6 or self._magnitude < 1e-6:
            return "%.4e %s" % (self._magnitude, self._abbrev)
        else:
            return "%.2f %s" % (self._magnitude, self._abbrev)

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._magnitude > 1e6 or self._magnitude < 1e-6:
            return "BioSimSpace.Types.Area(%.4e, '%s')" % (self._magnitude, self._abbrev)
        else:
            return "BioSimSpace.Types.Area(%f, '%s')" % (self._magnitude, self._abbrev)

    def __add__(self, other):
        """Addition operator."""

        # Add the magnitudes in a common unit.
        mag = self.angstroms2().magnitude() + other.angstroms2().magnitude()

        # Get new magnitude in the original unit.
        # Left-hand operand takes precedence.
        mag = Area(mag, "ANGSTROM2")._convert_to(self._unit).magnitude()

        # Return a new length object.
        return Area(mag, self._unit)

    def __sub__(self, other):
        """Subtraction operator."""

        # Subtract the magnitudes in a common unit.
        mag = self.angstroms2().magnitude() - other.angstroms2().magnitude()

        # Get new magnitude in the original unit.
        # Left-hand operand takes precedence.
        mag = Area(mag, "ANGSTROM2")._convert_to(self._unit).magnitude()

        # Return a new Area object.
        return Area(mag, self._unit)

    def __mul__(self, other):
        """Multiplication operator."""

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Multiplication by float.
        if type(other) is float:
            mag = self._magnitude * other
            return Area(mag, self._unit)

        # Multiplication by a Length.
        elif type(other) is _Length:
            mag = self.angstroms2().magnitude() * other.angstrom().magnitude()
            return _Volume(mag, "A3")

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
            return Area(mag, self._unit)

        # Multiplication by a Length.
        elif type(other) is _Length:
            mag = self.angstroms2().magnitude() * other.angstrom().magnitude()
            return _Volume(mag, "A3")

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
            return Area(mag, self._unit)

        # Division by another Area.
        elif type(other) is Area:
            return self.angstroms2().magnitude() / other.angstroms2().magnitude()

        # Division by a Length.
        elif type(other) is _Length:
            mag = self.angstroms2().magnitude() / other.angstroms().magnitude()
            return _Length(mag, "A")

        else:
            raise NotImplementedError

    def __lt__(self, other):
        """Less than operator."""
        return self.angstroms2().magnitude() < other.angstroms2().magnitude()

    def __le__(self, other):
        """Less than or equal to operator."""
        return self.angstroms2().magnitude() <= other.angstroms2().magnitude()

    def __eq__(self, other):
        """Equals to operator."""
        return self.angstroms2().magnitude() == other.angstroms2().magnitude()

    def __ne__(self, other):
        """Not equals to operator."""
        return self.angstroms2().magnitude() != other.angstroms2().magnitude()

    def __ge__(self, other):
        """Greater than or equal to operator."""
        return self.angstroms2().magnitude() >= other.angstroms2().magnitude()

    def __gt__(self, other):
        """Gretear than operator."""
        return self.angstroms2().magnitude() > other.angstroms2().magnitude()

    def magnitude(self):
        """Return the magnitude."""
        return self._magnitude

    def unit(self):
        """Return the unit."""
        return self._unit

    def meters2(self):
        """Return the area in square meters."""
        return Area((self._magnitude * self._supported_units[self._unit]).to(_Units.meter2), "METER2")

    def nanometers2(self):
        """Return the area in square nanometers."""
        return Area((self._magnitude * self._supported_units[self._unit]).to(_Units.nanometer2), "NANOMETER2")

    def angstroms2(self):
        """Return the area in square angstroms."""
        return Area((self._magnitude * self._supported_units[self._unit]).to(_Units.angstrom2), "ANGSTROM2")

    def picometers2(self):
        """Return the area in square picometers."""
        return Area((self._magnitude * self._supported_units[self._unit]).to(_Units.picometer2), "PICOMETER2")

    def _convert_to(self, unit):
        """Return the area in a different unit.

           Positional arguments:

           unit -- The unit to convert to.
        """
        if unit == "METER2":
            return self.meters2()
        elif unit == "NANOMETER2":
            return self.nanometers2()
        elif unit == "ANGSTROM2":
            return self.angstroms2()
        elif unit == "PICOMETER2":
            return self.picometers2()
        else:
            raise ValueError("Supported units are: '%s'" % list(self._supported_units.keys()))

    def _validate_unit(self, unit):
        """Validate that the unit are supported."""

        # Strip whitespace and convert to upper case.
        unit = unit.replace(" ", "").upper()

        # Replace any occurence of squared with 2.
        unit = unit.replace("SQUARED", "2").upper()
        unit = unit.replace("SQUARE", "2").upper()

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
from ._length import Length as _Length
from ._volume import Volume as _Volume
