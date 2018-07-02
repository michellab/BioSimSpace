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

        # Don't support negative lengths.
        if magnitude < 0:
            raise ValueError("The length cannot be negative!")

        # Check that the unit is supported.
        self._unit = self._validate_unit(unit)

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
        if self._magnitude > 1e6 or self._magnitude < 1e-6:
            return "%.4e %s" % (self._magnitude, self._abbrev)
        else:
            return "%.2f %s" % (self._magnitude, self._abbrev)

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._magnitude > 1e6 or self._magnitude < 1e-6:
            return "BioSimSpace.Types.Length(%.4e, '%s')" % (self._magnitude, self._abbrev)
        else:
            return "BioSimSpace.Types.Length(%f, '%s')" % (self._magnitude, self._abbrev)

    def __add__(self, other):
        """Addition operator."""

        # Add the magnitudes in a common unit.
        mag = self.angstroms().magnitude() + other.angstroms().magnitude()

        # Get new magnitude in the original unit.
        # Left-hand operand takes precedence.
        mag = Length(mag, "ANGSTROM")._convert_to(self._unit).magnitude()

        # Return a new length object.
        return Length(mag, self._unit)

    def __sub__(self, other):
        """Subtraction operator."""

        # Subtract the magnitudes in a common unit.
        mag = self.angstroms().magnitude() - other.angstroms().magnitude()

        # Get new magnitude in the original unit.
        # Left-hand operand takes precedence.
        mag = Length(mag, "ANGSTROM")._convert_to(self._unit).magnitude()

        # Return a new length object.
        return Length(mag, self._unit)

    def __mul__(self, other):
        """Multiplication operator."""

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Only support multiplication by float.
        if type(other) is float:
            mag = self._magnitude * other
            return Length(mag, self._unit)

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
            return Length(mag, self._unit)

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
            return Length(mag, self._unit)

        else:
            raise NotImplementedError

    def __lt__(self, other):
        """Less than operator."""
        return self.angstroms().magnitude() < other.angstroms().magnitude()

    def __le__(self, other):
        """Less than or equal to operator."""
        return self.angstroms().magnitude() <= other.angstroms().magnitude()

    def __eq__(self, other):
        """Equals to operator."""
        return self.angstroms().magnitude() == other.angstroms().magnitude()

    def __ne__(self, other):
        """Not equals to operator."""
        return self.angstroms().magnitude() != other.angstroms().magnitude()

    def __ge__(self, other):
        """Greater than or equal to operator."""
        return self.angstroms().magnitude() >= other.angstroms().magnitude()

    def __gt__(self, other):
        """Gretear than operator."""
        return self.angstroms().magnitude() > other.angstroms().magnitude()

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
