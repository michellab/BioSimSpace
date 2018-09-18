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
An area type.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

import Sire.Units as _Units

from ._type import Type as _Type

__all__ = ["Area"]

class Area(_Type):
    # Dictionary of allowed units.
    _supported_units = { "METER2"      : _Units.meter2,
                         "NANOMETER2"  : _Units.nanometer2,
                         "ANGSTROM2"   : _Units.angstrom2,
                         "PICOMETER2"  : _Units.picometer2 }

    # Map unit abbreviations to the full name.
    _abbreviations = { "M2"  : "METER2",
                       "NM2" : "NANOMETER2",
                       "A2"  : "ANGSTROM2",
                       "PM2" : "PICOMETER2" }

    # Print format.
    _print_format = { "METER2"      : "m^2",
                      "NANOMETER2"  : "nm^2",
                      "ANGSTROM2"   : "A^2",
                      "PICOMETER2"  : "pm^2" }

    def __init__(self, *args):
        """Constructor.

           Positional arguments
           --------------------

           magnitude : float
               The magnitude.

           unit : str
               The unit.

           or

           string : str
               A string representation of the area.
        """

        # Call the base class constructor.
        super().__init__(*args)

        # Don't support negative areas.
        if self._magnitude < 0:
            raise ValueError("The area cannot be negative!")

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
            mag = self.angstroms2().magnitude() * other.angstroms().magnitude()
            return _Volume(mag, "A3")

        # Multiplication by a string.
        elif type(other) is str:
            try:
                length = _Length(other)
                return self * length
            except:
                raise ValueError("Could not convert the string to a 'BioSimSpace.Types.Length'")

        else:
            raise TypeError("unsupported operand type(s) for *: '%s' and '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__))

    def __rmul__(self, other):
        """Multiplication operator."""

        # Multipliation is commutative: a*b = b*a
        return self.__mul__(other)

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

        # Division by a string.
        elif type(other) is str:
            try:
                length = _Length(other)
                return self / length
            except:
                try:
                    area = Area(other)
                    return self / area
                except:
                    raise ValueError("Could not convert the string to a "
                                     "'BioSimSpace.Types.Length' or a 'BioSimSpace.Types.Area'.")

        else:
            raise TypeError("unsupported operand type(s) for /: '%s' and '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__))

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

    def _default_unit(self, mag=None):
        """Internal method to return an object of the same type in the default unit.

           Positional argument:

           mag -- The magnitude (optional).
        """
        if mag is None:
            return self.angstroms2()
        else:
            return Area(mag, "ANGSTROM2")

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
