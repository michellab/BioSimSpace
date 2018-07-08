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

__all__ = ["Energy"]

class Energy:
    # Dictionary of allowed units.
    _supported_units = { "KILO CALORIES PER MOL" : _Units.kcal_per_mol,
                         "KILO JOULES PER MOL"   : _Units.kJ_per_mol,
                         "KT"                    : 2.479 * _Units.kJ_per_mol }

    # Map unit abbreviations to the full name.
    _abbreviations = { "KCAL/MOL" : "KILO CALORIES PER MOL",
                       "KJ/MOL"   : "KILO JOULES PER MOL",
                       "KT"       : "KT" }

    # Print formatting.
    _print_format = { "KILO CALORIES PER MOL" : "kcal/mol",
                      "KILO JOULES PER MOL"   : "kJ/mol",
                      "KT"                    : "KT" }

    def __init__(self, *args):
        """Constructor.

           Positional arguments:

           magnitude -- The magnitude.
           unit      -- The unit.

           or

           string    -- A string representation of the energy.
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

            # Convert the string to a Energy object.
            nrg = self._from_string(args[0])

            # Store the magnitude and unit.
            self._magnitude = nrg._magnitude
            self._unit = nrg._unit

        # No arguments.
        else:
            raise TypeError("__init__() missing positional argument(s): 'magnitude' and 'unit', or 'string'")

    def __str__(self):
        """Return a human readable string representation of the object."""
        if self._magnitude > 1e4 or self._magnitude < 1e-4:
            return "%.4e %s" % (self._magnitude, self._print_format[self.unit])
        else:
            return "%5.4f %s" % (self._magnitude, self._print_format[self._unit])

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._magnitude > 1e4 or abs(self._magnitude) < 1e-4:
            return "BioSimSpace.Types.Energy(%.4e, '%s')" % (self._magnitude, self._unit)
        else:
            return "BioSimSpace.Types.Energy(%5.4f, '%s')" % (self._magnitude, self._unit)

    def __add__(self, other):
        """Addition operator."""

        # Addition of another Energy object.
        if type(other) is Energy:
            # Add the magnitudes in a common unit.
            mag = self.kcal_per_mol().magnitude() + other.kcal_per_mol().magnitude()

            # Get new magnitude in the original unit.
            # Left-hand operand takes precedence.
            mag = Energy(mag, "kcal_per_mol")._convert_to(self._unit).magnitude()

            # Return a new temperature object.
            return Energy(mag, self._unit)

        # Addition of a string.
        elif type(other) is str:
            temp = self._from_string(other)
            return self + temp

        else:
            raise NotImplementedError

    def __sub__(self, other):
        """Subtraction operator."""

        # Subtraction of another Energy object.
        if type(other) is Energy:
            # Subtract the magnitudes in a common unit.
            mag = self.kcal_per_mol().magnitude() - other.kcal_per_mol().magnitude()

            # Get new magnitude in the original unit.
            # Left-hand operand takes precedence.
            mag = Energy(mag, "kcal_per_mol")._convert_to(self._unit).magnitude()

            # Return a new temperature object.
            return Energy(mag, self._unit)

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
            # Convert to kcal_per_mol and multiply.
            mag = self.kcal_per_mol().magnitude() * other

            # Get new magnitude in the original unit.
            mag = Energy(mag, "kcal_per_mol")._convert_to(self._unit).magnitude()

            # Return the new temperature.
            return Energy(mag, self._unit)

        else:
            raise NotImplementedError

    def __rmul__(self, other):
        """Multiplication operator."""

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Only support multiplication by float.
        if type(other) is float:
            # Convert to kcal_per_mol and multiply.
            mag = self.kcal_per_mol().magnitude() * other

            # Get new magnitude in the original unit.
            mag = Energy(mag, "kcal_per_mol")._convert_to(self._unit).magnitude()

            # Return the new temperature.
            return Energy(mag, self._unit)

        else:
            raise NotImplementedError

    def __truediv__(self, other):
        """Division operator."""

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Float division.
        if type(other) is float:
            # Convert to kcal_per_mol and divide.
            mag = self.kcal_per_mol().magnitude() / other

            # Get new magnitude in the original unit.
            mag = Energy(mag, "kcal_per_mol")._convert_to(self._unit).magnitude()

            # Return the new temperature.
            return Energy(mag, self._unit)

        # Division by another temperature.
        elif type(other) is Energy:
            return self.kcal_per_mol().magnitude() / other.kcal_per_mol().magnitude()

        # Division by a string.
        elif type(other) is str:
            temp = self._from_string(other)
            return self / temp

        else:
            raise NotImplementedError

    def __lt__(self, other):
        """Less than operator."""

        # Compare to another Energy object.
        if type(other) is Energy:
            return self.kcal_per_mol().magnitude() < other.kcal_per_mol().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.kcal_per_mol().magnitude() < self._from_string(other).kcal_per_mol().magnitude()

        else:
            raise NotImplementedError

    def __le__(self, other):
        """Less than or equal to operator."""

        # Compare to another Energy object.
        if type(other) is Energy:
            return self.kcal_per_mol().magnitude() <= other.kcal_per_mol().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.kcal_per_mol().magnitude() <= self._from_string(other).kcal_per_mol().magnitude()

        else:
            raise NotImplementedError

    def __eq__(self, other):
        """Equals to operator."""

        # Compare to another Energy object.
        if type(other) is Energy:
            return self.kcal_per_mol().magnitude() == other.kcal_per_mol().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.kcal_per_mol().magnitude() == self._from_string(other).kcal_per_mol().magnitude()

        else:
            raise NotImplementedError

    def __ne__(self, other):
        """Not equals to operator."""

        # Compare to another Energy object.
        if type(other) is Energy:
            return self.kcal_per_mol().magnitude() != other.kcal_per_mol().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.kcal_per_mol().magnitude() != self._from_string(other).kcal_per_mol().magnitude()

        else:
            raise NotImplementedError

    def __ge__(self, other):
        """Greater than or equal to operator."""

        # Compare to another Energy object.
        if type(other) is Energy:
            return self.kcal_per_mol().magnitude() >= other.kcal_per_mol().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.kcal_per_mol().magnitude() >= self._from_string(other).kcal_per_mol().magnitude()

        else:
            raise NotImplementedError

    def __gt__(self, other):
        """Gretear than operator."""

        # Compare to another Energy object.
        if type(other) is Energy:
            return self.kcal_per_mol().magnitude() > other.kcal_per_mol().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self.kcal_per_mol().magnitude() > self._from_string(other).kcal_per_mol().magnitude()

        else:
            raise NotImplementedError

    def magnitude(self):
        """Return the magnitude."""
        return self._magnitude

    def unit(self):
        """Return the unit."""
        return self._unit

    def kcal_per_mol(self):
        """Return the energy in kcal per mol."""
        return Energy((self._magnitude * self._supported_units[self._unit]).to(_Units.kcal_per_mol), "KILO CALORIES PER MOL")

    def kj_per_mol(self):
        """Return the energy in kJ per mol."""
        return Energy((self._magnitude * self._supported_units[self._unit]).to(_Units.kJ_per_mol), "KILO JOULES PER MOL")

    def kt(self):
        """Return the energy in KT."""
        return Energy((self._magnitude * self._supported_units[self._unit]).to(2.479 * _Units.kJ_per_mol), "KT")

    def _from_string(self, string):
        """Convert a string to a Energy object.

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

            # Create and return a new Energy object.
            return Energy(value, unit)

        else:
            raise TypeError("'string' must be of type 'str'")

    def _convert_to(self, unit):
        """Return the temperature in a different unit.

           Positional arguments:

           unit -- The unit to convert to.
        """
        if unit == "KILO CALORIES PER MOL":
            return self.kcal_per_mol()
        elif unit == "KILO JOULES PER MOL":
            return self.kj_per_mol()
        elif unit == "KT":
            return self.kt()
        else:
            raise ValueError("Supported units are: '%s'" % list(self._supported_units.keys()))

    def _validate_unit(self, unit):
        """Validate that the unit are supported."""

        # Strip whitespace and convert to upper case.
        unit = unit.replace(" ", "").upper()

        # Replace all instances of KILO with K.
        unit = unit.replace("KILO", "K")

        # Replace all instances of PER with the / character.
        unit = unit.replace("PER", "/")

        # Replace all instance of MOLE with MOL.
        unit = unit.replace("MOLE", "MOL")

        # Replace all instances of CALORIES with CAL.
        unit = unit.replace("CALORIES", "CAL")

        # Replace all instances of JOULES with J.
        unit = unit.replace("JOULES", "J")

        # Check that the unit is supported.
        if unit in self._abbreviations:
            return self._abbreviations[unit]
        else:
            raise ValueError("Supported units are: '%s'" % list(self._supported_units.keys()))
