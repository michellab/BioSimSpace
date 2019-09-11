######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2019
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
A collection of physical unit types.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["Type"]

import re as _re

class Type():
    """A base class for custom types."""
    def __init__(self, *args):
        """Constructor.

           ``*args`` can be a magnitude and unit, or a string representation
           of the type.

           Parameters
           ----------

           magnitude : float
               The magnitude.

           unit : str
               The unit.

           string : str
               A string representation of the type.
        """

	# Don't allow user to create an instance of this base class.
        if type(self) is Type:
            raise Exception("<Type> must be subclassed.")

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

            # Convert the string to an object of this type.
            obj = self._from_string(args[0])

            # Store the magnitude and unit.
            self._magnitude = obj._magnitude
            self._unit = obj._unit

        # No arguments.
        else:
            raise TypeError("__init__() missing positional argument(s): 'magnitude' and 'unit', or 'string'")

        # Set the documentation string.
        self.__doc__ = self._doc_strings[self._unit]

    def __str__(self):
        """Return a human readable string representation of the object."""
        if abs(self._magnitude) > 1e4 or abs(self._magnitude) < 1e-4:
            return "%.4e %s" % (self._magnitude, self._print_format[self._unit])
        else:
            return "%5.4f %s" % (self._magnitude, self._print_format[self._unit])

    def __repr__(self):
        """Return a human readable string representation of the object."""
        return self.__str__()

    def __pos__(self):
        """Unary + operator."""
        return type(self)(self.magnitude(), self.unit())

    def __neg__(self):
        """Unary - operator."""
        return type(self)(-self.magnitude(), self.unit())

    def __add__(self, other):
        """Addition operator."""

        # Addition of another object of the same type.
        if type(other) is type(self):
            # Add the magnitudes in a common unit.
            mag = self._default_unit().magnitude() + other._default_unit().magnitude()

            # Return a new object of the same type with the original unit.
            return self._default_unit(mag)._convert_to(self._unit)

        # Addition of a string.
        elif type(other) is str:
            temp = self._from_string(other)
            return self + temp

        else:
            raise TypeError("unsupported operand type(s) for +: '%s' and '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__))

    def __sub__(self, other):
        """Subtraction operator."""

        # Subtraction of another object of the same type.
        if type(other) is type(self):
            # Subtract the magnitudes in a common unit.
            mag = self._default_unit().magnitude() - other._default_unit().magnitude()

            # Return a new object of the same type with the original unit.
            return self._default_unit(mag)._convert_to(self._unit)

        # Addition of a string.
        elif type(other) is str:
            temp = self._from_string(other)
            return self - temp

        else:
            raise TypeError("unsupported operand type(s) for -: '%s' and '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__))

    def __mul__(self, other):
        """Multiplication operator."""

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Only support multiplication by float.
        if type(other) is float:
            # Convert to default unit and multiply.
            mag = self._default_unit().magnitude() * other

            # Return a new object of the same type with the original unit.
            return self._default_unit(mag)._convert_to(self._unit)

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
            # Convert to default unit and divide.
            mag = self._default_unit().magnitude() / other

            # Return a new object of the same type with the original unit.
            return self._default_unit(mag)._convert_to(self._unit)

        # Division by another object of the same type.
        elif type(other) is type(self):
            return self._default_unit().magnitude() / other._default_unit().magnitude()

        # Division by a string.
        elif type(other) is str:
            obj = self._from_string(other)
            return self / obj

        else:
            raise TypeError("unsupported operand type(s) for /: '%s' and '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__))

    def __lt__(self, other):
        """Less than operator."""

        # Compare to another object of the same type.
        if type(other) is type(self):
            return self._default_unit().magnitude() < other._default_unit().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self._default_unit().magnitude() < self._from_string(other)._default_unit().magnitude()

        else:
            raise TypeError("unorderable types: '%s' < '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__))

    def __le__(self, other):
        """Less than or equal to operator."""

        # Compare to another object of the same type.
        if type(other) is type(self):
            return self._default_unit().magnitude() <= other._default_unit().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self._default_unit().magnitude() <= self._from_string(other)._default_unit().magnitude()

        else:
            raise TypeError("unorderable types: '%s' <= '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__))

    def __eq__(self, other):
        """Equals to operator."""

        # Compare to another object of the same type.
        if type(other) is type(self):
            return self._default_unit().magnitude() == other._default_unit().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self._default_unit().magnitude() == self._from_string(other)._default_unit().magnitude()

        else:
            return False

    def __ne__(self, other):
        """Not equals to operator."""

        # Compare to another object of the same type.
        if type(other) is type(self):
            return self._default_unit().magnitude() != other._default_unit().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self._default_unit().magnitude() != self._from_string(other)._default_unit().magnitude()

        else:
            return True

    def __ge__(self, other):
        """Greater than or equal to operator."""

        # Compare to another object of the same type.
        if type(other) is type(self):
            return self._default_unit().magnitude() >= other._default_unit().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self._default_unit().magnitude() >= self._from_string(other)._default_unit().magnitude()

        else:
            raise TypeError("unorderable types: '%s' >= '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__))

    def __gt__(self, other):
        """Gretear than operator."""

        # Compare to another object of the same type.
        if type(other) is type(self):
            return self._default_unit().magnitude() > other._default_unit().magnitude()

        # Compare with a string.
        elif type(other) is str:
            return self._default_unit().magnitude() > self._from_string(other)._default_unit().magnitude()

        else:
            raise TypeError("unorderable types: '%s' > '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__))

    def magnitude(self):
        """Return the magnitude.

           Returns
           -------

           magnitude : float
               The magnitude of the type.
        """
        return self._magnitude

    def unit(self):
        """Return the unit.

           Returns
           -------

           unit : str
               The unit of the type.
        """
        return self._unit

    def _from_string(self, string):
        """Convert a string to an object of the same type.

           Parameters
           ----------

           string : str
               The string to interpret.

           Returns
           -------

           type : :class:`Type <BioSimSpace.Types>`
               The type object.
        """

        if string == "==SUPPRESS==":
            return type(self)(0, self._null_unit)

        string_copy = string

        if type(string) is str:
            # Strip white space from the string.
            string = string.replace(" ", "")

            # Try to match scientific format.
            match = _re.search(r"(\-?\d+\.?\d*e\-?\d+)(.*)", string, _re.IGNORECASE)

            # Try to match decimal format.
            if match is None:
                match = _re.search(r"(\-?\d+\.?\d*)(.*)", string, _re.IGNORECASE)

                # No matches, raise an error.
                if match is None:
                    raise ValueError("Could not interpret '%s' as '%s'"
                        % (string_copy, self.__class__.__qualname__))

            # Extract the value and unit.
            value, unit = match.groups()

            # Convert the value to a float.
            value = float(value)

            # Create and return a new object.
            return type(self)(value, unit)

        else:
            raise TypeError("'string' must be of type 'str'")
