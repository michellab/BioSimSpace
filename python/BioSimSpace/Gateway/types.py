"""
@package biosimspace
@author  Lester Hedges
@brief   A collection of variable types.
"""

class Boolean():
    """A boolean type."""

    def __init__(self, value):
        """Constructor.

           Positional arguments:

           value -- The boolean value: True or False.
        """

        # Set the value.
        self.value = value

    @property
    def value(self):
        """Get the boolean value."""
        return self._value

    @value.setter
    def value(self, value):
        """Set the boolean value."""
        if type(value) is bool:
            self._value = value
        else:
            raise TypeError("Argument is not of type 'bool'")

class Number():
    """A number type. Handles integers and floats."""

    def __init__(self, value):
        """Constructor.

           Positional arguments:

           value -- The value of the number.
        """

        # Set the value.
        self.value = value

    @property
    def value(self):
        """Get the boolean value."""
        return self._value

    @value.setter
    def value(self, value):
        """Set the value of the number."""

        # Integer.
        if type(value) is int:
            self._value = value
            self._is_float = False
            self._is_whole_number = True

        # Float.
        elif type(value) is float:
            self._value = value
            self._is_float = True

            # If the remainder is very small, the this can be treated
            # as a whole number.
            if value % 1 < 1e-6:
                self._is_whole_number = True

        else:
            raise TypeError("Argument is not of type 'int' or 'float'")

    def isFloat(self):
        """Whether the number is a float."""
        return self._is_float

    def isWholeNumber(self):
        """Whether this is a whole number."""
        return self._is_whole_number:
