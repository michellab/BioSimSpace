"""
@package biosimspace
@author  Lester Hedges
@brief   A collection of variable types.
"""

from math import floor

class Boolean():
    """A boolean type."""

    # The type of the value.
    value_type = bool

    def __init__(self, *value):
        """Constructor.

           Positional arguments:

           value -- The boolean value.
        """

        # Set the value.
        if len(value) == 0:
            self.value = False
        else:
            self.value = value[0]

    @property
    def value(self):
        """Get the boolean value."""
        return self._value

    @value.setter
    def value(self, value):
        """Set the boolean value."""

        # Python bool.
        if type(value) is bool:
            self._value = value

        # BioSimSpace bool.
        elif type(value) is Boolean:
            self._value = value.value

        # Python int.
        elif type(value) is int:
            if value == 0:
                self._value = False
            elif value == 1:
                self._value = True
            else:
                raise ValueError("Integer argument is not 0 or 1.")

        # BioSimSpace int.
        elif type(value) is Integer:
            if value.value == 0:
                self._value = False
            elif value.value == 1:
                self._value = True
            else:
                raise ValueError("Integer argument is not 0 or 1.")

        # Python string.
        elif type(value) is str:

            # Strip whitespace and convert to upper case.
            value = value.strip().upper()

            if value == 'TRUE':
                self._value = True
            elif value == 'FALSE':
                self._value = False
            else:
                raise ValueError("String argument is not 'True' or 'False'")

        # BioSimSpace string.
        elif type(value) is String:

            # Strip whitespace and convert to upper case.
            value = value.value.strip().upper()

            if value == 'TRUE':
                self._value = True
            elif value == 'FALSE':
                self._value = False
            else:
                raise ValueError("String argument is not 'True' or 'False'")

        else:
            raise ValueError("Cannot convert %s to %s" % (type(value), type(self)))

class Integer():
    """An integer type."""

    # The type of the value.
    value_type = int

    def __init__(self, *value):
        """Constructor.

           Positional arguments:

           value -- The value of the integer.
        """

        # Set the value.
        if len(value) == 0:
            self.value = 0
        else:
            self.value = value[0]

    @property
    def value(self):
        """Get the integer value."""
        return self._value

    @value.setter
    def value(self, value):
        """Set the value of the float."""

        # Python int.
        if type(value) is int:
            self._value = value

        # Python float.
        elif type(value) is float:
            self._value = floor(value)

        # BioSimSpace float.
        elif type(value) is Float:
            self._value = floor(value.value)

        # Python string.
        elif type(value) is str:
            self._value = floor(float(value))

        # BioSimSpace string.
        elif type(value) is String:
            self._value = floor(float(value.value))

        # Python bool.
        elif type(value) is bool:
            self._value = int(value)

        # BioSimSpace bool.
        elif type(value) is Boolean:
            self._value = int(value.value)

        else:
            raise ValueError("Cannot convert %s to %s" % (type(value), type(self)))

class Float():
    """An floating point type."""

    # The type of the value.
    value_type = float

    def __init__(self, *value):
        """Constructor.

           Positional arguments:

           value -- The value of the float.
        """

        # Set the value.
        if len(value) == 0:
            self.value = 0.0
        else:
            self.value = value[0]

    @property
    def value(self):
        """Get the float value."""
        return self._value

    @value.setter
    def value(self, value):
        """Set the value of the float."""

        # Python float.
        if type(value) is float:
            self._value = value

        # BioSimSpace float.
        elif type(value) is Float:
            self._value = value.value

        # Python int.
        elif type(value) is int:
            self._value = float(value)

        # BioSimSpace int.
        elif type(value) is Integer:
            self._value = float(value.value)

        # Python string.
        elif type(value) is str:
            self._value = float(value)

        # BioSimSpace string.
        elif type(value) is String:
            self._value = float(value.value)

        else:
            raise ValueError("Cannot convert %s to %s" % (type(value), type(self)))

class String():
    """A string type."""

    # The type of the value.
    value_type = str

    def __init__(self, *value):
        """Constructor.

           Positional arguments:

           value -- The string value.
        """

        # Set the value.
        if len(value) == 0:
            self.value = ''
        else:
            self.value = value[0]

    @property
    def value(self):
        """Get the string value."""
        return self._value

    @value.setter
    def value(self, value):
        """Set the string value."""

        # BioSimSpace type.
        if value.__class__.__module__ == 'BioSimSpace.Gateway.types':
            self._value = str(value.value)
        else:
            self._value = str(value)
