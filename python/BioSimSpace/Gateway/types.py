"""
@package biosimspace
@author  Lester Hedges
@brief   A collection of requirement classes.
"""

from os import path

class Requirement():
    """Base class for BioSimSpace Node requirements."""

    # Set the argparse argument type.
    _arg_type = None

    # Default to single arguments.
    _is_multi = False

    def __init__(self, name=None, help=None, default=None,
            units=None, minimum=None, maximum=None, allowed=None):
        """Constructor.

           Keyword arguments:

           name    -- The name of the requirement.
           help    -- The help string.
           default -- The default value.
           units   -- The units.
           minimum -- The minimum allowed value.
           maximum -- The maximum allowed value.
           allowed -- A list of allowed values.
        """

	# Don't allow user to create an instance of this base class.
        if type(self) == Requirement:
            raise Exception("<Requirement> must be subclassed.")

        # Required keyword arguments.

        if name is None:
            raise ValueError("Missing 'name' keyword argument!")
        elif type(name) is not str:
            raise ValueError("'name' keyword argument must be of type 'str'")
        else:
            self._name = name

        if help is None:
            raise ValueError("Missing 'help' keyword argument!")
        elif type(help) is not str:
            raise ValueError("'help' keyword argument must be of type 'str'")
        else:
            self._help = help

        # Set defaults.
        self._value = None

        # Set member data.
        self._default = default
        self._units = None
        self._min = minimum
        self._max = maximum
        self._allowed = allowed

    def setValue(self, value):
        """Validate and set the value."""

        # Validate the value.
        value = self._validate(value)

        # Now check value against any constraints.

        # Minimum.
        if self._min is not None and value < self._min:
            raise ValueError("The value (%d) is less than the allowed "
                "minimum (%d)" % (value, self._min))

        # Maximum.
        if self._max is not None and value < self._max:
            raise ValueError("The value (%d) is less than the allowed "
                "maximum (%d)" % (value, self._max))

        # Allowed values.
        if self._allowed is not None and value not in self._allowed:
            raise ValueError("The value (%d) is not in the list of allowed values: "
                "%s" % (value, str(self._allowed)))

        # All is okay. Set the value.
        self._value = value

    def value(self):
        """Return the value."""
        return self._value

    def name(self):
        """Return the name of the requirement."""
        return self._name

    def default(self):
        """Return the default value."""
        return self._default

    def units(self):
        """Return the units."""
        return self._units

    def helpText(self):
        """Return the documentation string."""
        return self._help

    def isMulti(self):
        """Whether the requirement has multiple values."""
        return self._is_multi

    def argType(self):
        """The command-line argument type."""
        return self._arg_type

    def min(self):
        """Return the minimum allowed value."""
        return self._min

    def max(self):
        """Return the maximum allowed value."""
        return self._max

    def allowed(self):
        """Return the allowed values."""
        return self._allowed

class Boolean(Requirement):
    """A boolean requirement."""

    # Set the argparse argument type.
    _arg_type = bool

    def __init__(self, name=None, help=None, default=None):
        """Constructor.

           Keyword arguments:

           name    -- The name of the requirement.
           help    -- The help string.
           default -- The default value.
        """

        # Call the base class constructor.
        super().__init__(name=name, help=help, default=default)

    def _validate(self, value):
        """Validate the value."""

        # Python bool.
        if type(value) is bool:
            return value

        # BioSimSpace bool.
        elif type(value) is Boolean:
            return value.value()

        else:
            raise ValueError("Cannot convert '%s' to '%s'" % (type(value), type(self)))

class Integer(Requirement):
    """An integer requirement."""

    # Set the argparse argument type.
    _arg_type = int

    def __init__(self, name=None, help=None, default=None,
            minimum=None, maximum=None, allowed=None):
        """Constructor.

           Keyword arguments:

           name    -- The name of the requirement.
           help    -- The help string.
           default -- The default value.
           min     -- The minimum allowed value.
           max     -- The maximum allowed value.
           allowed -- A list of allowed values.
        """

        # Call the base class constructor.
        super().__init__(name=name, help=help)

        # Set the minimum value.
        if minimum is not None:
            self._min = self._validate(minimum)

        # Set the maximum value.
        if maximum is not None:
            self._max = self._validate(maximum)

        # Set the allowed values.
        if allowed is not None:
            self._allowed = [self._validate(x) for x in allowed]

        # Set the default value.
        if default is not None:
            self._default = self._validate(default)

    def _validate(self, value):
        """Validate that the value is of the correct type."""

        # Python int.
        if type(value) is int:
            return value

        # BioSimSpace Integer.
        elif type(value) is Integer:
            return value.value()

        else:
            raise ValueError("Cannot convert '%s' to '%s'" % (type(value), type(self)))

class Float(Requirement):
    """A floating point requirement."""

    # Set the argparse argument type.
    _arg_type = float

    def __init__(self, name=None, help=None, default=None,
            minimum=None, maximum=None, allowed=None):
        """Constructor.

           Keyword arguments:

           name    -- The name of the requirement.
           help    -- The help string.
           default -- The default value.
           min     -- The minimum allowed value.
           max     -- The maximum allowed value.
           allowed -- A list of allowed values.
        """

        # Call the base class constructor.
        super().__init__(name=name, help=help)

        # Set the minimum value.
        if minimum is not None:
            self._min = self._validate(minimum)

        # Set the maximum value.
        if maximum is not None:
            self._max = self._validate(maximum)

        # Set the allowed values.
        if allowed is not None:
            self._allowed = [self._validate(x) for x in allowed]

        # Set the default value.
        if default is not None:
            self._default = self._validate(default)

    def _validate(self, value):
        """Validate that the value is of the correct type."""

        # Python float.
        if type(value) is float:
            return value

        # BioSimSpace Float.
        elif type(value) is Float:
            return value.value()

        else:
            raise ValueError("Cannot convert '%s' to '%s'" % (type(value), type(self)))

class String(Requirement):
    """A string requirement."""

    # Set the argparse argument type.
    _arg_type = str

    def __init__(self, name=None, help=None, default=None, allowed=None):
        """Constructor.

           Keyword arguments:

           name    -- The name of the requirement.
           help    -- The help string.
           default -- The default value.
           allowed -- A list of allowed values.
        """

        # Call the base class constructor.
        super().__init__(name=name, help=help)

        # Set the allowed values.
        if allowed is not None:
            self._allowed = [self._validate(x) for x in allowed]

        # Set the default value.
        if default is not None:
            self._default = self._validate(default)

    def _validate(self, value):
        """Validate that the value is of the correct type."""

        # Python str.
        if type(value) is str:
            return value

        # BioSimSpace String.
        elif type(value) is String:
            return value.value()

        else:
            raise ValueError("Cannot convert '%s' to '%s'" % (type(value), type(self)))

class File(Requirement):
    """A file set requirement."""

    # Set the argparse argument type.
    _arg_type = str

    def __init__(self, name=None, help=None):
        """Constructor.

           Keyword arguments:

           name    -- The name of the requirement.
           help    -- The help string.
        """

        # Call the base class constructor.
        super().__init__(name=name, help=help)

    def _validate(self, value):
        """Validate that the value is of the correct type."""

        # Check the type.
        if type(value) is str:
            file = value
        elif type(value) is String:
            value = value.value()
        else:
            raise ValueError("Cannot convert '%s' to '%s'" % (type(file), type(self)))

        # Make sure the file exists.
        if not path.isfile(file):
            raise IOError(('File doesn\'t exist: "{x}"').format(x=file))
        else:
            return file

class FileSet(Requirement):
    """A file requirement."""

    # Set the argparse argument type.
    _arg_type = str

    # Multiple files can be passed.
    _is_multi = True

    def __init__(self, name=None, help=None):
        """Constructor.

           Keyword arguments:

           name    -- The name of the requirement.
           help    -- The help string.
        """

        # Call the base class constructor.
        super().__init__(name=name, help=help)

    def _validate(self, value):
        """Validate that the value is of the correct type."""

        # We should receive a list of strings.
        if type(value) is list:

            # Loop over all strings.
            for file in value:

                # Check the types.
                if type(file) is str:
                    file = file
                elif type(file) is String:
                    file = file.value()
                else:
                    raise ValueError("Cannot convert '%s' to '%s'" % (type(file), type(self)))

            # Make sure the file exists.
            if not path.isfile(file):
                raise IOError(('File doesn\'t exist: "{x}"').format(x=file))

        # All is okay. Return the value.
        return value
