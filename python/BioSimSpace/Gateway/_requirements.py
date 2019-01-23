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
Functionality for defining and validating BioSimSpace input and output requirements.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

import BioSimSpace.Types as _Types

import bz2 as _bz2
import copy as _copy
import gzip as _gzip
import os as _os
import re as _re
import shutil as _shutil
import sys as _sys
import tarfile as _tarfile
import zipfile as _zipfile

__all__ = ["Boolean", "Integer", "Float", "String",     # Regular types.
           "File", "FileSet",                           # File types.
           "Length", "Area", "Volume",                  # Length types.
           "Charge",
           "Energy",
           "Pressure",
           "Temperature",
           "Time"]

class Requirement():
    """Base class for BioSimSpace Node requirements."""

    # Set the argparse argument type.
    _arg_type = None

    # Default to single arguments.
    _is_multi = False

    def __init__(self, help=None, default=None, unit=None, minimum=None,
            maximum=None, allowed=None, optional=False):
        """Constructor.


           Parameters
           ----------

           help : str
               The help string.

           default :
               The default value.

           unit : str
               The unit.

           minimum :
               The minimum allowed value.

           maximum :
               The maximum allowed value.

           allowed : list
               A list of allowed values.

           optional : bool
               Whether the requirement is optional.
        """

	# Don't allow user to create an instance of this base class.
        if type(self) is Requirement:
            raise Exception("<Requirement> must be subclassed.")

        # Required keyword arguments.

        if help is None:
            raise ValueError("Missing 'help' keyword argument!")
        elif type(help) is not str:
            raise TypeError("'help' keyword argument must be of type 'str'")
        else:
            self._help = help

        # Optional keywords aguments.

        if type(optional) is not bool:
            raise TypeError("'optional' keyword argument must be of type 'bool'")

        # Set defaults
        self._value = None
        self._default = None
        self._unit = unit
        self._min = None
        self._max = None
        self._allowed = None
        self._is_optional = optional

        # Set the minimum value.
        if minimum is not None:
            self._min = self._validate(minimum)

        # Set the maximum value.
        if maximum is not None:
            self._max = self._validate(maximum)

        # Min and max.
        if minimum is not None and maximum is not None:
            if self._max < self._min:
                raise ValueError("The maximum value '%s' is less than the minimum '%s'"
                    % (maximum, minimum))

        # Set the allowed values.
        if allowed is not None:
            if type(allowed) is not list:
                allowed = [allowed]
            self._allowed = [self._validate(x) for x in allowed]

            # Conflicting requirements.
            if self._min is not None or self._max is not None:
                raise ValueError("Conflicting requirement: cannot have allowed values and min/max.")

        # Set the default value.
        if default is not None:
            # Make sure the default is the right type.
            self._default = self._validate(default)
            # Make sure the default satisfies the constraints.
            self._validate_default()
            # Flag the the requirement is optional.
            self._is_optional = True

    def setValue(self, value, name):
        """Validate and set the value.


           Parameters
           ----------

           value :
               The value of the input requirement.

           name : str
               The name of the requirement.
        """

        if value is None and not self._is_optional:
            if type(name) is str:
                raise ValueError("Value is unset for requirement '%s'!" % name)
            else:
                raise ValueError("Value is unset!")

        # Validate the value.
        value = self._validate(value)

        # Now check value against any constraints.

        # Minimum.
        if self._min is not None and value < self._min:
            raise ValueError("The value (%s) is less than the allowed "
                "minimum (%s)" % (value, self._min))

        # Maximum.
        if self._max is not None and value > self._max:
            raise ValueError("The value (%s) is less than the allowed "
                "maximum (%s)" % (value, self._max))

        # Allowed values.
        if self._allowed is not None and value not in self._allowed:
            raise ValueError("The value (%s) is not in the list of allowed values: "
                "%s" % (value, str(self._allowed)))

        # All is okay. Set the value.
        self._value = value

    def getValue(self):
        """Return the value."""
        return self._value

    def getDefault(self):
        """Return the default value."""
        return self._default

    def getUnit(self):
        """Return the unit."""
        return self._unit

    def getHelp(self):
        """Return the documentation string."""
        return self._help

    def isMulti(self):
        """Whether the requirement has multiple values."""
        return self._is_multi

    def isOptional(self):
        """Whether the requirement is optional."""
        return self._is_optional

    def getArgType(self):
        """The command-line argument type."""
        return self._arg_type

    def getMin(self):
        """Return the minimum allowed value."""
        return self._min

    def getMax(self):
        """Return the maximum allowed value."""
        return self._max

    def getAllowedValues(self):
        """Return the allowed values."""
        return self._allowed

    def _validate_default(self):
        """Validate the default value."""

        if self._min is not None and self._default < self._min:
            raise ValueError("The default '%s' is less than the minimum "
                "allowed value '%s'" % (self._default, self._min))

        if self._max is not None and self._default > self._max:
            raise ValueError("The default '%s' is greater than the maximum "
                "allowed value '%s'" % (self._default, self._max))

        if self._allowed is not None and self._default not in self._allowed:
            raise ValueError("The default '%s' is not one of the allowed values %s"
                % (self._default, str(self._allowed)))

class Boolean(Requirement):
    """A boolean requirement."""

    # Set the argparse argument type.
    _arg_type = bool

    def __init__(self, help=None, default=None):
        """Constructor.


           Parameters
           ----------

           help : str
               The help string.

           default :
               The default value.
        """

        # Call the base class constructor.
        super().__init__(help=help, default=default)

    def _validate(self, value):
        """Validate the value."""

        if type(value) is bool:
            return value
        else:
            raise TypeError("The value should be of type 'bool'.")

class Integer(Requirement):
    """An integer requirement."""

    # Set the argparse argument type.
    _arg_type = int

    def __init__(self, help=None, default=None,
            minimum=None, maximum=None, allowed=None):
        """Constructor.


           Parameters
           ----------

           help : str
               The help string.

           default : int
               The default value.

           minimum : int
               The minimum allowed value.

           maximum : int
               The maximum allowed value.

           allowed : [ int ]
               A list of allowed values.
        """

        # Call the base class constructor.
        super().__init__(help=help, default=default, minimum=minimum,
            maximum=maximum, allowed=allowed)

    def _validate(self, value):
        """Validate that the value is of the correct type."""

        if type(value) is int:
            return value
        else:
            raise TypeError("The value should be of type 'int'.")

class Float(Requirement):
    """A floating point requirement."""

    # Set the argparse argument type.
    _arg_type = float

    def __init__(self, help=None, default=None,
            minimum=None, maximum=None, allowed=None):
        """Constructor.


           Parameters
           ----------

           help : str
               The help string.

           default : float
               The default value.

           minimum : float
               The minimum allowed value.

           maximum : float
               The maximum allowed value.

           allowed : [ float ]
               A list of allowed values.
        """

        # Call the base class constructor.
        super().__init__(help=help, default=default, minimum=minimum,
            maximum=maximum, allowed=allowed)

    def _validate(self, value):
        """Validate that the value is of the correct type."""

        if type(value) is float:
            return value
        elif type(value) is int:
            return float(value)
        else:
            raise TypeError("The value should be of type 'float' or 'int'.")

class String(Requirement):
    """A string requirement."""

    # Set the argparse argument type.
    _arg_type = str

    def __init__(self, help=None, default=None, allowed=None):
        """Constructor.


           Parameters
           ----------

           help : str
               The help string.

           default : str
               The default value.

           allowed : [ str ]
               A list of allowed values.
        """

        # Call the base class constructor.
        super().__init__(help=help, default=default, allowed=allowed)

    def _validate(self, value):
        """Validate that the value is of the correct type."""

        if type(value) is str:
            return value
        else:
            raise TypeError("The value should be of type 'str'")

class File(Requirement):
    """A file set requirement."""

    # Set the argparse argument type.
    _arg_type = str

    def __init__(self, help=None, optional=False):
        """Constructor.


           Parameters
           ----------

           help : str
               The help string.

           optional : bool
               Whether the file is optional.
        """

        # Call the base class constructor.
        super().__init__(help=help, optional=optional)

    def _validate(self, value):
        """Validate that the value is of the correct type."""

        # Handle optional requirement.
        if self._is_optional and value is None:
            return None

        # Check the type.
        if type(value) is str:
            file = _unarchive(value)
            if file is None:
                file = value
        else:
            raise TypeError("The value should be of type 'str'")

        # Make sure the file exists.
        if not _os.path.isfile(file):
            raise IOError("File doesn't exist: '%s'" % file)
        else:
            return file

class FileSet(Requirement):
    """A file requirement."""

    # Set the argparse argument type.
    _arg_type = str

    # Multiple files can be passed.
    _is_multi = True

    def __init__(self, help=None, optional=False):
        """Constructor.


           Parameters
           ----------

           help : str
               The help string.

           optional : bool
               Whether the file set is optional.
        """

        # Call the base class constructor.
        super().__init__(help=help, optional=optional)

    def getValue(self):
        """Return the value."""
        if self._value is None:
            return None
        else:
            return self._value.copy()

    def _validate(self, value):
        """Validate that the value is of the correct type."""

        # Handle optional requirement.
        if self._is_optional and value is None:
            return None

        # Handle single strings.
        if type(value) is str:
            value = [value]

        # The user can pass a list of compressed files so we need to keep
        # track of the names of the uncompressed files.
        uncompressed_files = []

        # We should receive a list of strings.
        if type(value) is list:

            # A single file was passed.
            if len(value) == 1:
                # Remove whitespace and split on commas.
                value = _re.sub(r"\s+", "", value[0]).split(',')

            # Loop over all strings.
            for file in value:

                # Check the types.
                if type(file) is not str:
                    raise TypeError("The value should be of type 'str'")

                # Check whether this is an archive.
                files = _unarchive(file)

                if files is not None:
                    if type(files) is list:
                        uncompressed_files += files
                    else:
                        uncompressed_files.append(files)
                else:
                    # Make sure the file exists.
                    if not _os.path.isfile(file):
                        raise IOError("File doesn't exist: '%s'" % file)

        if len(uncompressed_files) > 0:
            return uncompressed_files
        else:
            return value

class Length(Requirement):
    """A length requirement."""

    # Set the argparse argument type.
    _arg_type = str

    def __init__(self, help=None, default=None, unit=None,
            minimum=None, maximum=None, allowed=None):
        """Constructor.


           Parameters
           ----------

           help : str
               The help string.

           default : BioSimSpace.Types.Length
               The default value.

           unit : str
               The unit.

           minimum : BioSimSpace.Types.Length
               The minimum allowed value.

           maximum : BioSimSpace.Types.Length
               The maximum allowed value.

           allowed : [ BioSimSpace.Types.Length ]
        """

        # Validate the unit.
        if unit is not None:
            length = _Types.Length("1 %s" % unit)
            self._unit = length.unit()
            self._print_unit = length._print_format[length.unit()]
        else:
            raise ValueError("No unit has been specified!")

        # Call the base class constructor.
        super().__init__(help=help, default=default, unit=self._unit, minimum=minimum,
            maximum=maximum, allowed=allowed)

    def getValue(self):
        """Return the value."""
        if self._value is None:
            return None
        else:
            return _copy.deepcopy(self._value)

    def _validate(self, value):
        """Validate that the value is of the correct type."""

        if type(value) is _Types.Length:
            return value._convert_to(self._unit)

        else:
            # Extract the value and unit from the argument string.
            value, unit = _validate_unit_requirement(value, "length")

            if unit is None:
                return _Types.Length(value, self._unit)
            else:
                return _Types.Length(value, unit)._convert_to(self._unit)

class Area(Requirement):
    """An area requirement."""

    # Set the argparse argument type.
    _arg_type = str

    def __init__(self, help=None, default=None, unit=None,
            minimum=None, maximum=None, allowed=None):
        """Constructor.


           Parameters
           ----------

           help : str
               The help string.

           default : BioSimSpace.Types.Area
               The default value.

           unit : str
               The unit.

           minimum : BioSimSpace.Types.Area
               The minimum allowed value.

           maximum : BioSimSpace.Types.Area
               The maximum allowed value.

           allowed : [ BioSimSpace.Types.Area ]
        """

        # Validate the unit.
        if unit is not None:
            area = _Types.Area("1 %s" % unit)
            self._unit = area.unit()
            self._print_unit = area._print_format[area.unit()]
        else:
            raise ValueError("No unit has been specified!")

        # Call the base class constructor.
        super().__init__(help=help, default=default, unit=self._unit,
            minimum=minimum, maximum=maximum, allowed=allowed)

    def getValue(self):
        """Return the value."""
        if self._value is None:
            return None
        else:
            return _copy.deepcopy(self._value)

    def _validate(self, value):
        """Validate that the value is of the correct type."""

        if type(value) is _Types.Area:
            return value._convert_to(self._unit)

        else:
            # Extract the value and unit from the argument string.
            value, unit = _validate_unit_requirement(value, "area")

            if unit is None:
                return _Types.Area(value, self._unit)
            else:
                return _Types.Area(value, unit)._convert_to(self._unit)

class Volume(Requirement):
    """A volume requirement."""

    # Set the argparse argument type.
    _arg_type = str

    def __init__(self, help=None, default=None, unit=None,
            minimum=None, maximum=None, allowed=None):
        """Constructor.


           Parameters
           ----------

           help : str
               The help string.

           default : BioSimSpace.Types.Volume
               The default value.

           unit : str
               The unit.

           minimum : BioSimSpace.Types.Volume
               The minimum allowed value.

           maximum : BioSimSpace.Types.Volume
               The maximum allowed value.

           allowed : [ BioSimSpace.Types.Volume ]
        """

        # Validate the unit.
        if unit is not None:
            volume = _Types.Volume("1 %s" % unit)
            self._unit = volume.unit()
            self._print_unit = volume._print_format[volume.unit()]
        else:
            raise ValueError("No unit has been specified!")

        # Call the base class constructor.
        super().__init__(help=help, default=default, unit=self._unit,
            minimum=minimum, maximum=maximum, allowed=allowed)

    def getValue(self):
        """Return the value."""
        if self._value is None:
            return None
        else:
            return _copy.deepcopy(self._value)

    def _validate(self, value):
        """Validate that the value is of the correct type."""

        if type(value) is _Types.Volume:
            return value._convert_to(self._unit)

        else:
            # Extract the value and unit from the argument string.
            value, unit = _validate_unit_requirement(value, "volume")

            if unit is None:
                return _Types.Volume(value, self._unit)
            else:
                return _Types.Volume(value, unit)._convert_to(self._unit)

class Charge(Requirement):
    """A charge requirement."""

    # Set the argparse argument type.
    _arg_type = str

    def __init__(self, help=None, default=None, unit=None,
            minimum=None, maximum=None, allowed=None):
        """Constructor.


           Parameters
           ----------

           help : str
               The help string.

           default : BioSimSpace.Types.Charge
               The default value.

           unit : str
               The unit.

           minimum : BioSimSpace.Types.Charge
               The minimum allowed value.

           maximum : BioSimSpace.Types.Charge
               The maximum allowed value.

           allowed : [ BioSimSpace.Types.Charge ]
        """

        # Validate the unit.
        if unit is not None:
            nrg = _Types.Charge("1 %s" % unit)
            self._unit = nrg.unit()
            self._print_unit = nrg._print_format[nrg.unit()]
        else:
            raise ValueError("No unit has been specified!")

        # Call the base class constructor.
        super().__init__(help=help, default=default, unit=self._unit,
            minimum=minimum, maximum=maximum, allowed=allowed)

    def getValue(self):
        """Return the value."""
        if self._value is None:
            return None
        else:
            return _copy.deepcopy(self._value)

    def _validate(self, value):
        """Validate that the value is of the correct type."""

        if type(value) is _Types.Charge:
            return value._convert_to(self._unit)

        else:
            # Extract the value and unit from the argument string.
            value, unit = _validate_unit_requirement(value, "energy")

            if unit is None:
                return _Types.Charge(value, self._unit)
            else:
                return _Types.Charge(value, unit)._convert_to(self._unit)

class Energy(Requirement):
    """An energy requirement."""

    # Set the argparse argument type.
    _arg_type = str

    def __init__(self, help=None, default=None, unit=None,
            minimum=None, maximum=None, allowed=None):
        """Constructor.


           Parameters
           ----------

           help : str
               The help string.

           default : BioSimSpace.Types.Energy
               The default value.

           unit : str
               The unit.

           minimum : BioSimSpace.Types.Energy
               The minimum allowed value.

           maximum : BioSimSpace.Types.Energy
               The maximum allowed value.

           allowed : [ BioSimSpace.Types.Energy ]
        """

        # Validate the unit.
        if unit is not None:
            nrg = _Types.Energy("1 %s" % unit)
            self._unit = nrg.unit()
            self._print_unit = nrg._print_format[nrg.unit()]
        else:
            raise ValueError("No unit has been specified!")

        # Call the base class constructor.
        super().__init__(help=help, default=default, unit=self._unit,
            minimum=minimum, maximum=maximum, allowed=allowed)

    def getValue(self):
        """Return the value."""
        if self._value is None:
            return None
        else:
            return _copy.deepcopy(self._value)

    def _validate(self, value):
        """Validate that the value is of the correct type."""

        if type(value) is _Types.Energy:
            return value._convert_to(self._unit)

        else:
            # Extract the value and unit from the argument string.
            value, unit = _validate_unit_requirement(value, "energy")

            if unit is None:
                return _Types.Energy(value, self._unit)
            else:
                return _Types.Energy(value, unit)._convert_to(self._unit)

class Pressure(Requirement):
    """A pressure requirement."""

    # Set the argparse argument type.
    _arg_type = str

    def __init__(self, help=None, default=None, unit=None,
            minimum=None, maximum=None, allowed=None):
        """Constructor.


           Parameters
           ----------

           help : str
               The help string.

           default : BioSimSpace.Types.Pressure
               The default value.

           unit : str
               The unit.

           minimum : BioSimSpace.Types.Pressure
               The minimum allowed value.

           maximum : BioSimSpace.Types.Pressure
               The maximum allowed value.

           allowed : [ BioSimSpace.Types.Pressure ]
        """

        # Validate the unit.
        if unit is not None:
            press = _Types.Pressure("1 %s" % unit)
            self._unit = press.unit()
            self._print_unit = press._print_format[press.unit()]
        else:
            raise ValueError("No unit has been specified!")

        # Call the base class constructor.
        super().__init__(help=help, default=default, unit=self._unit,
            minimum=minimum, maximum=maximum, allowed=allowed)

    def getValue(self):
        """Return the value."""
        if self._value is None:
            return None
        else:
            return _copy.deepcopy(self._value)

    def _validate(self, value):
        """Validate that the value is of the correct type."""

        if type(value) is _Types.Pressure:
            return value._convert_to(self._unit)

        else:
            # Extract the value and unit from the argument string.
            value, unit = _validate_unit_requirement(value, "pressure")

            if unit is None:
                return _Types.Pressure(value, self._unit)
            else:
                return _Types.Pressure(value, unit)._convert_to(self._unit)

class Temperature(Requirement):
    """A temperature requirement."""

    # Set the argparse argument type.
    _arg_type = str

    def __init__(self, help=None, default=None, unit=None,
            minimum=None, maximum=None, allowed=None):
        """Constructor.


           Parameters
           ----------

           help : str
               The help string.

           default : BioSimSpace.Types.Temperature
               The default value.

           unit : str
               The unit.

           minimum : BioSimSpace.Types.Temperature
               The minimum allowed value.

           maximum : BioSimSpace.Types.Temperature
               The maximum allowed value.

           allowed : [ BioSimSpace.Types.Temperature ]
        """

        # Validate the unit.
        if unit is not None:
            temp = _Types.Temperature("1 %s" % unit)
            self._unit = temp.unit()
            self._print_unit = temp._print_format[temp.unit()]
        else:
            raise ValueError("No unit has been specified!")

        # Call the base class constructor.
        super().__init__(help=help, default=default, unit=self._unit,
            minimum=minimum, maximum=maximum, allowed=allowed)

    def getValue(self):
        """Return the value."""
        if self._value is None:
            return None
        else:
            return _copy.deepcopy(self._value)

    def _validate(self, value):
        """Validate that the value is of the correct type."""

        if type(value) is _Types.Temperature:
            return value._convert_to(self._unit)

        else:
            # Extract the value and unit from the argument string.
            value, unit = _validate_unit_requirement(value, "temperature")

            if unit is None:
                return _Types.Temperature(value, self._unit)
            else:
                return _Types.Temperature(value, unit)._convert_to(self._unit)

class Time(Requirement):
    """A time requirement."""

    # Set the argparse argument type.
    _arg_type = str

    def __init__(self, help=None, default=None, unit=None,
            minimum=None, maximum=None, allowed=None):
        """Constructor.


           Parameters
           ----------

           help : str
               The help string.

           default : BioSimSpace.Types.Time
               The default value.

           unit : str
               The unit.

           minimum : BioSimSpace.Types.Time
               The minimum allowed value.

           maximum : BioSimSpace.Types.Time
               The maximum allowed value.

           allowed : [ BioSimSpace.Types.Time ]
        """

        # Validate the unit.
        if unit is not None:
            time = _Types.Time("1 %s" % unit)
            self._unit = time.unit()
            self._print_unit = time._print_format[time.unit()]
        else:
            raise ValueError("No unit has been specified!")

        # Call the base class constructor.
        super().__init__(help=help, default=default, unit=self._unit,
            minimum=minimum, maximum=maximum, allowed=allowed)

    def getValue(self):
        """Return the value."""
        if self._value is None:
            return None
        else:
            return _copy.deepcopy(self._value)

    def _validate(self, value):
        """Validate that the value is of the correct type."""

        if type(value) is _Types.Time:
            return value._convert_to(self._unit)

        else:
            # Extract the value and unit from the argument string.
            value, unit = _validate_unit_requirement(value, "time")

            if unit is None:
                return _Types.Time(value, self._unit)
            else:
                return _Types.Time(value, unit)._convert_to(self._unit)

def _validate_unit_requirement(value, unit_type):
    """Helper function to validate input requirements with units.


       Parameters
       ----------

       value : str
           The value of the input requirement.

       unit_type: str
           The unit type.


        Returns
        -------

        (value, unit) : tuple
            The value and unit of the requirement.
    """

    # No unit by default.
    unit = None

    # Float.
    if type(value) is float:
        pass

    # Integer.
    elif type(value) is int:
        value = float(value)

    # String.
    elif type(value) is str:
        # First try to directly convert to a float.
        try:
            value = float(value)

        # Use a regular expression to extract the value and unit.
        except ValueError:

            # Strip white space from the string.
            value = value.replace(" ", "")

            # Try to match scientific format.
            match = _re.search("(\-?\d+\.?\d*e\-?\d+)(.*)", value, _re.IGNORECASE)

            # Try to match decimal format.
            if match is None:
                match = _re.search("(\-?\d+\.?\d*)(.*)", value, _re.IGNORECASE)

                # No matches, raise an error.
                if match is None:
                    raise ValueError("Could not interpret %s: '%s'" % (unit_type, value))

            # Extract the value and unit.
            value, unit = match.groups()

            # Convert the value to a float.
            value = float(value)

    # Unsupported.
    else:
        raise TypeError("The value should be of type 'float', 'int', or 'str'")

    return (value, unit)

def _unarchive(name):
    """Decompress an archive and return a list of files.


       Parameters
       ----------

       name : str
           The name of the archive (full path).


       Returns
       -------

       files : [ str ]
           A list of file names.
    """

    # Get the directory name.
    dir = _os.path.dirname(name)

    # If the file compressed file has been passed on the command-line, then
    # we'll extract it to a directory to avoid littering the current workspace.
    if dir == "uploads":
        dir += "/"
    else:
        dir = "uncompressed/"

    # List of supported tar file formats.
    tarfiles = ["tar.gz", "tar.bz2", "tar"]

    # Check whether this is a tar compressed file.
    for tar_name in tarfiles:

        # Found a match.
        if tar_name in name.lower():

            # The list of decompressed files.
            files = []

            # Decompress the archive.
            with _tarfile.open(name) as tar:
                # We need to call tar.list(), otherwise the tar object will not know
                # about nested directories, i.e. it will appear as if ther is a single
                # member.
                print("Decompressing...")
                tar.list(verbose=False)

                # Loop over all of the members and get the file names.
                # If the name has no extension, then we assume that it's a directory.
                for file in tar.members:
                    if _os.path.splitext(file.name)[1] is not "":
                        files.append(dir + file.name)

                # Now extract all of the files.
                tar.extractall(dir)

            return files

    # Get the file name and extension.
    file, ext = _os.path.splitext(name)

    # This is a zip file.
    if ext.lower() == ".zip":
        # The list of decompressed files.
        files = None

        with _zipfile.ZipFile(name) as zip:
            files = zip.namelist()
            print("Decompressing...")
            for file in files:
                print(file)
            zip.extractall(dir)

        return files

    # This is a gzip file.
    if ext.lower() == ".gz" or ext == ".gzip":
        with _gzip.open(name, "rb") as f_in:
            with open(file, "wb") as f_out:
                print("Decompressing...\n%s" % name)
                _shutil.copyfileobj(f_in, f_out)

        return file

    # This is a bzip2 file.
    if ext.lower() == ".bz2" or ext == ".bzip2":
        with _bz2.open(name, "rb") as f_in:
            with open(file, "wb") as f_out:
                print("Decompressing...\n%s" % name)
                _shutil.copyfileobj(f_in, f_out)

        return file

    # If we get this far, then this is not a supported archive.
    return None
