######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2018
#
# Authors: Lester Hedges

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
Author: Lester Hedges
"""

import bz2
import gzip
import os
import re
import shutil
import sys
import tarfile
import zipfile

__all__ = ["Boolean", "Integer", "Float", "String", "File", "FileSet"]

class Requirement():
    """Base class for BioSimSpace Node requirements."""

    # Set the argparse argument type.
    _arg_type = None

    # Default to single arguments.
    _is_multi = False

    def __init__(self, help=None, default=None, units=None, minimum=None,
            maximum=None, allowed=None, optional=False):
        """Constructor.

           Keyword arguments:

           help     -- The help string.
           default  -- The default value.
           units    -- The units.
           minimum  -- The minimum allowed value.
           maximum  -- The maximum allowed value.
           allowed  -- A list of allowed values.
           optional -- Whether the requirement is optional.
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

        # Set defaults.
        self._value = None

        # Set member data.
        self._default = default
        self._units = None
        self._min = minimum
        self._max = maximum
        self._allowed = allowed
        self._is_optional = optional

    def setValue(self, value):
        """Validate and set the value."""

        if value is None and not self._is_optional:
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

    def getUnits(self):
        """Return the units."""
        return self._units

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

           Keyword arguments:

           help    -- The help string.
           default -- The default value.
        """

        # Call the base class constructor.
        super().__init__(help=help)

        # Set the default value.
        if default is not None:
            # Make sure the default is the right type.
            self._default = self._validate(default)
            # Flag the the requirement is optional.
            self._is_optional = True

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

           Keyword arguments:

           help    -- The help string.
           default -- The default value.
           min     -- The minimum allowed value.
           max     -- The maximum allowed value.
           allowed -- A list of allowed values.
        """

        # Call the base class constructor.
        super().__init__(help=help)

        # Set the minimum value.
        if minimum is not None:
            self._min = self._validate(minimum)

        # Set the maximum value.
        if maximum is not None:
            self._max = self._validate(maximum)

        # Min and max.
        if minimum is not None and maximum is not None:
            if maximum < minimum:
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

           Keyword arguments:

           help    -- The help string.
           default -- The default value.
           min     -- The minimum allowed value.
           max     -- The maximum allowed value.
           allowed -- A list of allowed values.
        """

        # Call the base class constructor.
        super().__init__(help=help)

        # Set the minimum value.
        if minimum is not None:
            self._min = self._validate(minimum)

        # Set the maximum value.
        if maximum is not None:
            self._max = self._validate(maximum)

        # Min and max.
        if minimum is not None and maximum is not None:
            if maximum < minimum:
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

           Keyword arguments:

           help    -- The help string.
           default -- The default value.
           allowed -- A list of allowed values.
        """

        # Call the base class constructor.
        super().__init__(help=help)

        # Set the allowed values.
        if allowed is not None:
            if type(allowed) is not list:
                allowed = [allowed]
            self._allowed = [self._validate(x) for x in allowed]

        # Set the default value.
        if default is not None:
            # Make sure the default is the right type.
            self._default = self._validate(default)
            # Make sure the default satisfies the constraints.
            self._validate_default()
            # Flag the the requirement is optional.
            self._is_optional = True

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

           Keyword arguments:

           help     -- The help string.
           optional -- Whether the file is optional.
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
        if not os.path.isfile(file):
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

           Keyword arguments:

           help     -- The help string.
           optional -- Whether the requirement is optional.
        """

        # Call the base class constructor.
        super().__init__(help=help, optional=optional)

    def getValue(self):
        """Return the value."""
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
                value = re.sub(r"\s+", "", value[0]).split(',')

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
                    if not os.path.isfile(file):
                        raise IOError("File doesn't exist: '%s'" % file)

        if len(uncompressed_files) > 0:
            return uncompressed_files
        else:
            return value

def _unarchive(name):
    """Decompress an archive and return a list of files."""

    # Get the directory name.
    dir = os.path.dirname(name)

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
            with tarfile.open(name) as tar:
                # We need to call tar.list(), otherwise the tar object will not know
                # about nested directories, i.e. it will appear as if ther is a single
                # member.
                print("Decompressing...")
                tar.list(verbose=False)

                # Loop over all of the members and get the file names.
                # If the name has no extension, then we assume that it's a directory.
                for file in tar.members:
                    if os.path.splitext(file.name)[1] is not "":
                        files.append(dir + file.name)

                # Now extract all of the files.
                tar.extractall(dir)

            return files

    # Get the file name and extension.
    file, ext = os.path.splitext(name)

    # This is a zip file.
    if ext.lower() == ".zip":
        # The list of decompressed files.
        files = None

        with zipfile.ZipFile(name) as zip:
            files = zip.namelist()
            print("Decompressing...")
            for file in files:
                print(file)
            zip.extractall(dir)

        return files

    # This is a gzip file.
    if ext.lower() == ".gz" or ext == ".gzip":
        with gzip.open(name, "rb") as f_in:
            with open(file, "wb") as f_out:
                print("Decompressing...\n%s" % name)
                shutil.copyfileobj(f_in, f_out)

        return file

    # This is a bzip2 file.
    if ext.lower() == ".bz2" or ext == ".bzip2":
        with bz2.open(name, "rb") as f_in:
            with open(file, "wb") as f_out:
                print("Decompressing...\n%s" % name)
                shutil.copyfileobj(f_in, f_out)

        return file

    # If we get this far, then this is not a supported archive.
    return None
