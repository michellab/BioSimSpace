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
Functionality for reading/writing molecular systems.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

from Sire.Base import wrap as _wrap
from Sire.IO import MoleculeParser as _MoleculeParser

from .._SireWrappers import System as _System

from collections import OrderedDict as _OrderedDict
from glob import glob
from io import StringIO as _StringIO
import sys as _sys

# Context manager for capturing stdout.
# Taken from:
# https://stackoverflow.com/questions/16571150/how-to-capture-stdout-output-from-a-python-function-call
class _Capturing(list):
    def __enter__(self):
        self._stdout = _sys.stdout
        _sys.stdout = self._stringio = _StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio
        _sys.stdout = self._stdout

# Capture the supported format information
with _Capturing() as format_info:
    print(r"%s" % _MoleculeParser.supportedFormats())

# Create a list of the supported formats.
_formats = []

# Create a dictionary of format-description key-value pairs.
_formats_dict = _OrderedDict()

# Loop over the format information to populate the dictionary.
for index, line in enumerate(format_info):
    if "Parser" in line:
        format = line.split()[2]
        extensions = format_info[index+1]
        description = format_info[index+2]

        if format != "SUPPLEMENTARY":
            _formats.append(format)
            _formats_dict[format.replace(" ", "").upper()] = (format, description)

# Delete the redundant variables.
del format_info, index, line, format, extensions, description

def fileFormats():
    """Return a list of the supported formats."""
    return _formats

def formatInfo(format):
    """Return information for the specified file format.

       Positional arguments:

       format -- The file format.
    """

    try:
        return _formats_dict[format.replace(" ", "").upper()][1]
    except KeyError:
        print("Unsupported format: '%s'" % format)
        return None

def readMolecules(files):
    """Read a molecular system from file.

       Positional arguments:

       files -- A file name, or a list of file names.
    """

    # Convert to a list.
    if type(files) is str:
        files = [files]

    # Check that all arguments are of type 'str'.
    if type(files) is list:
        if not all(isinstance(x, str) for x in files):
            raise TypeError("'files' must be a list of 'str' types.")
        if len(files) == 0:
            raise ValueError("The list of input files is empty!")
    else:
        raise TypeError("'files' must be of type 'str', or a list of 'str' types.")

    # Try to read the files and return a molecular system.
    try:
        system = _MoleculeParser.read(files)
    except:
        raise IOError("Failed to read molecules from: %s" % files)

    return _System(system)

def saveMolecules(filebase, system, fileformat):
    """Save a molecular system to file.

       Positional arguments:

       filebase   -- The base name of the output file.
       system     -- The molecular system.
       fileformat -- The file format (or formats) to save to.
    """

    # Check that the filebase is a string.
    if type(filebase) is not str:
        raise TypeError("'filebase' must be of type 'str'")

    # Check that that the system is of the correct type.
    if type(system) is not _System:
        raise TypeError("'system' must be of type 'BioSimSpace.System'")

    # Check that fileformat argument is of the correct type.

    # Convert to a list if a single string is passed.
    # We split on ',' since the user might pass system.fileFormat() as the argument.
    if type(fileformat) is str:
        fileformat = fileformat.split(",")
    # Lists and tuples are okay!
    elif type(fileformat) is list:
        pass
    elif type(fileformat) is tuple:
        pass
    # Invalid.
    else:
        raise TypeError("'fileformat' must be a 'str' or a 'list' of 'str' types.")

    # Make sure all items in list or tuple are strings.
    if not all(isinstance(x, str) for x in fileformat):
        raise TypeError("'fileformat' must be a 'str' or a 'list' of 'str' types.")

    # Make a list of the matched file formats.
    formats = []

    # Make sure that all of the formats are valid.
    for format in fileformat:
        try:
            f = _formats_dict[format.replace(" ", "").upper()][0]
            formats.append(f)
        except KeyError:
            raise ValueError("Unsupported file format '%s'. Supported formats "
                "are: %s." % (format, str(_formats)))

    # A list of the files that have been written.
    files = []

    # Save the system using each file format.
    for format in formats:
        file = _MoleculeParser.save(system._getSireSystem(), filebase, \
                {"fileformat":_wrap(format)})
        files += file

    return files
