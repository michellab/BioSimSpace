"""Tools for reading/writing molecular systems."""

from Sire.Base import wrap as _wrap
from Sire.IO import MoleculeParser as _MoleculeParser

from collections import OrderedDict as _OrderedDict
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
            _formats_dict[format.upper()] = (extensions, description)

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
        return _formats_dict[format.upper()][1]
    except KeyError:
        print("Unsupported format: '%s'" % format)
        return None

def readMolecules(files):
    """Read a molecular system from file.

       Positional arguments:

       files -- A file name, or a list of file names.
    """

    return _MoleculeParser.read(files)

def saveMolecules(filebase, system, fileformat):
    """Save a molecular system to file.

       Positional arguments:

       filebase   -- The base name of the output file.
       system     -- The molecular system.
       fileformat -- The file format to save as.
    """

    return _MoleculeParser.save(system, filebase, \
                      {"fileformat":_wrap(fileformat)})
