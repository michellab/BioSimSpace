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

import Sire.Base as _SireBase
import Sire.IO as _SireIO
import Sire.Mol as _SireMol
import Sire.System as _SireSystem

from .._SireWrappers import Molecule as _Molecule
from .._SireWrappers import System as _System

from collections import OrderedDict as _OrderedDict
from glob import glob
from io import StringIO as _StringIO

import os.path as _path
import subprocess as _subprocess
import sys as _sys

# Set the bundled GROMACS topology file directory.
_gromacs_path = _path.dirname(_SireBase.getBinDir()) + "/share/gromacs/top"

# The directory is missing. GROMACS must not be installed.
if not _path.isdir(_gromacs_path):
    print("Missing GROMACS topology file directory: '%s'" % _gromacs_path)

    # Attempt to install GROMACS.
    print("Trying to install GROMACS.")
    command = "%s/conda install -y -q -c bioconda gromacs" % _SireBase.getBinDir()
    proc = _subprocess.run(command, shell=True, stdout=_subprocess.PIPE)

    # The installation failed.
    if proc.returncode != 0:
        raise RuntimeError("GROMACS installation failed: '%s'" % command)

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
    print(r"%s" % _SireIO.MoleculeParser.supportedFormats())

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

def readMolecules(files, map={}):
    """Read a molecular system from file.

       Positional arguments:

       files -- A file name, or a list of file names.

       Keyword arguments:

       map   -- A dictionary that maps system "properties" to their user defined
                values. This allows the user to refer to properties with their
                own naming scheme, e.g. { "charge" : "my-charge" }
    """

    # Add the GROMACS topology file path.
    if "GROMACS_PATH" not in map:
        map["GROMACS_PATH"] = _gromacs_path

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
        system = _SireIO.MoleculeParser.read(files, map)
    except:
        raise IOError("Failed to read molecules from: %s" % files)

    return _System(system)

def saveMolecules(filebase, system, fileformat, map={}):
    """Save a molecular system to file.

       Positional arguments:

       filebase   -- The base name of the output file.
       system     -- The molecular system.
       fileformat -- The file format (or formats) to save to.

       Keyword arguments:

       map   -- A dictionary that maps system "properties" to their user defined
                values. This allows the user to refer to properties with their
                own naming scheme, e.g. { "charge" : "my-charge" }
    """

    # Add the GROMACS topology file path.
    if "GROMACS_PATH" not in map:
        map["GROMACS_PATH"] = _gromacs_path

    # Check that the filebase is a string.
    if type(filebase) is not str:
        raise TypeError("'filebase' must be of type 'str'")

    # Check that that the system is of the correct type.

    # A Mystem object.
    if type(system) is _System:
        pass
    # A Molecule object.
    elif type(system) is _Molecule:
        system = [system]
    # A list of Molecule objects.
    elif type(system) is list and all(isinstance(x, _Molecule) for x in system):
        pass
    # Invalid type.
    else:
        raise TypeError("'system' must be of type 'BioSimSpace.System', "
            + "'BioSimSpace.Molecule, or a list of 'BiSimSpace.Molecule' types.")

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

    # We have a list of molecules. Create a new system and add each molecule.
    if type(system) is list:

        # Create a Sire system and molecule group.
        s = _SireSystem.System("BioSimSpace System")
        m = _SireMol.MoleculeGroup("all")

        # Add all of the molecules to the group.
        for molecule in system:
            m.add(molecule._getSireMolecule())

        # Add the molecule group to the system.
        s.add(m)

        # Wrap the system.
        system = _System(s)

    # A list of the files that have been written.
    files = []

    # Save the system using each file format.
    for format in formats:
        # Add the file format to the property map.
        _map = map
        map["fileformat"] = _SireBase.wrap(format)

        # Write the file.
        file = _SireIO.MoleculeParser.save(system._getSireSystem(), filebase, map)
        files += file

    return files
