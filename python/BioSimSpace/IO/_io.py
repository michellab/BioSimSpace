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
Functionality for reading/writing molecular systems.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["fileFormats", "formatInfo", "readMolecules", "readPDB", "saveMolecules"]

from collections import OrderedDict as _OrderedDict
from io import StringIO as _StringIO
from warnings import warn as _warn

import os as _os
import sys as _sys
import tempfile as _tempfile

# Wrap the import of PyPDB since it imports Matplotlib, which will fail if
# we don't have a display running.
try:
    import pypdb as _pypdb
    _has_pypdb = True
except:
    _has_pypdb = False

# Flag that we've not yet raised a warning about GROMACS not being installed.
_has_gmx_warned = False

from Sire import Base as _SireBase
from Sire import IO as _SireIO
from Sire import Mol as _SireMol
from Sire import System as _SireSystem

from BioSimSpace import _gromacs_path
from BioSimSpace import _isVerbose
from BioSimSpace._SireWrappers import Molecule as _Molecule
from BioSimSpace._SireWrappers import Molecules as _Molecules
from BioSimSpace._SireWrappers import System as _System

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

# Create a dictionary of format-description key:value pairs.
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
    """Return a list of the supported formats.

       Returns
       -------

       formats : [str]
           A list of the support file formats.
    """
    return _formats

def formatInfo(format):
    """Return information for the specified file format.

       Parameters
       ----------

       format : str
           The file format.

       Returns
       -------

       info : str
           A description of the named file format.

       Examples
       --------

       Display information regarding the PDB format.

       >>> import BioSimSpace as BSS
       >>> BSS.formatInfo("PDB")

       Print information for each of the supported file formats.

       >>> import BioSimSpace as BSS
       >>> for format in BSS.IO.fileFormats:
       ...     BSS.IO.formatInfo(format)

    """

    try:
        return _formats_dict[format.replace(" ", "").upper()][1]
    except KeyError:
        print("Unsupported format: '%s'" % format)
        return None

def readPDB(id, property_map={}):
    """Read a molecular system from a Protein Data Bank (PDBP) ID in the RSCB PDB
       website.

       Parameters
       ----------

       id : str
           The PDB ID string.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       system : :class:`System <BioSimSpace._SireWrappers.System>`
           A molecular system.

       Examples
       --------

       Create a molecular system from the deoxy human haemoglobin Protein
       Data Bank (PDB) record.

       >>> import BioSimSpace as BSS
       >>> system = BSS.readPDB("1a3n")
    """

    if not _has_pypdb:
        _warn("BioSimSpace.IO: PyPDB could not be imported on this system.")
        return None

    if type(id) is not str:
        raise TypeError("'id' must be of type 'str'")

    # Strip any whitespace from the PDB ID and convert to upper case.
    id = id.replace(" ", "").upper()

    # Create a temporary directory to write the PDB file.
    tmp_dir = _tempfile.TemporaryDirectory()

    # Attempt to download the PDB file. (Compression is currently broken!)
    try:
        pdb_string = _pypdb.get_pdb_file(id, filetype="pdb", compression=False)
    except:
        raise IOError("Invalid PDB ID: '%s'" % id)

    # Create the name of the PDB file.
    pdb_file = "%s/%s.pdb" % (tmp_dir.name, id)

    # Now write the PDB string to file.
    with open(pdb_file, "w") as file:
        file.write(pdb_string)

    # Read the file and return a molecular system.
    return readMolecules(pdb_file, property_map)

def readMolecules(files, property_map={}):
    """Read a molecular system from file.

       Parameters
       ----------

       files : str, [str]
           A file name, or a list of file names.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       system : :class:`System <BioSimSpace._SireWrappers.System>`
           A molecular system.

       Examples
       --------

       Load a molecular system from AMBER coordinate and topology files.

       >>> import BioSimSpace as BSS
       >>> system = BSS.IO.readMolecules(["ala.rst7", "ala.prm7"])

       Load the same system, but map the "charge" property to the key "my-charge".

       >>> import BioSimSpace as BSS
       >>> system = BSS.IO.readMolecules(["ala.rst7", "ala.prm7"], property_map={"charge" : "my-charge"})

       >>> import BioSimSpace as BSS
       >>> system = BSS.IO.readMolecules(["ala.rst7", "ala.prm7"])

       Load a molecular system from all of the files contained within a directory.

       >>> import BioSimSpace as BSS
       >>> system = BSS.IO.readMolecules(BSS.IO.glob("dir/*"))

       Load a molecular system from GROMACS coordinate and topology files using
       a custom GROMACS topology directory.

       >>> import BioSimSpace as BSS
       >>> system = BSS.IO.readMolecules(["mol.gro87", "mol.grotop"], property_map={"GROMACS_PATH" : "/path/to/gromacs/topology"})
    """

    global _has_gmx_warned
    if _gromacs_path is None and not _has_gmx_warned:
        _warn("BioSimSpace.IO: Please install GROMACS (http://www.gromacs.org) "
              "for GROMACS topology file support.")
        _has_gmx_warned = True

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

    # Validate the map.
    if type(property_map) is not dict:
        raise TypeError("'property_map' must be of type 'dict'")

    # Add the GROMACS topology file path.
    if _gromacs_path is not None and ("GROMACS_PATH" not in property_map):
        property_map["GROMACS_PATH"] = _gromacs_path

    # Check that the files exist.
    for file in files:
        if not _os.path.isfile(file):
            raise IOError("Missing input file: '%s'" % file)

    # Try to read the files and return a molecular system.
    try:
        system = _SireIO.MoleculeParser.read(files, property_map)
    except Exception as e:
        if "There are no lead parsers!" in str(e):
            msg = ("Failed to read molecules from %s. "
                   "It looks like you failed to include a topology file."
                  ) % files
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None
        else:
            if "Incompatibility" in str(e):
                msg = "Incompatibility between molecular information in files: %s" % files
                if _isVerbose():
                    raise IOError(msg) from e
                else:
                    raise IOError(msg) from None
            else:
                msg = "Failed to read molecules from: %s" % files
                if _isVerbose():
                    raise IOError(msg) from e
                else:
                    raise IOError(msg) from None

    return _System(system)

def saveMolecules(filebase, system, fileformat, property_map={}):
    """Save a molecular system to file.

       Parameters
       ----------

       filebase : str
           The base name of the output file.

       system : :class:`System <BioSimSpace._SireWrappers.System>`, \
                :class:`Molecule< BioSimSpace._SireWrappers.Molecule>` \
                :class:`Molecule< BioSimSpace._SireWrappers.Molecules>`
           The molecular system.

       fileformat : str, [str]
           The file format (or formats) to save to.

       property_map : dict
           A dictionary that maps system "properties" to their user
           defined values. This allows the user to refer to properties
           with their own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       files : [str]
           The list of files that were generated.

       Examples
       --------

       Load a molecular system from AMBER coordinate and topology files then
       try to save it to all supported file formats.

       >>> import BioSimSpace as BSS
       >>> system = BSS.IO.readMolecules(["ala.rst7", "ala.prm7"])
       >>> for format in BSS.IO.fileFormats():
       ...     try:
       ...         BSS.IO.saveMolecules("test", system, format)
       ...     except:
       ...         print("Could not convert to format: '%s'" % format)

       Load a molecular system from AMBER coordinate and topology files then
       try to save it to GROMACS format, mapping and un-mapping the charge
       property along the way.

       >>> import BioSimSpace as BSS
       >>> system = BSS.IO.readMolecules(["ala.rst7", "ala.prm7"], property_map={"charge" : "my-charge"})
       >>> BSS.IO.saveMolecules("test", system, ["gro87", "grotop"], property_map={"charge" : "my-charge"})
    """

    global _has_gmx_warned
    if _gromacs_path is None and not _has_gmx_warned:
        _warn("BioSimSpace.IO: Please install GROMACS (http://www.gromacs.org) "
              "for GROMACS topology file support.")
        _has_gmx_warned = True

    # Check that the filebase is a string.
    if type(filebase) is not str:
        raise TypeError("'filebase' must be of type 'str'")

    # Check that that the system is of the correct type.

    # A System object.
    if type(system) is _System:
        pass
    # A Molecule object.
    elif type(system) is _Molecule:
        system = _System(system)
    elif type(system) is _Molecules:
        system = system.toSystem()
    # A list of Molecule objects.
    elif type(system) is list and all(isinstance(x, _Molecule) for x in system):
        system = _System(system)
    # Invalid type.
    else:
        raise TypeError("'system' must be of type 'BioSimSpace.SireWrappers.System', "
                        "'BioSimSpace._SireWrappers.Molecule, 'BioSimSpace._SireWrappers.Molecules' "
                        "or a list of 'BiSimSpace._SireWrappers.Molecule' types.")

    # Check that fileformat argument is of the correct type.

    # Convert to a list if a single string is passed.
    # We split on ',' since the user might pass system.fileFormats() as the argument.
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

    # Validate the map.
    if type(property_map) is not dict:
        raise TypeError("'property_map' must be of type 'dict'")

    # Copy the map.
    _property_map = property_map.copy()

    # Add the GROMACS topology file path.
    if _gromacs_path is not None and ("GROMACS_PATH" not in _property_map):
        _property_map["GROMACS_PATH"] = _gromacs_path

    # Get the directory name.
    dirname = _os.path.dirname(filebase)

    # If the user has passed a directory, make sure that is exists.
    if _os.path.basename(filebase) != filebase:
        # Create the directory if it doesn't already exist.
        if not _os.path.isdir(dirname):
            _os.makedirs(dirname, exist_ok=True)

    # Store the current working directory.
    dir = _os.getcwd()

    # Change to the working directory for the process.
    # This avoid problems with relative paths.
    if dirname != "":
        _os.chdir(dirname)

    # A list of the files that have been written.
    files = []

    # Save the system using each file format.
    for format in formats:
        # Add the file format to the property map.
        _property_map["fileformat"] = _SireBase.wrap(format)

        # Write the file.
        try:
            file = _SireIO.MoleculeParser.save(system._getSireObject(), filebase, _property_map)
            files += file
        except Exception as e:
            if dirname != "":
                _os.chdir(dir)
            msg = "Failed to save system to format: '%s'" % format
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

    # Change back to the original directory.
    if dirname != "":
        _os.chdir(dir)

    # Return the list of files.
    return files
