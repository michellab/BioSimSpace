######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2022
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
__email__  = "lester.hedges@gmail.com"

__all__ = ["fileFormats",
           "formatInfo",
           "readMolecules",
           "readPDB",
           "readPerturbableSystem",
           "saveMolecules",
           "savePerturbableSystem"]

from collections import OrderedDict as _OrderedDict
from glob import glob as _glob
from io import StringIO as _StringIO

import os as _os
import shlex as _shlex
import sys as _sys
import subprocess as _subprocess
import tempfile as _tempfile
import warnings as _warnings

# Wrap the import of PyPDB since it imports Matplotlib, which will fail if
# we don't have a display running.
try:
    import pypdb as _pypdb
    _has_pypdb = True
except:
    _has_pypdb = False

# Flag that we've not yet raised a warning about GROMACS not being installed.
_has_gmx_warned = False

from sire.legacy import Base as _SireBase
from sire.legacy import IO as _SireIO
from sire.legacy import Mol as _SireMol
from sire.legacy import System as _SireSystem

from .. import _amber_home
from .. import _gmx_path
from .. import _isVerbose
from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from .._SireWrappers import Molecule as _Molecule
from .._SireWrappers import Molecules as _Molecules
from .._SireWrappers import System as _System
from .. import _Utils

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

def readPDB(id, pdb4amber=False, work_dir=None, property_map={}):
    """Read a molecular system from a Protein Data Bank (PDBP) ID in the RSCB PDB
       website.

       Parameters
       ----------

       id : str
           The PDB ID string, or path to a PDB file.

       pdb4amber : bool
           Whether to process the PDB file using pdb4amber. This reformats the file
           such that it can be handled by the AMBER suite of tools.

       work_dir : str
           The working directory used to run pdb4amber.

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
       >>> system = BSS.IO.readPDB("1a3n")

       Create a molecular system from a PDB file on disk and re-format so that
       it is compatible with the AmberTools suite.
       Data Bank (PDB) record.

       >>> import BioSimSpace as BSS
       >>> system = BSS.IO.readPDB("file.pdb", pdb4amber=True)
    """

    if not _has_pypdb:
        _warnings.warn("BioSimSpace.IO: PyPDB could not be imported on this system.")
        return None

    if not isinstance(id, str):
        raise TypeError("'id' must be of type 'str'")

    if not isinstance(pdb4amber, bool):
        raise TypeError("'pdb4amber' must be of type 'bool'")

    if work_dir and not isinstance(work_dir, str):
        raise TypeError("'work_dir' must be of type 'str'")

    # Create a temporary working directory and store the directory name.
    if work_dir is None:
        tmp_dir = _tempfile.TemporaryDirectory()
        work_dir = tmp_dir.name

    # User specified working directory.
    else:
        # Use full path.
        if work_dir[0] != "/":
            work_dir = _os.getcwd() + "/" + work_dir

        # Create the directory if it doesn't already exist.
        if not _os.path.isdir(work_dir):
            _os.makedirs(work_dir, exist_ok=True)

    # Path to a PDB file.
    if _os.path.isfile(id):
        pdb_file = _os.path.abspath(id)

    # ID from the Protein Data Bank.
    else:
        if not _has_pypdb:
            _warnings.warn("BioSimSpace.IO: PyPDB could not be imported on this system.")
            return None

        # Strip any whitespace from the PDB ID and convert to upper case.
        id = id.replace(" ", "").upper()

        # Attempt to download the PDB file. (Compression is currently broken!)
        with _warnings.catch_warnings(record=True) as w:
            pdb_string = _pypdb.get_pdb_file(id, filetype="pdb", compression=False)
            if w:
                raise IOError("Retrieval failed, invalid PDB ID: %s" % id)

        # Create the name of the PDB file.
        pdb_file = "%s/%s.pdb" % (work_dir, id)

        # Now write the PDB string to file.
        with open(pdb_file, "w") as file:
            file.write(pdb_string)

        # Store the absolute path of the file.
        pdb_file = _os.path.abspath(pdb_file)

    # Process the file with pdb4amber.
    if pdb4amber:
        # Check that pdb4amber exists.
        if _amber_home is None:
            raise _MissingSoftwareError("Please install AmberTools for pdb4amber support: http://ambermd.org")
        else:
            _pdb4amber_exe = "%s/bin/pdb4amber" % _amber_home
            if not _os.path.isfile(_pdb4amber_exe):
                raise IOError("Missing pdb4amber executable: '%s'" % _pdb4amber_exe)

                # Create the file prefix.
        prefix = work_dir + "/"

        # Create the pdb4amber command.
        command = "%s %s -o pdb4amber.pdb" % (_pdb4amber_exe, pdb_file)

        # Create files for stdout/stderr.
        stdout = open(prefix + "pdb4amber.out", "w")
        stderr = open(prefix + "pdb4amber.err", "w")

        # Run pdb4amber as a subprocess.
        proc = _subprocess.run(_Utils.command_split(command), cwd=work_dir,
            shell=False, stdout=stdout, stderr=stderr)
        stdout.close()
        stderr.close()

        # Check that the output PDB file was generated.
        # the expected output was generated.
        if _os.path.isfile("%s/pdb4amber.pdb" % work_dir):
            pdb_file = "%s/pdb4amber.pdb" % work_dir
        else:
            raise IOError("pdb4amber failed!")

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
    if _gmx_path is None and not _has_gmx_warned:
        _warnings.warn("BioSimSpace.IO: Please install GROMACS (http://www.gromacs.org) "
                       "for GROMACS topology file support.")
        _has_gmx_warned = True

    # Glob string to catch wildcards and convert to list.
    if isinstance(files, str):
        files = _glob(files)

    # Check that all arguments are of type 'str'.
    if isinstance(files, (list, tuple)):
        if not all(isinstance(x, str) for x in files):
            raise TypeError("'files' must be a list of 'str' types.")
        if len(files) == 0:
            raise ValueError("The list of input files is empty!")
    else:
        raise TypeError("'files' must be of type 'str', or a list of 'str' types.")

    # Validate the map.
    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'")

    # Add the GROMACS topology file path.
    if _gmx_path is not None and ("GROMACS_PATH" not in property_map):
        property_map["GROMACS_PATH"] = _gmx_path

    # Check that the files exist.
    for file in files:
        if not _os.path.isfile(file):
            raise IOError("Missing input file: '%s'" % file)

    # Try to read the files and return a molecular system.
    try:
        system = _SireIO.MoleculeParser.read(files, property_map)
    except Exception as e0:
        if "There are no lead parsers!" in str(e0):
            # First check to see if the failure was due to the presence
            # of CMAP records.
            for file in files:
                try:
                    amber_prm = _SireIO.AmberPrm(file)
                except Exception as e1:
                    if "CMAP" in str(e1).upper():
                        msg = ("Unable to parse AMBER topology file. CMAP "
                              "records are currently unsupported.")
                        if _isVerbose():
                            raise IOError(msg) from e1
                        else:
                            raise IOError(msg) from None
                    elif "CHAMBER" in str(e1).upper():
                        msg = ("Unable to parse AMBER topology file. "
                              "CHAMBER files are currently unsupported.")
                        if _isVerbose():
                            raise IOError(msg) from e1
                        else:
                            raise IOError(msg) from None

            msg = ("Failed to read molecules from %s. "
                   "It looks like you failed to include a topology file."
                  ) % files
            if _isVerbose():
                raise IOError(msg) from e0
            else:
                raise IOError(msg) from None
        else:
            if "Incompatibility" in str(e0):
                msg = "Incompatibility between molecular information in files: %s" % files
                if _isVerbose():
                    raise IOError(msg) from e0
                else:
                    raise IOError(msg) from None
            else:
                msg = "Failed to read molecules from: %s" % files
                if _isVerbose():
                    raise IOError(msg) from e0
                else:
                    raise IOError(msg) from None

    return _System(system)

def saveMolecules(filebase, system, fileformat, property_map={}):
    """Save a molecular system to file.

       Parameters
       ----------

       filebase : str
           The base name of the output files.

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
    if _gmx_path is None and not _has_gmx_warned:
        _warnings.warn("BioSimSpace.IO: Please install GROMACS (http://www.gromacs.org) "
                       "for GROMACS topology file support.")
        _has_gmx_warned = True

    # Check that the filebase is a string.
    if not isinstance(filebase, str):
        raise TypeError("'filebase' must be of type 'str'")

    # Check that that the system is of the correct type.

    # A System object.
    if isinstance(system, _System):
        pass
    # A Molecule object.
    elif isinstance(system, _Molecule):
        system = _System(system)
    elif isinstance(system, _Molecules):
        system = system.toSystem()
    # A list of Molecule objects.
    elif isinstance(system, list) and all(isinstance(x, _Molecule) for x in system):
        system = _System(system)
    # Invalid type.
    else:
        raise TypeError("'system' must be of type 'BioSimSpace.SireWrappers.System', "
                        "'BioSimSpace._SireWrappers.Molecule, 'BioSimSpace._SireWrappers.Molecules' "
                        "or a list of 'BiSimSpace._SireWrappers.Molecule' types.")

    # Check that fileformat argument is of the correct type.

    # Convert to a list if a single string is passed.
    # We split on ',' since the user might pass system.fileFormats() as the argument.
    if isinstance(fileformat, str):
        fileformat = fileformat.split(",")
    # Lists and tuples are okay!
    elif isinstance(fileformat, (list, tuple)):
        pass
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
    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'")

    # Copy the map.
    _property_map = property_map.copy()

    # Add the GROMACS topology file path.
    if _gmx_path is not None and ("GROMACS_PATH" not in _property_map):
        _property_map["GROMACS_PATH"] = _gmx_path

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

        # Warn the user if any molecules are parameterised with a force field
        # that uses geometric combining rules. While we can write this to file
        # the information is lost on read.
        if format == "PRM7":

            # Get the name of the "forcefield" property.
            forcefield = _property_map.get("forcefield", "forcefield")

            # Loop over all molecules in the system.
            for mol in system.getMolecules():
                if mol._sire_object.hasProperty(forcefield):
                    if mol._sire_object.property(forcefield).combiningRules() == "geometric":
                        _warnings.warn("AMBER topology files do not support force fields that "
                                       "use geometric combining rules, as this cannot be specified "
                                       "in the file. When this file is re-read, then arithmetic "
                                       "combining rules will be assumed.")
                        # Exit after the first non-arithmetic molecule we encounter.
                        break

        # Write the file.
        try:
            # Add CONECT record for single molecule PDB files.
            if format == "PDB" and system.nMolecules() == 1:
                # Generate a PDB parser object.
                pdb = _SireIO.PDB2(system._sire_object, _property_map)

                # Get the lines.
                lines = pdb.toLines()

                # For molecules that were perturbable, we must create the
                # CONECT record from the bonding information.
                if system[0]._sire_object.hasProperty("is_perturbable") or \
                   system[0]._sire_object.hasProperty("was_perturbable"):
                    bond = _property_map.get("bond", "bond0")
                    conect = _bond_to_conect(system[0]._sire_object.property(bond),
                                             system[0]._sire_object.info())

                # Create a connectivty object and generate the CONECT record.
                else:
                    conect = _SireMol.Connectivity(system[0]._sire_object,
                                                   _SireMol.CovalentBondHunter(),
                                                   _property_map).toCONECT()

                # Create the updated PDB file.
                pdb_records = "\n".join(lines[:-2]) \
                            + "\n" + conect + "\n"  \
                            + "\n".join(lines[-2:])

                # Write the default file to get the full path.
                file = _SireIO.MoleculeParser.save(system._sire_object, filebase, _property_map)

                # Now overwrite the file the PDB file with the updated records.
                with open(file[0], "w") as pdb_file:
                    pdb_file.write(pdb_records)

            else:
                # Make sure AMBER and GROMACS files have the expected water topology
                # and save GROMACS files with an extension such that they can be run
                # directly by GROMACS without needing to be renamed.
                if format == "PRM7" or format == "RST7":
                    system = system.copy()
                    system._set_water_topology("AMBER", _property_map)
                    file = _SireIO.MoleculeParser.save(system._sire_object, filebase, _property_map)
                elif format == "GroTop":
                    system = system.copy()
                    system._set_water_topology("GROMACS")
                    file = _SireIO.MoleculeParser.save(system._sire_object, filebase, _property_map)[0]
                    new_file = file.replace("grotop", "top")
                    _os.rename(file, new_file)
                    file = [new_file]
                elif format == "Gro87":
                    # Write to 3dp by default, unless greater precision is
                    # requested by the user.
                    if "precision" not in _property_map:
                        _property_map["precision"] = _SireBase.wrap(3)
                    file = _SireIO.MoleculeParser.save(system._sire_object, filebase, _property_map)[0]
                    new_file = file.replace("gro87", "gro")
                    _os.rename(file, new_file)
                    file = [new_file]
                else:
                    file = _SireIO.MoleculeParser.save(system._sire_object, filebase, _property_map)

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

def savePerturbableSystem(filebase, system, property_map={}):
    """Save a system containing a perturbable molecule. This will be written in
       AMBER format, with a topology file for each end state of the perturbation,
       i.e. 'filebase0.prm7' and 'filebase1.prm7'. Coordinates for the end
       states are written to 'filebase0.rst7' and 'filebase1.rst7'.

       Parameters
       ----------

       filebase : str
           The base name of the output files.

       system : :class:`System <BioSimSpace._SireWrappers.System>`
           The molecular system.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }
    """

    # Check that the filebase is a string.
    if not isinstance(filebase, str):
        raise TypeError("'filebase' must be of type 'str'")

    # Check that the system is valid.
    if isinstance(system, _System):
        pass
    # A Molecule object.
    elif isinstance(system, _Molecule):
        system = _System(system)
    # A Molecules object.
    elif isinstance(system, _Molecules):
        system = system.toSystem()
    # A list of Molecule objects.
    elif isinstance(system, list) and all(isinstance(x, _Molecule) for x in system):
        system = _System(system)
    # Invalid type.
    else:
        raise TypeError("'system' must be of type 'BioSimSpace.SireWrappers.System', "
                        "'BioSimSpace._SireWrappers.Molecule, 'BioSimSpace._SireWrappers.Molecules' "
                        "or a list of 'BiSimSpace._SireWrappers.Molecule' types.")

    # Validate the map.
    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'")

    # Validate that there is a single perturbable molecule in the system.
    pert_mols = system.getPerturbableMolecules()
    if len(pert_mols) != 1:
        raise ValueError("The 'system' must contain a single perturbable molecule. "
                        f"Found {len(pert_mols)}!")

    # Extract the molecule
    pert_mol = pert_mols[0]

    # Create a copy of the system for the lambda=0 and lambda=1 end states.
    system0 = system.copy()
    system1 = system.copy()

    # Update the perturbable molecule in each system.
    system0.updateMolecules(pert_mol._toRegularMolecule(property_map=property_map,
                                                        is_lambda1=False))
    system1.updateMolecules(pert_mol._toRegularMolecule(property_map=property_map,
                                                        is_lambda1=True))

    # Save the topology files.
    saveMolecules(filebase + "0", system0, "prm7")
    saveMolecules(filebase + "1", system1, "prm7")

    # Save the coordinate files.
    saveMolecules(filebase + "0", system0, "rst7")
    saveMolecules(filebase + "1", system1, "rst7")

def readPerturbableSystem(top0, coords0, top1, coords1, property_map={}):
    """Read a perturbable system from file.

       Parameters
       ----------

       top0 : str
           The path to the topology file for the lambda=0 end state.

       coords0 : str
           The path to the coordinate file for the lambda=0 end state.

       top1 : str
           The path to the topology file for the lambda=1 end state.

       coords1 : str
           The path to the coordinate file for the lambda=1 end state.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       system : :class:`System <BioSimSpace._SireWrappers.System>`
           A molecular system.
    """

    if not isinstance(top0, str):
        raise TypeError("'top0' must be of type 'str'.")

    if not isinstance(coords0, str):
        raise TypeError("'coords0' must be of type 'str'.")

    if not isinstance(top1, str):
        raise TypeError("'top1' must be of type 'str'.")

    if not isinstance(coords1, str):
        raise TypeError("'coords1' must be of type 'str'.")

    # Check that the coordinate and topology files can be parsed.

    # lamba = 0 coordinates.
    try:
        _SireIO.AmberRst7(coords0)
    except Exception as e:
        msg = f"Unable to read lambda=0 coordinate file: {coords0}"
        if _isVerbose():
            raise IOError(msg) from e
        else:
            raise IOError(msg) from None

    # lamba = 1 coordinates.
    try:
        _SireIO.AmberRst7(coords1)
    except Exception as e:
        msg = f"Unable to read lambda=1 coordinate file: {coords1}"
        if _isVerbose():
            raise IOError(msg) from e
        else:
            raise IOError(msg) from None

    # lamba = 0 topology.
    try:
        parser = _SireIO.AmberPrm(top0)
    except Exception as e:
        msg = f"Unable to read lambda=0 topology file: {top0}"
        if _isVerbose():
            raise IOError(msg) from e
        else:
            raise IOError(msg) from None
    if parser.isEmpty():
        raise ValueError(f"Unable to read topology file for lamba=0 end state: {top0}")

    # lamba = 1 topology.
    try:
        parser = _SireIO.AmberPrm(top1)
    except Exception as e:
        msg = f"Unable to read lambda=1 topology file: {top1}"
        if _isVerbose():
            raise IOError(msg) from e
        else:
            raise IOError(msg) from None
    if parser.isEmpty():
        raise ValueError(f"Unable to read topology file for lamba=1 end state: {top1}")

    # Try loading the two end states.
    system0 = readMolecules([coords0, top0], property_map=property_map)
    system1 = readMolecules([coords1, top1], property_map=property_map)

    # Make sure the systems have the same number of molecules.
    if system0.nMolecules() != system1.nMolecules():
        raise ValueError("The two topologies contain a different number of molecules!")

    # Now loop through the molecules in each system to work out which
    # is the perturbable molecule. This will differ in the  'ambertype'
    # 'LJ' or 'charge' property.
    ambertype = property_map.get("ambertype", "ambertype")
    LJ = property_map.get("LJ", "LJ")
    charge = property_map.get("charge", "charge")
    has_pert = False
    for idx, (mol0, mol1) in enumerate(zip(system0.getMolecules(), system1.getMolecules())):
        for atom0, atom1 in zip(mol0.getAtoms(), mol1.getAtoms()):
            if atom0._sire_object.property(ambertype) != atom1._sire_object.property(ambertype) or \
               atom0._sire_object.property(LJ)        != atom1._sire_object.property(LJ)        or \
               atom0._sire_object.property(charge)    != atom1._sire_object.property(charge):
                   has_pert = True
                   break
        if has_pert:
            break

    if not has_pert:
        raise ValueError("No perturbable molecule was found?")

    # Extract the perturbable molecule.
    pert_mol = system0[idx]

    # Extract and copy the Sire molecule.
    mol = pert_mol._sire_object.__deepcopy__()

    # Make the molecule editable.
    mol = mol.edit()

    # Rename all properties in the molecule for the lambda=0 end state,
    # e.g.: "prop" --> "prop0". Then delete all properties named "prop"
    # and "prop1".
    for prop in mol.propertyKeys():
        # See if this property exists in the user map.
        new_prop = property_map.get(prop, prop) + "0"

        # Copy the property using the updated name.
        mol = mol.setProperty(new_prop, mol.property(prop)).molecule()

        # Delete the redundant property.
        mol = mol.removeProperty(prop).molecule()

    # Now add the properties for the lambda=1 end state.
    mol1 = system1[idx]._sire_object
    for prop in mol1.propertyKeys():
        # See if this property exists in the user map.
        new_prop = property_map.get(prop, prop) + "1"

        # Copy the property using the updated name.
        mol = mol.setProperty(new_prop, mol1.property(prop)).molecule()

    # Flag that the molecule is perturbable.
    mol.setProperty("is_perturbable", _SireBase.wrap(True))

    # Add the molecule0 and molecule1 properties.
    mol.setProperty("molecule0", system0[idx]._sire_object)
    mol.setProperty("molecule1", system1[idx]._sire_object)

    # Commit the changes.
    mol = _Molecule(mol.commit())

    # Update the molecule in the original system.
    system0.updateMolecules(mol)

    return system0

def _bond_to_conect(bonds, info):
    """Helper function to create PDB CONECT records from bonding
       information from a molecular potential.

       Parameters
       ----------

       bonds : Sire.MM._MM.TwoAtomFunctions
           The bond potential property.

       info : Sire.Mol._Mol.MoleculeInfo
           The molecule info object.

       Returns
       -------

       conect : str
           The PDB CONECT records.
    """

    # Initialise the CONECT dictionary.
    conect_dict = {}

    # Loop over each potential term and the bonded atoms to the dictionary,
    # remembering that PDB atoms are one-indexed.
    for b in bonds.potentials():
        idx0 = info.atomIdx(b.atom0()).value() + 1
        idx1 = info.atomIdx(b.atom1()).value() + 1
        if idx0 in conect_dict:
            conect_dict[idx0].append(idx1)
        else:
            conect_dict[idx0] = [idx1]

    # Now create the CONECT record string.

    # Initialise the string.
    conect = ""

    # Sort keys in numerical order.
    k = list(conect_dict.keys())
    k.sort()

    # Store the number of keys.
    num_k = len(k)

    # Add a record for each root index.
    for x, idx0 in enumerate(k):
        # Add bonds in numerical order.
        v = conect_dict[idx0]
        v.sort()
        conect += f"CONECT{idx0:>5}"
        for idx1 in v:
            conect += f"{idx1:>5}"
        if x < num_k - 1:
            conect += "\n"

    return conect
