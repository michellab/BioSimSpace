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
Functionality for running simulation processes.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

import Sire as _Sire

from BioSimSpace import _is_interactive, _is_notebook

from ..Protocol._protocol import Protocol as _Protocol
from .._SireWrappers import System as _System

import collections as _collections
import operator as _operator
import os as _os
import pygtail as _pygtail
import timeit as _timeit
import warnings as _warnings
import tempfile as _tempfile
import zipfile as _zipfile

if _is_notebook():
    from IPython.display import FileLink as _FileLink

__all__ = ["Process"]

class _MultiDict(dict):
    """A multi-valued dictionary."""
    def __setitem__(self, key, value):
        """Add the given value to the list of values for this key."""
        self.setdefault(key, []).append(value)

class Process():
    """Base class for running different biomolecular simulation processes."""

    def __init__(self, system, protocol, name=None, work_dir=None, seed=None):
        """Constructor.

           Positional arguments:

           system    -- The molecular system.
           protocol  -- The protocol for the process.

           Keyword arguments:

           name      -- The name of the process.
           work_dir  -- The working directory for the process.
           seed      -- A random number seed.
        """

	# Don't allow user to create an instance of this base class.
        if type(self) is Process:
            raise Exception("<Process> must be subclassed.")

        # Check that the system is valid.
        if type(system) is not _System:
            raise TypeError("'system' must be of type 'BioSimSpace.System'")

        # Check that the protocol is valid.
        if not isinstance(protocol, _Protocol):
            raise TypeError("'protocol' must be of type 'BioSimSpace.Protocol'")

        # Check that the working directory is valid.
        if work_dir is not None and type(work_dir) is not str:
            raise TypeError("'work_dir' must be of type 'str'")

        # Check that the seed is valid.
        if seed is not None and type(seed) is not int:
            raise TypeError("'seed' must be of type 'int'")

        # Set the process to None.
        self._process = None

        # Is the process running interactively? If so, don't block
        # when a get method is called.
        self._is_blocked = not _is_interactive()

        # Whether this process can generate trajectory data.
        # Even if a process can generate a trajectory, whether it does
        # will depend on the chosen protocol.
        self._has_trajectory = False

	# Copy the passed system, protocol, and process name.
        self._system = system._getSireSystem()
        self._protocol = protocol

        # Set the name
        if name is None:
            self._name = None
        else:
            self.setName(name)

        # Set the random number seed.
        if seed is None:
            self._is_seeded = False
            self._seed = None
        else:
            self._is_seeded = True
            self.setSeed(seed)

        # Set the timer and running time None.
        self._timer = None
        self._runtime = None

        # Set the command-line string to None
        self._command = None

        # Set the list of input files to None.
        self._input_files = None

        # Create a temporary working directory and store the directory name.
        if work_dir is None:
            self._tmp_dir = _tempfile.TemporaryDirectory()
            self._work_dir = self._tmp_dir.name

        # User specified working directory.
        else:
            self._work_dir = work_dir

            # Create the directory if it doesn't already exist.
            if not _os.path.isdir(work_dir):
                _os.makedirs(work_dir)

        # Files for redirection of stdout and stderr.
        self._stdout_file = "%s/%s.out" % (self._work_dir, name)
        self._stderr_file = "%s/%s.err" % (self._work_dir, name)

        # Create the files. This makes sure that the 'stdout' and 'stderr'
        # methods can be called when the files are empty.
        open(self._stdout_file, "a").close()
        open(self._stderr_file, "a").close()

        # Initialise lists to store the contents of stdout and stderr.
        self._stdout = []
        self._stderr = []

        # Clean up any existing offset files.

        stdout_offset = "%s.offset" % self._stdout_file
        stderr_offset = "%s.offset" % self._stderr_file

        if _os.path.isfile(stdout_offset):
            _os.remove(stdout_offset)

        if _os.path.isfile(stderr_offset):
            _os.remove(stderr_offset)

        # Initialise the configuration file string list.
        self._config = []

        # Initaliae the command-line argument dictionary.
        self._args = _collections.OrderedDict()

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.Process.%s: system=%s, protocol=%s, exe='%s', name='%s', work_dir='%s' seed=%s>" \
            % (self.__class__.__name__, str(_System(self._system)), self._protocol.__repr__(),
               self._exe, self._name, self._work_dir, self._seed)

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return "BioSimSpace.Process.%s(%s, %s, exe='%s', name='%s', work_dir='%s', seed=%s)" \
            % (self.__class__.__name__, str(_System(self._system)), self._protocol.__repr__(),
               self._exe, self._name, self._work_dir, self._seed)

    def run(self, system=None, protocol=None, autostart=True, restart=False):
        """Create and run a new process.

           Keyword arguments:

           system    -- The molecular system.
           protocol  -- The simulation protocol.
           autostart -- Whether to start the process automatically.
           restart   -- Whether to restart the simulation, i.e. use the original system.

           return    -- The new process object.
        """

        # Try to get the current system.
        if not restart:
            system = self.getSystem()

        # Use the existing system.
        if system is None:
            system = self._system

        # Check that the new system is valid.
        else:
            if type(system) is not _System:
                raise TypeError("'system' must be of type 'BioSimSpace.System'")

        # Use the existing protocol.
        if protocol is None:
            protocol = self._protocol

        # Check that the new protocol is valid.
        else:
            if not isinstance(protocol, _Protocol):
                raise TypeError("'protocol' must be of type 'BioSimSpace.Protocol'")

        # Create the new process.
        process = type(self)(system, protocol)

        # Return the new process object.
        if autostart:
            return process.start()
        else:
            return process

    def getPackageName(self):
        """Return the package name."""
        return self._package_name

    def getName(self):
        """Return the process name."""
        return self._name

    def setName(self, name):
        """Set the process name."""

        if type(name) is not str:
            raise TypeError("'name' must be of type 'str'")
        else:
            self._name = name

    def getSeed(self):
        """Return the random number seed."""
        return self._seed

    def setSeed(self, seed):
        """Set the random number seed."""

        if type(seed) is not int:
            _warnings.warn("The seed must be an integer. Disabling seeding.")
            self._seed = None
        else:
            self._seed = seed

    def wait(self, max_time=None):
        """Wait for the process to finish.

           Keyword arguments:

           max_time -- The maximimum time to wait (in minutes).
        """

        # The process isn't running.
        if not self.isRunning():
            return

        if not max_time is None:
            if max_time <= 0:
                _warnings.warn("Maximum running time must be greater than zero. Using default (60 mins).")
                max_time = 60

            # Convert the time to milliseconds.
            max_time = max_time * 60 * 1000

            # Wait for the desired amount of time.
            self._process.wait(max_time)

        else:
            # Wait for the process to finish.
            self._process.wait()

    def isRunning(self):
        """Return whether the process is running"""
        try:
            return self._process.isRunning()
        except AttributeError:
            return False

    def isError(self):
        """Return whether the process errored."""
        try:
            return self._process.isError()
        except AttributeError:
            print("The process hasn't been started!")
            return None

    def kill(self):
        """Kill the running process."""
        if not self._process is None and self._process.isRunning():
            self._process.kill()

    def stdout(self, n=10):
        """Print the last n lines of the stdout buffer.

           Keyword arguments:

           n -- The number of lines to print.
        """

        # Ensure that the number of lines is positive.
        if n < 0:
            raise ValueError("The number of lines must be positive!")

        # Append any new lines to the stdout list.
        for line in _pygtail.Pygtail(self._stdout_file):
            self._stdout.append(line.rstrip())

        # Get the current number of lines.
        num_lines = len(self._stdout)

        # Set the line from which to start printing.
        if num_lines < n:
            start = 0
        else:
            start = num_lines - n

        # Print the lines.
        for x in range(start, num_lines):
            print(self._stdout[x])

    def stderr(self, n=10):
        """Print the last n lines of the stderr buffer.

           Keyword arguments:

           n -- The number of lines to print.
        """

        # Ensure that the number of lines is positive.
        if n < 0:
            raise ValueError("The number of lines must be positive!")

        # Append any new lines to the stdout list.
        for line in _pygtail.Pygtail(self._stderr_file):
            self._stderr.append(line.rstrip())

        # Get the current number of lines.
        num_lines = len(self._stderr)

        # Set the line from which to start printing.
        if num_lines < n:
            start = 0
        else:
            start = num_lines - n

        # Print the lines.
        for x in range(start, num_lines):
            print(self._stderr[x])

    def exe(self):
        """Return the executable."""
        return self._exe

    def inputFiles(self):
        """Return the list of input files."""
        return self._input_files.copy()

    def workDir(self):
        """Return the working directory."""
        return self._work_dir

    def getStdout(self, block="AUTO"):
        """Return the entire stdout for the process as a list of strings.

           Keyword arguments:

           block -- Whether to block until the process has finished running.
        """

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        # Append any new lines to the stdout list.
        for line in _pygtail.Pygtail(self._stdout_file):
            self._stdout.append(line.rstrip())

        return self._stdout.copy()

    def getStderr(self, block="AUTO"):
        """Return the entire stderr for the process as a list of strings.

           Keyword arguments:

           block -- Whether to block until the process has finished running.
        """

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        # Append any new lines to the stdout list.
        for line in _pygtail.Pygtail(self._stderr_file):
            self._stderr.append(line.rstrip())

        return self._stderr.copy()

    def getOutput(self, name=None, block="AUTO"):
        """Return a link to a zip file of the working directory.

           Keyword arguments:

           name  -- The name of the zip file.
           block -- Whether to block until the process has finished running.
        """

        if name is None:
            name = self._name
        else:
            if type(name) is not str:
                raise TypeError("'name' must be of type 'str'")

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        # Generate the zip file name.
        zipname = "%s.zip" % name

        # Glob all of the output files.
        output = _glob.glob("%s/*" % self._work_dir)

        with _zipfile.ZipFile(zipname, "w") as zip:
            # Loop over all of the file outputs.
            for file in output:
                zip.write(file, arcname=_os.path.basename(file))


        # Return a link to the archive.
        if _is_notebook():
            return _FileLink(zipname)
        # Return the path to the archive.
        else:
            return zipname

    def command(self):
        """Return the command-line string used to run the process."""
        return self._command

    def getConfig(self):
        """Get the list of configuration file strings."""
        return self._config.copy()

    def setConfig(self, config):
        """Set the list of configuration file strings."""

        # Check that the passed configuration is a list of strings.
        if _is_list_of_strings(config):
            self._config = config
            self.writeConfig(self._config_file)

        # The user has passed a path to a file.
        elif _os.path.isfile(config):

            # Clear the existing config.
            self._config = []

            # Read the contents of the file.
            with open(config, "r") as f:
                for line in f:
                    self._config.append(line.rstrip())

            # Write the new configuration file.
            self.writeConfig(self._config_file)

        else:
            raise ValueError("'config' must be a list of strings, or a file path.")

    def addToConfig(self, config):
        """Add a string to the configuration list."""
        # Append a single string.
        if type(config) is str:
            self._config.append(config)
            self.writeConfig(self._config_file)

        # Extend the list with the additional strings.
        elif _is_list_of_strings(config):
            self._config.extend(config)
            self.writeConfig(self._config_file)

        # A path to a file.
        elif _os.path.isfile(config):

            # Read the contents of the file.
            with open(file, "r") as f:
                for line in f:
                    self._config.append(line)

            # Write the new configuration file.
            self.writeConfig(self._config_file)

        else:
            raise ValueError("'config' must be a string, list of strings, or a file path.")

    def resetConfig(self):
        """Reset the configuration parameters."""
        self._generate_config()

        # Write the new configuration file.
        self.writeConfig(self._config_file)

    def writeConfig(self, file):
        """Write the configuration to file."""
        with open(file, "w") as f:
            for line in self._config:
                f.write("%s\n" % line)

    def getArgs(self):
        """Get the dictionary of command-line arguments."""
        return self._args.copy()

    def getArgString(self):
        """Get the command-line arguments string."""
        return " ".join(self.getArgStringList())

    def getArgStringList(self):
        """Convert the argument dictionary into a list of strings."""

        # Create an empty list.
        args = []

        # Add the arguments to the list.
        for key, value in self._args.items():
            # Boolean flag.
            if type(value) is bool:
                if value:
                    args.append(str(key))
            else:
                args.append(str(key))
                args.append(str(value))

        return args

    def setArgs(self, args):
        """Set the dictionary of command-line arguments."""
        if isinstance(args, _collections.OrderedDict):
            self._args = args

        elif isinstance(args, dict):
            self._args = _collections.OrderedDict(args)

    def setArg(self, arg, value):
        """Set a specific command-line argument.

           Keyword arguments:

           arg   -- The argument to set.
           value -- The value of the argument.

           For command-line flags, i.e. boolean arguments, the key should
           specify whether the flag is enabled (True) or not (False).
        """
        self._args[arg] = value

    def insertArg(self, arg, value, index):
        """Insert a command-line argument at a specific index.

           Keyword arguments:

           arg   -- The argument to set.
           value -- The value of the argument.
           index -- The index in the dictionary.
        """
        _odict_insert(self._args, arg, value, index)

    def addArgs(self, args):
        """Append additional command-line arguments.

           Keyword arguments:

           args -- A dictionary of arguments.
        """
        if isinstance(args, dict) or isinstance(args, _collections.OrderedDict):
            for arg, value in args.items():
                self._args[arg] = value

    def deleteArg(self, arg):
        """Delete an argument from the dictionary."""
        try:
            del self._args[arg]

        except KeyError:
            pass

    def clearArgs(self):
        """Clear all of the command-line arguments."""
        self._args.clear()

    def resetArgs(self):
        """Reset the command-line arguments."""
        self._generate_args()

    def runTime(self):
        """Return the running time for the process (in minutes)."""

        # The process is still running.
        if self._process.isRunning():
            self._runtime = (_timeit.default_timer() - self._timer) / 60
            return self._runtime

        # The process has finished.
        else:
            # Return the runtime and reset the timer.
            if self._timer is not None:
                self._runtime = (_timeit.default_timer() - self._timer) / 60
                self._timer = None
                return self._runtime

            # The process has finished. Return the previous run time.
            else:
                return self._runtime

def _getAABox(system):
    """Get the axis-aligned bounding box for the molecular system.

       Keyword arguments:

       system -- A Sire molecular system.
    """

    # Initialise the coordinates vector.
    coord = []

    # Loop over all of the molecules.
    for n in system.molNums():

        # Extract the atomic coordinates and append them to the vector.
        try:
            coord.extend(system[n].property("coordinates").toVector())

        except UserWarning:
            raise

    # Return the AABox for the coordinates.
    return _Sire.Vol.AABox(coord)

def _get_box_size(system):
    """Get the size of the periodic box."""

    try:
        box = system.property("space")
        return box.dimensions()

    except UserWarning:
        return None

def _restrain_backbone(system):
    """Restrain protein backbone atoms.

        Keyword arguments:

        system -- A Sire molecular system.
    """

    # Copy the original system.
    s = system

    # A list of amino acid name abbreviations.
    # Since we only want to restrain atoms in protein backbones, we compare
    # molecule residue names against this list in order to determine whether
    # the molecule is a protein.
    amino_acids = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE",
        "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "SEC",
        "VAL", "TRP", "TYR"]

    # Loop over all molecules by number.
    for n in s.molNums():

        # Extract the molecule and make it editable.
        m = s.molecule(n).edit()

        # Protein flag.
        is_protein = False

        # Compare each residue name against the amino acid list.
        for res in m.residues():

            # This residue is an amino acid.
            if res.name().value().upper() in amino_acids:
                is_protein = True
                break

        # Restrain the protein backbone.
        if is_protein:

            # Select the backbone.
            backbone = m.atoms(_Sire.Mol.AtomName("CA", _Sire.ID.CaseInsensitive) *
                               _Sire.Mol.AtomName("N",  _Sire.ID.CaseInsensitive) *
                               _Sire.Mol.AtomName("C",  _Sire.ID.CaseInsensitive) *
                               _Sire.Mol.AtomName("O",  _Sire.ID.CaseInsensitive))

            # Set the restrained property for each atom in the backbone.
            for atom in backbone:
                m = m.atom(atom.index()).setProperty("restrained", 1.0).molecule()

            # Update the system.
            s.update(m.commit())

    # Return the new system.
    return s

    def _generate_args(self):
        """Generate the dictionary of command-line arguments."""
        self.clearArgs()

    def _get_trajectory_files(self):
        """Get all files associated with the molecular trajectory."""
        return None

def _is_list_of_strings(lst):
    """Check whether the passed argument is a list of strings."""
    if lst and isinstance(lst, list):
        return all(isinstance(elem, str) for elem in lst)
    else:
        return False

def _odict_insert(dct, key, value, index):
    """Insert an item into an ordered dictionary."""

    # Store the original size of the dictionary.
    n = len(dct)

    # Make sure the index is within range.
    if index < 0 or index > n+1:
        raise IndexError("Dictionary index out of range!")

    # Insert the new item at the end.
    dct[key] = value

    # Now loop over the original dict, moving any items
    # beyond 'index' to the end.

    # Index counter
    i = 0

    for item in list(dct):
        i += 1

        if i > n:
            break
        elif i > index:
            dct.move_to_end(item)
