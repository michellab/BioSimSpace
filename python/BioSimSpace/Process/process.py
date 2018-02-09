"""
@package biosimspace
@author  Lester Hedges
@brief   Base class and helper functions for the various sample modules.
"""

from Sire.ID import CaseInsensitive
from Sire.Mol import AtomName
from Sire.Vol import AABox

from ..Protocol.protocol import Protocol

from collections import OrderedDict
from operator import add, sub
from os import makedirs, path, remove
from timeit import default_timer as timer
from warnings import warn

import tempfile

try:
    from Sire import try_import
    pygtail = try_import("pygtail")
except ImportError:
    raise ImportError("Pygtail is not installed. Please install pygtail in order to use BioSimSpace.")

class MultiDict(dict):
    """A multi-valued dictionary."""
    def __setitem__(self, key, value):
        """Add the given value to the list of values for this key."""
        self.setdefault(key, []).append(value)

class Process():
    """Base class for running different biomolecular simulation processes."""

    def __init__(self, system, protocol, name="process", work_dir=None, seed=None):
        """Constructor.

           Keyword arguments:

           system    -- The molecular system.
           protocol  -- The protocol for the process.
           name      -- The name of the process.
           work_dir  -- The working directory for the process.
           seed      -- A random number seed.
        """

	# Don't allow user to create an instance of this base class.
        if type(self) == Process:
            raise Exception("<Process> must be subclassed.")

        # Set the process to None.
        self._process = None

	# Copy the passed system, protocol, and process name.
        self._system = system
        self._protocol = protocol
        self._name = name

        # Set the random number seed.
        if seed is None:
            self._is_seeded = False
            self.seed = 0
        else:
            self._is_seeded = True
            self.seed = seed

        # Set the timer and running time None.
        self._timer = None
        self._runtime = None

        # Set the command-line string to None
        self._command = None

        # Set the list of input files to None.
        self._input_files = None

        # Check whether this is a Protocol object.
        # If not, we assume that the user has passed a custom configuration file.
        if isinstance(protocol, Protocol):
            self._is_custom = False
        else:
            self._is_custom = True

        # Create a temporary working directory and store the directory name.
        if work_dir is None:
            self._tmp_dir = tempfile.TemporaryDirectory()
            self._work_dir = self._tmp_dir.name

        # User specified working directory.
        else:
            self._work_dir = work_dir

            # Create the directory if it doesn't already exist.
            if not path.isdir(work_dir):
                makedirs(work_dir)

        # Files for redirection of stdout and stderr.
        self._stdout_file = "%s/%s.out" % (self._work_dir, name)
        self._stderr_file = "%s/%s.err" % (self._work_dir, name)

        # Create the files. This makes sure that the 'stdout' and 'stderr'
        # methods can be called when the files are empty.
        open(self._stdout_file, 'a').close()
        open(self._stderr_file, 'a').close()

        # Initialise lists to store the contents of stdout and stderr.
        self._stdout = []
        self._stderr = []

        # Clean up any existing offset files.

        stdout_offset = "%s.offset" % self._stdout_file
        stderr_offset = "%s.offset" % self._stderr_file

        if path.isfile(stdout_offset):
            remove(stdout_offset)

        if path.isfile(stderr_offset):
            remove(stderr_offset)

        # Initialise the configuration file string list.
        self._config = []

        # Initaliae the command-line argument dictionary.
        self._args = OrderedDict()

    @property
    def seed(self):
        """Return the random number seed."""
        return self._seed

    @seed.setter
    def seed(self, seed):
        """Set the random number seed."""

        if type(seed) is not int:
            warn("The seed must be an integer. Disabling seeding.")
            self._seed = None

        else:
            self._seed = seed

    def wait(self, max_time=None):
        """Wait for the process to finish.

           Keyword arguments:

           max_time -- The maximimum time to wait (in minutes).
        """

        if not max_time is None:
            if max_time <= 0:
                warn("Maximum running time must be greater than zero. Using default (60 mins).")
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
        return self._process.isRunning()

    def isError(self):
        """Return whether the process errored."""
        return self._process.isError()

    def kill(self):
        """Kill the running process."""
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
        for line in pygtail.Pygtail(self._stdout_file):
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
        for line in pygtail.Pygtail(self._stderr_file):
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
        return self._input_files

    def workDir(self):
        """Return the working directory."""
        return self._work_dir

    def getOutput(self):
        """Return the entire stdout for the process as a list of strings."""

        # Append any new lines to the stdout list.
        for line in pygtail.Pygtail(self._stdout_file):
            self._stdout.append(line.rstrip())

        return self._stdout

    def getError(self):
        """Return the entire stderr for the process as a list of strings."""

        # Append any new lines to the stdout list.
        for line in pygtail.Pygtail(self._stderr_file):
            self._stderr.append(line.rstrip())

        return self._stderr

    def command(self):
        """Return the command-line string used to run the process."""
        return self._command

    def getConfig(self):
        """Get the list of configuration file strings."""
        return self._config

    def setConfig(self, config):
        """Set the list of configuration file strings."""

        # Check that the passed configuration is a list of strings.
        if _is_list_of_strings(config):
            self._config = config
            self.writeConfig(self._config_file)

        # The user has passed a path to a file.
        elif path.isfile(config):

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
        elif path.isfile(config):

            # Read the contents of the file.
            with open(file, "r") as f:
                for line in f:
                    self._config.append(line)

            # Write the new configuration file.
            self.writeConfig(self._config_file)

        else:
            raise ValueError("'config' must be a string, list of strings, or a file path.")

    def writeConfig(self, file):
        """Write the configuration to file."""
        with open(file, "w") as f:
            for line in self._config:
                f.write("%s\n" % line)

    def getArgs(self):
        """Get the dictionary of command-line arguments."""
        return self._args

    def getArgString(self):
        """Get the command-line arguments string."""
        return ' '.join(self.getArgStringList())

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
        if isinstance(args, OrderedDict):
            self._args = args

        elif isinstance(args, dict):
            self._args = OrderedDict(args)

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
        if isinstance(args, dict) or isinstance(args, OrderedDict):
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
            self._runtime = (timer() - self._timer) / 60
            return self._runtime

        # The process has finished.
        else:
            # Return the runtime and reset the timer.
            if self._timer is not None:
                self._runtime = (timer() - self._timer) / 60
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
    return AABox(coord)

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
    amino_acids = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE',
        'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'SEC',
        'VAL', 'TRP', 'TYR']

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
            backbone = m.atoms(AtomName("CA", CaseInsensitive) *
                               AtomName("N",  CaseInsensitive) *
                               AtomName("C",  CaseInsensitive) *
                               AtomName("O",  CaseInsensitive))

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
