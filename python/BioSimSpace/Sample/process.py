"""
@package biosimspace
@author  Lester Hedges
@brief   Base class and helper functions for the various sample modules.
"""

import Sire.Mol

from ..Protocol.protocol import Protocol

from operator import add, sub
from os import makedirs, path, remove
from timeit import default_timer as timer

import tempfile

try:
    pygtail = Sire.try_import("pygtail")
except ImportError:
    raise ImportError("Pygtail is not installed. Please install pygtail in order to use BioSimSpace.")

class Process():
    """Base class for running different biomolecular simulation processes."""

    def __init__(self, system, protocol, name="process", work_dir=None, seed=None):
        """Constructor.

           Keyword arguments:

           system   -- The molecular system.
           protocol -- The protocol for the process.
           name     -- The name of the process.
           work_dir -- The working directory for the process.
           seed     -- A random number seed.
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
        self._seed = seed

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

            # Create the directory if it doesn't already exist..
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

def _compute_box_size(system, tol=0.05, buffer=0):
    """Compute the box size and origin from the atomic coordinates.

       Keyword arguments:

       system -- A Sire molecular system.
       tol    -- The tolerance for determining whether the box is square
                 and whether the origin lies at (0, 0, 0).
       buffer -- The percentage by which to expand the box to account for
                 periodic wrapping.
    """

    # Store the list of molecule indices.
    mol_nums = system.molNums()

    # Initialise the min and max box size for each dimension.
    box_min = [1000000]  * 3
    box_max = [-1000000] * 3

    # The number of water molecules.
    # We'll assume that all 3- or 4-atom molecules are water.
    num_water = 0

    # Loop over all of the molecules.
    for num in mol_nums:

        # Store the molecule.
        mol = system[num]

        # Get the number of atoms.
        num_atoms = mol.nAtoms()

        # Is this a water molecule?
        if num_atoms is 3 or num_atoms is 4:
            num_water += 1

        # Loop over all atoms in the molecule.
        for atom in mol.atoms():

            # Extract the atomic coordinates.
            try:
                coord = atom.property("coordinates")

            except UserWarning:
               raise

            # Check coordinates against the current min/max.
            for x in range(0, 3):

               if coord[x] < box_min[x]:
                   box_min[x] = coord[x]

               elif coord[x] > box_max[x]:
                   box_max[x] = coord[x]

    # Calculate the base length of the simulation box.
    box_size = list(map(sub, box_max, box_min))

    # Calculate the centre of the box.
    box_origin = [x * 0.5 for x in list(map(add, box_min, box_max))]

    # Store the base length with the maximum size and its index.
    max_size = 0
    for index, size in enumerate(box_size):
        if size > max_size:
            max_size, index = size, index

    # Loop over all box dimensions.
    for x in range(0, 3):

        # Assume the box is square if the base lengths are similar.
        if box_size[x] > (1 - tol) * max_size:
            box_size[x] = max_size
            box_origin[x] = box_origin[index]

        # Assume the origin is at zero if the centre of mass is
        # close to (0, 0, 0)
        if box_origin[x] / box_size[x] < tol:
            box_origin[x] = 0

    # Add a buffer to the box size to account for atom wrapping.
    box_size = [x * (1 + buffer) for x in box_size]

    # Return the box size, origin, and whether the system is solvated.
    return (tuple(box_size), tuple(box_origin), num_water > 0)
