"""
@package biosimspace
@author  Lester Hedges
@brief   A class for running simulations using AMBER.
"""

from Sire.Base import findExe, Process
from Sire.IO import AmberPrm, AmberRst7, MoleculeParser

from . import process
from ..Protocol.protocol import Protocol, ProtocolType

from math import ceil, floor
from os import chdir, getcwd, path
from timeit import default_timer as timer
from warnings import warn

import __main__ as main
import os

try:
    from Sire import try_import
    pygtail = try_import("pygtail")
except ImportError:
    raise ImportError("Pygtail is not installed. Please install pygtail in order to use BioSimSpace.")

class AmberProcess(process.Process):
    """A class for running simulations using AMBER."""

    def __init__(self, system, protocol, exe=None, name="amber",
            work_dir=None, seed=None):
        """Constructor.

           Keyword arguments:

           system        -- The molecular system.
           protocol      -- The protocol for the AMBER process.
           exe           -- The full path to the AMBER executable.
           name          -- The name of the process.
           work_dir      -- The working directory for the process.
           seed          -- A random number seed.
        """

        # Call the base class constructor.
        super().__init__(system, protocol, name, work_dir)

        # If the path to the executable wasn't specified, then search
        # for it in $PATH. For now, we'll just search for 'sander', which
        # is available free as part of AmberTools. In future, we will
        # look for all possible executables in order of preference: pmemd.cuda,
        # pmemd, sander, etc., as well as their variants, e.g. pmemd.MPI.
        if exe is None:
            self._exe = findExe("sander").absoluteFilePath()

        else:
            # Make sure executable exists.
            if path.isfile(exe):
                self._exe = protocol
            else:
                raise IOError(('AMBER executable doesn\'t exist: "{x}"').format(x=exe))

        # Initialise the energy dictionary and header.
        self._stdout_dict = process._MDict()
        self._stdout_title = None

        # The names of the input files.
        self._crd_file = "%s/%s.crd" % (self._work_dir, name)
        self._top_file = "%s/%s.top" % (self._work_dir, name)

        # Set the path for the AMBER configuration file.
        # The 'protocol' argument may contain the path to a custom file.

        # Set the config file name.
        if not self._is_custom:
            self._config_file = "%s/%s.amber" % (self._work_dir, name)

        # The user has supplied a custom config file.
        else:
            # Make sure the file exists.
            if path.isfile(protocol):
                self._config_file = protocol
            else:
                raise IOError(('AMBER configuration file doesn\'t exist: "{x}"').format(x=config_file))

        # Create the list of input files.
        self._input_files = [self._config_file, self._crd_file, self._top_file]

        # Now set up the working directory for the process.
        self._setup()

    def _setup(self):
        """Setup the input files and working directory ready for simulation."""

        # Create the input files...

        # CRD file (coordinates).
        crd = AmberRst7(self._system)
        crd.writeToFile(self._crd_file)

        # TOP file (topology).
        top = AmberPrm(self._system)
        top.writeToFile(self._top_file)

        # Generate the AMBER configuration file.
        # Skip if the user has passed a custom config.
        if not self._is_custom:
            self._generate_config_file()

        # Return the list of input files.
        return self._input_files

    def _generate_config_file(self):
        """Generate an AMBER configuration file."""

        # Check whether the system contains periodic box information.
        # For now, well not attempt to generate a box if the system property
        # is missing. If no box is present, we'll assume a non-periodic simulation.
        if 'space' in self._system.propertyKeys():
            has_box = True
        else:
            has_box = False

        # Open the configuration file for writing.
        f = open(self._config_file, "w")

        # Add configuration variables for a minimisation simulation.
        if self._protocol.type() == ProtocolType.MINIMISATION:
            f.write("Minimisation.\n")
            f.write(" &cntrl\n")
            f.write("  imin=1,\n")                               # Minisation simulation.
            f.write("  ntx=1,\n")                                # Only read coordinates from file.
            f.write("  ntxo=1,\n")                               # Output coordinates in ASCII.
            f.write("  ntpr=100,\n")                             # Output energies every 100 steps.
            f.write("  irest=0,\n")                              # Don't restart.
            f.write("  maxcyc=%s,\n" % self._protocol.steps)     # Set the number of steps.
            f.write("  cut=8.0,\n")                              # Non-bonded cut-off.
            f.write(" /\n")

        # Add configuration variables for an equilibration simulation.
        elif self._protocol.type() == ProtocolType.EQUILIBRATION:
            pass

        # Add configuration variables for a production simulation.
        elif self._protocol.type() == ProtocolType.PRODUCTION:
            pass

        # Close the configuration file.
        f.close()

    def start(self):
        """Start the AMBER simulation."""

        # Store the current working directory.
        dir = getcwd()

        # Change to the working directory for the process.
        # This avoid problems with relative paths.
        chdir(self._work_dir)

        # Create a list of the command-line arguments.
        args = ["-O",                                # Overwrite.
                "-i", "%s.amber" % self._name,       # Input file.
                "-p", "%s.top" % self._name,         # Topology file.
                "-c", "%s.crd" % self._name,         # Coordinate file.
                "-o", "stdout",                      # Redirect to stdout.
                "-r", "%s.restart.crd" % self._name, # Restart file.
                "-inf", "%s.nrg" % self._name]       # Energy info file.

        # Append a trajectory file if this is a production run.
        if self._protocol.type() == ProtocolType.PRODUCTION:
            args.append("-x %s.trajectory.crd" % self._name)

        # Write the command-line process to a README.txt file.
        with open("README.txt", "w") as f:

            # Set the command-line string.
            self._command = "%s " % self._exe + ' '.join(args)

            # Write the command to file.
            f.write("# AmberProcess was run with the following command:\n")
            f.write("%s\n" % self._command)

        # Start the timer.
        self._timer = timer()

        # Start the simulation.
        self._process = Process.run(self._exe, args,
            "%s.out"  % self._name, "%s.err"  % self._name)

        # Change back to the original working directory.
        chdir(dir)

    def getSystem(self):
        """Get the latest molecular configuration as a Sire system."""

        # Create the name of the restart CRD file.
        restart = "%s/%s.restart.crd" % (self._work_dir, self._name)

        # Check that the file exists.
        if path.isfile(restart):
            # Create and return the molecular system.
            return MoleculeParser.read(restart, self._top_file)

        else:
            return None
