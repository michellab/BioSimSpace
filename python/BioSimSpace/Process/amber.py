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

class Amber(process.Process):
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
        super().__init__(system, protocol, name, work_dir, seed)

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
        self._rst_file = "%s/%s.rst7" % (self._work_dir, name)
        self._prm_file = "%s/%s.prm7" % (self._work_dir, name)

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
        self._input_files = [self._config_file, self._rst_file, self._prm_file]

        # Now set up the working directory for the process.
        self._setup()

    def _setup(self):
        """Setup the input files and working directory ready for simulation."""

        # Create the input files...

        # RST file (coordinates).
        rst = AmberRst7(self._system)
        rst.writeToFile(self._rst_file)

        # PRM file (topology).
        prm = AmberPrm(self._system)
        prm.writeToFile(self._prm_file)

        # Generate the AMBER configuration file.
        # Skip if the user has passed a custom config.
        if not self._is_custom:
            self._generate_config()
            self.writeConfig(self._config_file)

        # Generate the dictionary of command-line arguments.
        self._generate_args()

        # Return the list of input files.
        return self._input_files

    def _generate_config(self):
        """Generate AMBER configuration file strings."""

        # Clear the existing configuration list.
        self._config = []

        # Check whether the system contains periodic box information.
        # For now, well not attempt to generate a box if the system property
        # is missing. If no box is present, we'll assume a non-periodic simulation.
        if 'space' in self._system.propertyKeys():
            has_box = True
        else:
            has_box = False

        # Add configuration variables for a minimisation simulation.
        if self._protocol.type() == ProtocolType.MINIMISATION:
            self.addToConfig("Minimisation")
            self.addToConfig(" &cntrl")
            self.addToConfig("  imin=1,")                   # Minisation simulation.
            self.addToConfig("  ntx=1,")                    # Only read coordinates from file.
            self.addToConfig("  ntxo=1,")                   # Output coordinates in ASCII.
            self.addToConfig("  ntpr=100,")                 # Output energies every 100 steps.
            self.addToConfig("  irest=0,")                  # Don't restart.
            self.addToConfig("  maxcyc=%d,"
                % self._protocol.steps)                     # Set the number of steps.
            self.addToConfig("  cut=8.0,")                  # Non-bonded cut-off.
            self.addToConfig(" /")

        # Add configuration variables for an equilibration simulation.
        elif self._protocol.type() == ProtocolType.EQUILIBRATION:
            # Work out the number of integration steps.
            steps = ceil(self._protocol.runtime / 2e-6)

            # Set the random number seed.
            if self._seed is None:
                seed = -1
            else:
                seed = self._seed

            self.addToConfig("Equilibration.")
            self.addToConfig(" &cntrl")
            self.addToConfig("  ig=%d," % seed)             # Random number seed.
            self.addToConfig("  ntx=1,")                    # Only read coordinates from file.
            self.addToConfig("  ntxo=1,")                   # Output coordinates in ASCII.
            self.addToConfig("  ntpr=100,")                 # Output energies every 100 steps.
            self.addToConfig("  irest=0,")                  # Don't restart.
            self.addToConfig("  dt=0.002,")                 # Time step (2fs).
            self.addToConfig("  nstlim=%d," % steps)        # Number of integration steps.
            self.addToConfig("  ntc=2,")                    # Enable SHAKE.
            self.addToConfig("  ntf=2,")                    # Don't calculate forces for constrained bonds.
            self.addToConfig("  ntt=3,")                    # Langevin dynamics.
            self.addToConfig("  gamma_ln=2,")               # Collision frequency (ps).
            self.addToConfig("  cut=8.0,")                  # Non-bonded cut-off.
            self.addToConfig("  ntp=1,")                    # Isotropic pressure scaling.
            self.addToConfig("  pres0=1.01325,")            # Atompspheric pressure.

            # Restrain the backbone.
            if self._protocol.is_restrained:
                self.addToConfig("  ntr=1,")
                self.addToConfig("  restraint_wt = 2,")
                self.addToConfig("  restraintmask = '@CA,C,O,N',")

            # Heating/cooling protocol.
            if not self._protocol.isConstantTemp():
                self.addToConfig("  tempi=%.2f," % self._protocol.temperature_start)
                self.addToConfig("  temp0=%.2f," % self._protocol.temperature_end)
                self.addToConfig("  nmropt=1,")
                self.addToConfig(" /")
                self.addToConfig("&wt TYPE='TEMP0', istep1=0, istep2=%d, value1=%.2f, value2=%.2f /"
                    % (steps, self._protocol.temperature_start, self._protocol.temperature_end))
                self.addToConfig("&wt TYPE='END' /")

            # Constant temperature equilibration.
            else:
                self.addToConfig("  temp0=%.2f," % self._protocol.temperature_start)
                self.addToConfig(" /")

        # Add configuration variables for a production simulation.
        elif self._protocol.type() == ProtocolType.PRODUCTION:
            # Work out the number of integration steps.
            steps = ceil(self._protocol.runtime / 2e-6)

            # Set the random number seed.
            if self._seed is None:
                seed = -1
            else:
                seed = self._seed

            self.addToConfig("Production.")
            self.addToConfig(" &cntrl")
            self.addToConfig("  ig=%d," % seed)             # Random number seed.
            self.addToConfig("  ntx=1,")                    # Only read coordinates from file.
            self.addToConfig("  ntxo=1,")                   # Output coordinates in ASCII.
            self.addToConfig("  ntpr=100,")                 # Output energies every 100 steps.
            self.addToConfig("  ntwx=%d,"                   # Trajectory sampling frequency.
                % floor(steps / self._protocol.frames))
            self.addToConfig("  irest=0,")                  # Don't restart.
            self.addToConfig("  dt=0.002,")                 # Time step (2fs).
            self.addToConfig("  nstlim=%d," % steps)        # Number of integration steps.
            self.addToConfig("  ntc=2,")                    # Enable SHAKE.
            self.addToConfig("  ntf=2,")                    # Don't calculate forces for constrained bonds.
            self.addToConfig("  ntt=3,")                    # Langevin dynamics.
            self.addToConfig("  gamma_ln=2,")               # Collision frequency (ps).
            self.addToConfig("  cut=8.0,")                  # Non-bonded cut-off.
            self.addToConfig("  temp0=%.2f,"                # Temperature.
                % self._protocol.temperature)

            # Constant pressure control.
            if self._protocol.ensemble is 'NPT':
                self.addToConfig("  ntp=1,")                # Isotropic pressure scaling.
                self.addToConfig("  pres0=1.01325,")        # Atompspheric pressure.

            self.addToConfig(" /")

    def _generate_args(self):
        """Generate the dictionary of command-line arguments."""

        # Clear the existing arguments.
        self.clearArgs()

        # Add the default arguments.
        self.setArg("-O", True)                             # Overwrite.
        self.setArg("-i", "%s.amber" % self._name)          # Input file.
        self.setArg("-p", "%s.prm7" % self._name)           # Topology file.
        self.setArg("-c", "%s.rst7" % self._name)           # Coordinate file.
        self.setArg("-o", "stdout")                         # Redirect to stdout.
        self.setArg("-r", "%s.restart.crd" % self._name)    # Restart file.
        self.setArg("-inf", "%s.nrg" % self._name)          # Energy info file.

        # Skip if the user has passed a custom config.
        if not self._is_custom:

            # Append a reference file if this a constrained equilibration.
            if self._protocol.type() == ProtocolType.EQUILIBRATION:
                if self._protocol.is_restrained:
                    self.setArg("-ref", "%s.crd" % self._name)

            # Append a trajectory file if this is a production run.
            elif self._protocol.type() == ProtocolType.PRODUCTION:
                self.setArg("-x", "%s.trajectory.crd" % self._name)

    def start(self):
        """Start the AMBER simulation. """

        # Store the current working directory.
        dir = getcwd()

        # Change to the working directory for the process.
        # This avoid problems with relative paths.
        chdir(self._work_dir)

        # Create the arguments string list.
        args = self._generate_args_string()

        # Write the command-line process to a README.txt file.
        with open("README.txt", "w") as f:

            # Set the command-line string.
            self._command = "%s " % self._exe + ' '.join(args)

            # Write the command to file.
            f.write("# AMBER was run with the following command:\n")
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
            return MoleculeParser.read(restart, self._prm_file)

        else:
            return None
