"""
@package biosimspace
@author  Lester Hedges
@brief   A class for running simulations using AMBER.
"""

from Sire import try_import
from Sire.Base import findExe, Process
from Sire.IO import AmberPrm, AmberRst7, MoleculeParser

from . import process
from ..Protocol.protocol import Protocol, ProtocolType
from ..Trajectory.trajectory import Trajectory

from math import ceil, floor
from os import chdir, environ, getcwd, path
from re import findall
from shutil import copyfile
from time import sleep
from timeit import default_timer as timer
from warnings import warn

try:
    pygtail = try_import("pygtail")
except ImportError:
    raise ImportError("Pygtail is not installed. Please install pygtail in order to use BioSimSpace.")

try:
    watchdog = try_import("watchdog")
    from watchdog.events import PatternMatchingEventHandler
    from watchdog.observers import Observer
except ImportError:
    raise ImportError("Watchdog is not installed. Please install watchdog in order to use BioSimSpace.")

class Watcher:
    """A class to watch for changes to the AMBER energy info file. An event handler
       is used trigger updates to the energy dictionary each time the file is modified.
    """
    def __init__(self, proc):
        """Constructor.

           Positional arguments:

           proc -- The Amber Process object.
        """

        self._process = proc
        self._observer = Observer()

    def start(self):
        """Start the file watcher."""

        # Setup the event handler and observer.
        event_handler = Handler(self._process)
        self._observer.schedule(event_handler, self._process._work_dir)
        self._observer.daemon = True
        self._observer.start()

class Handler(PatternMatchingEventHandler):
    """An event handler to trigger updates to the energy dictionary each time
       the log file is changed.
    """

    # Overwrite defaults.
    case_sensitive = False
    ignore_directories = True
    ignore_patterns = []
    patterns = "*.nrg"

    def __init__(self, proc):
        """Constructor.

           Positional arguments:

           proc -- The Amber Process object.
        """
        self._process = proc

    def on_any_event(self, event):
        """Update the dictionary when the file is modified.

           Positional arguments:

           event -- The file system event.
        """

        # N.B.
        #
        # Since the watchdog package is cross-platform it doesn't support
        # detection of "close-write" operations, so multiple "modified" events
        # can be triggered while the log file is being written. As such, we
        # check to see if the file has been updated by seeing whether the
        # NSTEP record is different to the most recent entry in the dictionary.
        # So far, no issues have been found with processing partially written
        # files, i.e. duplicate or missing records.

        if event.event_type == 'modified':
            # If this is the first time the file has been modified since the
            # process started, then wipe the dictionary and flag that the file
            # is now being watched.
            if not self._process._is_watching:
                self._process._stdout_dict = process.MultiDict()
                self._process._is_watching = True

            # Now update the dictionary with any new records.
            self._process._update_energy_dict()

class Amber(process.Process):
    """A class for running simulations using AMBER."""

    def __init__(self, system, protocol, exe=None, name="amber",
            work_dir=None, seed=None):
        """Constructor.

           Positional arguments:

           system   -- The molecular system.
           protocol -- The protocol for the AMBER process.

           Keyword arguments:

           exe      -- The full path to the AMBER executable.
           name     -- The name of the process.
           work_dir -- The working directory for the process.
           seed     -- A random number seed.
        """

        # Call the base class constructor.
        super().__init__(system, protocol, name, work_dir, seed)

        # This process can generate trajectory data.
        self._has_trajectory = True

        # If the path to the executable wasn't specified, then search
        # for it in $PATH. For now, we'll just search for 'sander', which
        # is available free as part of AmberTools. In future, we will
        # look for all possible executables in order of preference: pmemd.cuda,
        # pmemd, sander, etc., as well as their variants, e.g. pmemd.MPI.
        if exe is None:
            # Search AMBERHOME, if set.
            if "AMBERHOME" in environ:
                amber_home = environ.get("AMBERHOME")
                if path.isfile("%s/bin/sander" % amber_home):
                    self._exe = "%s/bin/sander" % amber_home

            # Search PATH.
            else:
                self._exe = findExe("sander").absoluteFilePath()

        else:
            # Make sure executable exists.
            if path.isfile(exe):
                self._exe = protocol
            else:
                raise IOError(('AMBER executable doesn\'t exist: "{x}"').format(x=exe))

        # Initialise the energy dictionary and header.
        self._stdout_dict = process.MultiDict()

        # Create the name of the energy output file and wipe the
        # contents of any existing file.
        self._nrg_file = "%s/%s.nrg" % (self._work_dir, name)
        open(self._nrg_file, 'w').close()

        # Initialise the energy watcher.
        self._watcher = None
        self._is_watching = False

        # The names of the input files.
        self._rst_file = "%s/%s.rst7" % (self._work_dir, name)
        self._top_file = "%s/%s.prm7" % (self._work_dir, name)

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
        self._input_files = [self._config_file, self._rst_file, self._top_file]

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
        prm.writeToFile(self._top_file)

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
            if self._protocol.gas_phase:
                has_box = False
            else:
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
            if not has_box:
                self.addToConfig("  ntb=0,")                # No periodic box.
                self.addToConfig("  cut=999.,")             # Non-bonded cut-off.
            else:
                self.addToConfig("  cut=8.0,")              # Non-bonded cut-off.
            self.addToConfig(" /")

        # Add configuration variables for an equilibration simulation.
        elif self._protocol.type() == ProtocolType.EQUILIBRATION:

            # Convert the timestep to nanoseconds.
            timestep = self._protocol.timestep * 1e-6

            # Work out the number of integration steps.
            steps = ceil(self._protocol.runtime / timestep)

            # Set the random number seed.
            if self._is_seeded:
                seed = self._seed
            else:
                seed = -1

            # Convert the timestep to picoseconds.
            timestep = self._protocol.timestep * 1e-3

            self.addToConfig("Equilibration.")
            self.addToConfig(" &cntrl")
            self.addToConfig("  ig=%d," % seed)             # Random number seed.
            self.addToConfig("  ntx=1,")                    # Only read coordinates from file.
            self.addToConfig("  ntxo=1,")                   # Output coordinates in ASCII.
            self.addToConfig("  ntpr=100,")                 # Output energies every 100 steps.
            self.addToConfig("  ntwr=500,")                 # Save restart configuration every 500 steps.
            self.addToConfig("  irest=0,")                  # Don't restart.
            self.addToConfig("  dt=%.3f," % timestep)       # Time step.
            self.addToConfig("  nstlim=%d," % steps)        # Number of integration steps.
            self.addToConfig("  ntc=2,")                    # Enable SHAKE.
            self.addToConfig("  ntf=2,")                    # Don't calculate forces for constrained bonds.
            self.addToConfig("  ntt=3,")                    # Langevin dynamics.
            self.addToConfig("  gamma_ln=2,")               # Collision frequency (ps).
            if not has_box:
                self.addToConfig("  ntb=0,")                # No periodic box.
                self.addToConfig("  cut=999.,")             # Non-bonded cut-off.
            else:
                self.addToConfig("  cut=8.0,")              # Non-bonded cut-off.
            self.addToConfig("  ntp=1,")                    # Isotropic pressure scaling.
            self.addToConfig("  pres0=1.01325,")            # Atompspheric pressure.

            # Restrain the backbone.
            if self._protocol.is_restrained:
                self.addToConfig("  ntr=1,")
                self.addToConfig("  restraint_wt=2,")
                self.addToConfig("  restraintmask='@CA,C,O,N',")

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

            # Convert the timestep to nanoseconds.
            timestep = self._protocol.timestep * 1e-6

            # Work out the number of integration steps.
            steps = ceil(self._protocol.runtime / timestep)

            # Set the random number seed.
            if self._seed is None:
                seed = -1
            else:
                seed = self._seed

            # Convert the timestep to picoseconds.
            timestep = self._protocol.timestep * 1e-3

            self.addToConfig("Production.")
            self.addToConfig(" &cntrl")
            self.addToConfig("  ig=%d," % seed)             # Random number seed.
            if self._protocol.restart:
                self.addToConfig("  ntx=5,")                # Read coordinates and velocities.
            else:
                self.addToConfig("  ntx=1,")                # Only read coordinates.
            self.addToConfig("  ntxo=1,")                   # Output coordinates in ASCII.
            self.addToConfig("  ntpr=100,")                 # Output energies every 100 steps.
            self.addToConfig("  ntwr=500,")                 # Save restart configuration every 500 steps.
            self.addToConfig("  ntwx=%d,"                   # Trajectory sampling frequency.
                % floor(steps / self._protocol.frames))
            if self._protocol.restart:
                self.addToConfig("  irest=1,")              # Restart using previous velocities.
            else:
                self.addToConfig("  irest=0,")              # Don't restart.
            self.addToConfig("  dt=%.3f," % timestep)       # Time step.
            self.addToConfig("  nstlim=%d," % steps)        # Number of integration steps.
            self.addToConfig("  ntc=2,")                    # Enable SHAKE.
            self.addToConfig("  ntf=2,")                    # Don't calculate forces for constrained bonds.
            self.addToConfig("  ntt=3,")                    # Langevin dynamics.
            self.addToConfig("  gamma_ln=2,")               # Collision frequency (ps).
            if not has_box:
                self.addToConfig("  ntb=0,")                # No periodic box.
                self.addToConfig("  cut=999.,")             # Non-bonded cut-off.
            else:
                self.addToConfig("  cut=8.0,")              # Non-bonded cut-off.
            if not self._protocol.restart:
                self.addToConfig("  tempi=%.2f,"            # Initial temperature.
                    % self._protocol.temperature)
            self.addToConfig("  temp0=%.2f,"                # Target temperature.
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
        self.setArg("-r", "%s.crd" % self._name)            # Restart file.
        self.setArg("-inf", "%s.nrg" % self._name)          # Energy info file.

        # Skip if the user has passed a custom config.
        if not self._is_custom:

            # Append a reference file if this a constrained equilibration.
            if self._protocol.type() == ProtocolType.EQUILIBRATION:
                if self._protocol.is_restrained:
                    self.setArg("-ref", "%s.rst7" % self._name)

            # Append a trajectory file if this is a production run.
            elif self._protocol.type() == ProtocolType.PRODUCTION:
                self.setArg("-x", "%s.nc" % self._name)

    def start(self):
        """Start the AMBER simulation. """

        # Process is already running.
        if self._process is not None:
            if self._process.isRunning():
                return

        # Reset the watcher.
        self._is_watching = False

        # Store the current working directory.
        dir = getcwd()

        # Change to the working directory for the process.
        # This avoid problems with relative paths.
        chdir(self._work_dir)

        # Create the arguments string list.
        args = self.getArgStringList()

        # Write the command-line process to a README.txt file.
        with open("README.txt", "w") as f:

            # Set the command-line string.
            self._command = "%s " % self._exe + self.getArgString()

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

	# Watch the energy info file for changes.
        self._watcher = Watcher(self)
        self._watcher.start()

    def getSystem(self, block='AUTO'):
        """Get the latest molecular configuration as a Sire system.

           Keyword arguments:

           block -- Whether to block until the process has finished running.
        """

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block is 'AUTO' and self._is_blocked:
            self.wait()

        # Create the name of the restart CRD file.
        restart = "%s/%s.crd" % (self._work_dir, self._name)

        # Check that the file exists.
        if path.isfile(restart):
            # Create and return the molecular system.
            return MoleculeParser.read(restart, self._top_file)

        else:
            return None

    def getCurrentSystem(self):
        """Get the latest molecular configuration as a Sire system."""
        return self.getSystem(block=False)

    def getTrajectory(self, block='AUTO'):
        """Return a trajectory object."""

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block is 'AUTO' and self._is_blocked:
            self.wait()

        return Trajectory(process=self)

    def getRecord(self, record, time_series=False, block='AUTO'):
        """Get a record from the stdout dictionary.

           Positional arguments:

           record      -- The record key.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block is 'AUTO' and self._is_blocked:
            self.wait()

        return self._get_stdout_record(record.strip().upper(), time_series)

    def getCurrentRecord(self, record, time_series=False):
        """Get a current record from the stdout dictionary.

           Positional arguments:

           record      -- The record key.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self._get_stdout_record(record.strip().upper(), time_series)

    def getRecords(self, block='AUTO'):
        """Return the dictionary of stdout time-series records.

           Keyword arguments:

           block -- Whether to block until the process has finished running.
        """

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block is 'AUTO' and self._is_blocked:
            self.wait()

        return self._stdout_dict

    def getCurrentRecords(self):
        """Return the current dictionary of stdout time-series records."""
        return getRecords(block=False)

    def getTime(self, time_series=False, block='AUTO'):
        """Get the time (in nanoseconds).

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """

        # No time records for minimisation protocols.
        if self._protocol.type() == ProtocolType.MINIMISATION:
            return None

        # Get the list of time steps.
        time_steps = self.getRecord('TIME(PS)', time_series, block)

        # Convert from picoseconds to nanoseconds.
        if time_steps is not None:
            if time_series:
                return [x * 1e-3 for x in time_steps]
            else:
                return 1e-3 * time_steps

    def getCurrentTime(self, time_series=False):
        """Get the current time (in nanoseconds).

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getTime(time_series, block=False)

    def getStep(self, time_series=False, block='AUTO'):
        """Get the number of integration steps.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord('NSTEP', time_series, block)

    def getCurrentStep(self, time_series=False):
        """Get the current number of integration steps.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getStep(time_series, block=False)

    def getBondEnergy(self, time_series=False, block='AUTO'):
        """Get the bond energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord('BOND', time_series, block)

    def getCurrentBondEnergy(self, time_series=False):
        """Get the current bond energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getBondEnergy(time_series, block=False)

    def getAngleEnergy(self, time_series=False, block='AUTO'):
        """Get the angle energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord('ANGLE', time_series, block)

    def getCurrentAngleEnergy(self, time_series=False):
        """Get the current angle energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getAngleEnergy(time_series, block=False)

    def getDihedralEnergy(self, time_series=False, block='AUTO'):
        """Get the dihedral energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord('DIHED', time_series, block)

    def getCurrentDihedralEnergy(self, time_series=False):
        """Get the current dihedral energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getDihedralEnergy(time_series, block=False)

    def getElectrostaticEnergy(self, time_series=False, block='AUTO'):
        """Get the electrostatic energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord('EELECT', time_series, block)

    def getCurrentElectrostaticEnergy(self, time_series=False):
        """Get the current dihedral energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getElectrostaticEnergy(time_series, block=False)

    def getElectrostaticEnergy14(self, time_series=False, block='AUTO'):
        """Get the electrostatic energy between atoms 1 and 4.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord('1-4 EEL', time_series, block)

    def getCurrentElectrostaticEnergy14(self, time_series=False):
        """Get the current electrostatic energy between atoms 1 and 4.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getElectrostaticEnergy14(time_series, block=False)

    def getVanDerWaalsEnergy(self, time_series=False, block='AUTO'):
        """Get the Van der Vaals energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord('VDWAALS', time_series, block)

    def getCurrentVanDerWaalsEnergy(self, time_series=False):
        """Get the current Van der Vaals energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getVanDerWaalsEnergy(time_series, block=False)

    def getHydrogenBondEnergy(self, time_series=False, block='AUTO'):
        """Get the hydrogen bond energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord('EHBOND', time_series, block)

    def getCurrentHydrogenBondEnergy(self, time_series=False):
        """Get the current hydrogen bond energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getHydrogenBondEnergy(time_series, block=False)

    def getRestraintEnergy(self, time_series=False, block='AUTO'):
        """Get the restraint energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord('RESTRAINT', time_series, block)

    def getCurrentRestraintEnergy(self, time_series=False):
        """Get the current restraint energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getRestraintEnergy(time_series, block=False)

    def getPotentialEnergy(self, time_series=False, block='AUTO'):
        """Get the potential energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord('EPTOT', time_series, block)

    def getCurrentPotentialEnergy(self, time_series=False):
        """Get the current potential energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getPotentialEnergy(time_series, block=False)

    def getKineticEnergy(self, time_series=False, block='AUTO'):
        """Get the kinetic energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord('EKTOT', time_series, block)

    def getCurrentKineticEnergy(self, time_series=False):
        """Get the current kinetic energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getKineticEnergy(time_series, block=False)

    def getNonBondedEnergy14(self, time_series=False, block='AUTO'):
        """Get the non-bonded energy between atoms 1 and 4.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord('1-4 NB', time_series, block)

    def getCurrentNonBondedEnergy14(self, time_series=False):
        """Get the current non-bonded energy between atoms 1 and 4.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getNonBondedEnergy14(time_series, block=False)

    def getTotalEnergy(self, time_series=False, block='AUTO'):
        """Get the total energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        if self._protocol.type() == ProtocolType.MINIMISATION:
            return self.getRecord('ENERGY', time_series, block)
        else:
            return self.getRecord('ETOT', time_series, block)

    def getCurrentTotalEnergy(self, time_series=False):
        """Get the current total energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getTotalEnergy(time_series, block=False)

    def getCentreOfMassKineticEnergy(self, time_series=False, block='AUTO'):
        """Get the kinetic energy of the centre of mass in translation.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord('EKCMT', time_series, block)

    def getCurrentCentreOfMassKineticEnergy(self, time_series=False):
        """Get the current kinetic energy of the centre of mass in translation.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getCentreOfMassKineticEnergy(time_series, block=False)

    def getVirial(self, time_series=False, block='AUTO'):
        """Get the virial.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord('VIRIAL', time_series, block)

    def getCurrentVirial(self, time_series=False):
        """Get the current virial.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getVirial(time_series, block=False)

    def getTemperature(self, time_series=False, block='AUTO'):
        """Get the temperature.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord('TEMP(K)', time_series, block)

    def getCurrentTemperature(self, time_series=False):
        """Get the current temperature.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getTemperature(time_series, block=False)

    def getPressure(self, time_series=False, block='AUTO'):
        """Get the temperature.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord('PRESS', time_series, block)

    def getCurrentPressure(self, time_series=False):
        """Get the current pressure.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getPressure(time_series, block=False)

    def getVolume(self, time_series=False, block='AUTO'):
        """Get the volume.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord('VOLUME', time_series, block)

    def getCurrentVolume(self, time_series=False):
        """Get the current volume.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getVolume(time_series, block=False)

    def getDensity(self, time_series=False):
        """Get the density.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord('DENSITY', time_series, block)

    def getCurrentDensity(self, time_series=False):
        """Get the current density.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getDensity(time_series, block=False)

    def _update_energy_dict(self):
        """Read the energy info file and update the dictionary."""

        # Flag that this isn't a header line.
        is_header = False

        # Open the file for reading.
        with open(self._nrg_file, "r") as f:

            # Loop over all of the lines.
            for line in f:

                # Skip empty lines and summary reports.
                if len(line) > 0 and line[0] is not '|':

                    # The output format is different for minimisation protocols.
                    if self._protocol.type() == ProtocolType.MINIMISATION:

                        # No equals sign in the line.
                        if "=" not in line:

                            # Split the line using whitespace.
                            data = line.upper().split()

                            # If we find a header, jump to the top of the loop.
                            if len(data) > 0:
                                if data[0] == 'NSTEP':
                                    is_header = True
                                    continue

                        # Process the header record.
                        if is_header:

                            # Split the line using whitespace.
                            data = line.upper().split()

                            # The file hasn't been updated.
                            if 'NSTEP' in self._stdout_dict and data[0] == self._stdout_dict['NSTEP'][-1]:
                                return

                            else:
                                # Add the timestep and energy records to the dictionary.
                                self._stdout_dict['NSTEP'] = data[0]
                                self._stdout_dict['ENERGY'] = data[1]

                                # Turn off the header flag now that the data has been recorded.
                                is_header = False

                    # All other protocols have output that is formatted as RECORD = VALUE.

                    # Use a regex search to split the line into record names and values.
                    records = findall("(\d*\-*\d*\s*[A-Z]+\(*[A-Z]*\)*)\s*=\s*(\-*\d+\.?\d*)", line.upper())

                    # Append each record to the dictionary.
                    for key, value in records:

                        # Strip whitespace from the record key.
                        key = key.strip()

                        # The file hasn't been updated.
                        if key == 'NSTEP':
                            if 'NSTEP' in self._stdout_dict and value == self._stdout_dict['NSTEP'][-1]:
                                return
                            else:
                                self._stdout_dict[key] = value
                        else:
                            self._stdout_dict[key] = value

    def kill(self):
        """Kill the running process."""

        # Stop and join the watchdog observer.
        if self._watcher is not None:
            self._watcher._observer.stop()
            self._watcher._observer.join()

        # Kill the process.
        self._process.kill()

    def wait(self):
        """Wait for the process to finish."""

        # Loop until the process is finished.
        # For some reason we can't use Sire.Base.Process.wait() since it
        # doesn't work properly with the background threads used for the
        # watchdog observer.
        while self._process.isRunning():
            sleep(1)

        # Stop and join the watchdog observer.
        self._watcher._observer.stop()
        self._watcher._observer.join()

    def _get_stdout_record(self, key, time_series=False):
        """Helper function to get a stdout record from the dictionary.

           Positional arguments:

           key        -- The record key.

           Keyword arguments:

           time_series -- Whether to return a time series of records.
        """

        # No data!
        if len(self._stdout_dict) is 0:
            return None

        if type(time_series) is not bool:
            warn("Non-boolean time-series flag. Defaulting to False!")
            time_series = False

        # Return the list of dictionary values.
        if time_series:
            try:
                if key is 'NSTEP':
                    return [int(x) for x in self._stdout_dict[key]]
                else:
                    return [float(x) for x in self._stdout_dict[key]]

            except KeyError:
                return None

        # Return the most recent dictionary value.
        else:
            try:
                if key is 'NSTEP':
                    return int(self._stdout_dict[key][-1])
                else:
                    return float(self._stdout_dict[key][-1])

            except KeyError:
                return None

    def _get_trajectory_files(self):
        """Get all files associated with the molecular trajectory."""

        # Name of the trajectory file.
        traj_file = "%s/%s.nc" % (self._work_dir, self._name)

        # Return the trajectory and topology file.
        if path.isfile("%s/%s.nc" % (self._work_dir, self._name)):
            return (traj_file, self._top_file)

        # No trajectory file.
        else:
            return None
