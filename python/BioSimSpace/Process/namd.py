"""
@package biosimspace
@author  Lester Hedges
@brief   A class for running simulations using NAMD.
"""

from Sire.Base import findExe, Process
from Sire.IO import CharmmPSF, MoleculeParser, PDB2

from . import process
from ..Protocol.protocol import Protocol, ProtocolType

from math import ceil, floor
from os import chdir, getcwd, path
from timeit import default_timer as timer
from warnings import warn

import __main__ as main

try:
    from Sire import try_import
    pygtail = try_import("pygtail")
except ImportError:
    raise ImportError("Pygtail is not installed. Please install pygtail in order to use BioSimSpace.")

class Namd(process.Process):
    """A class for running simulations using NAMD."""

    def __init__(self, system, protocol, exe=None, name="namd",
            work_dir=None, charmm_params=True, seed=None):
        """Constructor.

           Keyword arguments:

           system        -- The molecular system.
           protocol      -- The protocol for the NAMD process.
           exe           -- The full path to the NAMD executable.
           name          -- The name of the process.
           work_dir      -- The working directory for the process.
           charmm_params -- Whether the parameter file is in CHARMM format.
           seed          -- A random number seed.
        """

        # Call the base class constructor.
        super().__init__(system, protocol, name, work_dir, seed)

        # If the path to the executable wasn't specified, then search
        # for it in $PATH.
        if exe is None:
            self._exe = findExe("namd2").absoluteFilePath()

        else:
            # Make sure executable exists.
            if path.isfile(exe):
                self._exe = protocol
            else:
                raise IOError(('NAMD executable doesn\'t exist: "{x}"').format(x=exe))

        # Set the parameter type.
        self._is_charmm_params = charmm_params

        # Initialise the stdout dictionary and title header.
        self._stdout_dict = process.MultiDict()
        self._stdout_title = None

        # The names of the input files.
        self._psf_file = "%s/%s.psf" % (self._work_dir, name)
        self._pdb_file = "%s/%s.pdb" % (self._work_dir, name)
        self._param_file = "%s/%s.params" % (self._work_dir, name)
        self._velocity_file = None
        self._restraint_file = None

        # Set the path for the NAMD configuration file.
        # The 'protocol' argument may contain the path to a custom file.

        # Set the config file name.
        if not self._is_custom:
            self._config_file = "%s/%s.namd" % (self._work_dir, name)

        # The user has supplied a custom config file.
        else:
            # Make sure the file exists.
            if path.isfile(protocol):
                self._config_file = protocol
            else:
                raise IOError(('NAMD configuration file doesn\'t exist: "{x}"').format(x=config_file))

        # Create the list of input files.
        self._input_files = [self._config_file, self._psf_file, self._pdb_file, self._param_file]

        # Now set up the working directory for the process.
        self._setup()

    @property
    def is_charmm_params(self):
        """Return whether the parameters are in CHARMM format."""
        return self._is_charmm_params

    @is_charmm_params.setter
    def is_charmm_params(self, charmm_params):
        """Set the parameter type."""

        if type(charmm_params) is bool:
            self._is_charmm_params = charmm_params

        else:
            warn("Non-boolean CHARMM parameter flag. Defaulting to True!")
            self._is_charmm_params = True

    def _setup(self):
        """Setup the input files and working directory ready for simulation."""

        # Create the input files...

        # PSF and parameter files.
        psf = CharmmPSF(self._system)
        psf.writeToFile(self._psf_file)

        # PDB file.
        pdb = PDB2(self._system)
        pdb.writeToFile(self._pdb_file)

        # Try to write a PDB "velocity" restart file.
        # The file will only be generated if all atoms in self._system have
        # a "velocity" property.

        # First generate a name for the velocity file.
        velocity_file = path.splitext(self._pdb_file)[0] + ".vel"

        # Write the velocity file.
        has_velocities = pdb.writeVelocityFile(velocity_file)

        # If a file was written, store the name of the file and update the
        # list of input files.
        if has_velocities:
            self._velocity_file = velocity_file
            self._input_files.append(self._velocity_file)

        # NAMD requires donor, acceptor, and non-bonded exclusion record entries
        # in the PSF file. We check that these are present and append blank
        # records if not. The records are actually redundant (they have no
        # affect on the MD) so could be stripped (or zeroed) by the CharmmPSF
        # parser.
        has_donors = False
        has_acceptors = False
        has_non_bonded = False

        # Open the PSF file for reading.
        with open(self._psf_file) as f:

            # Read the next line.
            line = f.readline()

            # There are donor records.
            if "!NDON" in line:
                has_donors = True

            # There are acceptor records.
            elif "NACC" in line:
                has_acceptors = True

            # There are non-bonded exclusion records.
            elif "NNB" in line:
                has_non_bonded = True

        # Append empty donor record.
        if not has_donors:
            f = open(self._psf_file, "a")
            f.write("\n%8d !NDON: donors\n" % 0)
            f.close()

        # Append empty acceptor record.
        if not has_acceptors:
            f = open(self._psf_file, "a")
            f.write("\n%8d !NACC: acceptors\n" % 0)
            f.close()

        # Append empty non-bonded exclusion record.
        if not has_acceptors:
            f = open(self._psf_file, "a")
            f.write("\n%8d !NNB: excluded\n" % 0)
            f.close()

        # Generate the NAMD configuration file.
        # Skip if the user has passed a custom config.
        if not self._is_custom:
            self._generate_config()
            self.writeConfig(self._config_file)

        # Return the list of input files.
        return self._input_files

    def _generate_config(self):
        """Generate NAMD configuration file strings."""

        # Clear the existing configuration list.
        self._config = []

        # Flag that the system doesn't contain a box.
        has_box = False

        # Check whether the system contains periodic box information.
        if 'space' in self._system.propertyKeys():
            # Flag that we have found a box.
            has_box = True

            # Get the box size.
            box_size = self._system.property('space').dimensions()

            # Since the box is translationally invariant, we set the cell
            # origin to be the average of the atomic coordinates. This
            # ensures a consistent wrapping for coordinates in the  NAMD
            # output files.
            origin = tuple(process._getAABox(self._system).center())

        # Work out the box from the atomic coordinates.
        elif not self._protocol.gas_phase:
            # Get the axis-aligned bounding box for the molecular system.
            aabox = process._getAABox(self._system)

            # Work out the box size and origin.
            box_size = 2 * aabox.halfExtents()
            origin = tuple(aabox.center())

        # Append generic configuration variables.

        # Topology.
        self.addToConfig("structure             %s" % path.basename(self._psf_file))
        self.addToConfig("coordinates           %s" % path.basename(self._pdb_file))

        # Velocities.
        if not self._velocity_file is None:
            self.addToConfig("velocities            %s" % path.basename(self._velocity_file))

        # Parameters.
        if self._is_charmm_params:
            self.addToConfig("paraTypeCharmm        on")
        self.addToConfig("parameters            %s" % path.basename(self._param_file))

        # Random number seed.
        if not self._seed is None:
            self.addToConfig("seed                  %d" % self._seed)

        # Exclusion policy.
        self.addToConfig("exclude               scaled1-4")

        # Non-bonded potential parameters.

        # Gas phase.
        if self._protocol.gas_phase:
            self.addToConfig("cutoff                999.")
            self.addToConfig("zeroMomentum          yes")
            self.addToConfig("switching             off")

        # Solvated.
        else:
            # Only use a cutoff if the box is large enough.
            if min(box_size) > 26:
                self.addToConfig("cutoff                12.")
                self.addToConfig("pairlistdist          14.")
                self.addToConfig("switching             on")
                self.addToConfig("switchdist            10.")

            # Load the XSC file.
            if has_box:
                self.addToConfig("extendedSystem        %s.xsc" % self._name)
                # We force the cell origin to be located at the system's centre
                # of geometry. This ensures a consistent periodic wrapping for
                # all NAMD output.
                self.addToConfig("cellOrigin            %.1f   %.1f   %.1f" % origin)

            # Set the cell using the axis-aligned bounding box.
            else:
                self.addToConfig("cellBasisVector1     %.1f   0.    0." % box_size[0])
                self.addToConfig("cellBasisVector2      0.   %.1f   0." % box_size[1])
                self.addToConfig("cellBasisVector3      0.    0.   %.1f" % box_size[2])
                self.addToConfig("cellOrigin            %.1f   %.1f   %.1f" % origin)

            # Wrap all molecular coordinates to the periodic box.
            self.addToConfig("wrapAll               on")

            # Periodic electrostatics.
            self.addToConfig("PME                   yes")
            self.addToConfig("PMEGridSpacing        1.")

        # Output file parameters.
        self.addToConfig("outputName            %s_out" % self._name)
        self.addToConfig("binaryOutput          no")
        self.addToConfig("binaryRestart         no")

        # Output frequency.
        self.addToConfig("restartfreq           500")
        self.addToConfig("xstFreq               500")

        # Printing frequency.
        self.addToConfig("outputEnergies        100")
        self.addToConfig("outputTiming          1000")

        # Add configuration variables for a minimisation simulation.
        if self._protocol.type() == ProtocolType.MINIMISATION:
            self.addToConfig("temperature           300")

            # Work out the number of steps. This must be a multiple of
            # stepspercycle, which is set the default of 20.
            steps = 20 * ceil(self._protocol.steps / 20)
            self.addToConfig("minimize              %d" % steps)

        # Add configuration variables for an equilibration simulation.
        elif self._protocol.type() == ProtocolType.EQUILIBRATION:
            # Set the Tcl temperature variable.
            if self._protocol.isConstantTemp():
                self.addToConfig("set temperature       %.2f" % self._protocol.temperature_start)
            else:
                # Cannot have a target temperature of exactly zero.
                if self._protocol.temperature_end == approx(0):
                    self.addToConfig("set temperature       0.01")
                else:
                    self.addToConfig("set temperature       %.2f" % self._protocol.temperature_end)
            self.addToConfig("temperature           $temperature")

            # Integrator parameters.
            self.addToConfig("timestep              2.")
            self.addToConfig("rigidBonds            all")
            self.addToConfig("nonbondedFreq         1")
            self.addToConfig("fullElectFrequency    2")

            # Constant temperature control.
            self.addToConfig("langevin              on")
            self.addToConfig("langevinDamping       1.")
            self.addToConfig("langevinTemp          $temperature")
            self.addToConfig("langevinHydrogen      no")

            # Constant pressure control.
            if not self._protocol.isConstantTemp():
                self.addToConfig("langevinPiston        on")
                self.addToConfig("langevinPistonTarget  1.01325")
                self.addToConfig("langevinPistonPeriod  100.")
                self.addToConfig("langevinPistonDecay   50.")
                self.addToConfig("langevinPistonTemp    $temperature")
                self.addToConfig("useGroupPressure      yes")
                self.addToConfig("useFlexibleCell       no")
                self.addToConfig("useConstantArea       no")

            # Restrain the backbone.
            if self._protocol.is_restrained:
                # Create a restrained system.
                restrained = process._restrain_backbone(self._system)

                # Create a PDB object, mapping the "occupancy" property to "restrained".
                p = PDB2(restrained, {"occupancy" : "restrained"})

                # File name for the restraint file.
                self._restraint_file = "%s/%s.restrained" % (self._work_dir, self._name)

                # Write the PDB file.
                p.writeToFile(self._restraint_file)

                # Update the configuration file.
                self.addToConfig("fixedAtoms            yes")
                self.addToConfig("fixedAtomsFile        %s.restrained" % self._name)

            # Work out number of steps needed to exceed desired running time,
            # rounded up to the nearest 20.
            steps = ceil(self._protocol.runtime / 2e-6)
            steps = 20 * ceil(steps / 20)

            # Heating/cooling simulation.
            if not self._protocol.isConstantTemp():
                # Work out temperature step size (assuming a unit increment).
                denom = abs(self._protocol.temperature_end - self._protocol.temperature_start)
                freq = floor(steps / denom)

                self.addToConfig("reassignFreq          %d" % freq)
                self.addToConfig("reassignTemp          %.2f" % self._protocol.temperature_start)
                self.addToConfig("reassignIncr          1.")
                self.addToConfig("reassignHold          %.2f" % self._protocol.temperature_end)

            # Run the simulation.
            self.addToConfig("run                   %d" % steps)

        # Add configuration variables for a production simulation.
        elif self._protocol.type() == ProtocolType.PRODUCTION:
            # Set the Tcl temperature variable.
            self.addToConfig("set temperature       %.2f" % self._protocol.temperature)
            self.addToConfig("temperature           $temperature")

            # Integrator parameters.
            self.addToConfig("timestep              2.")
            if self._protocol.first_step is not 0:
                self.addToConfig("firsttimestep         %d" % self._protocol.first_step)
            self.addToConfig("rigidBonds            all")
            self.addToConfig("nonbondedFreq         1")
            self.addToConfig("fullElectFrequency    2")

            # Constant temperature control.
            self.addToConfig("langevin              on")
            self.addToConfig("langevinDamping       1.")
            self.addToConfig("langevinTemp          $temperature")
            self.addToConfig("langevinHydrogen      no")

            # Constant pressure control.
            if self._protocol.ensemble is 'NPT':
                self.addToConfig("langevinPiston        on")
                self.addToConfig("langevinPistonTarget  1.01325")
                self.addToConfig("langevinPistonPeriod  100.")
                self.addToConfig("langevinPistonDecay   50.")
                self.addToConfig("langevinPistonTemp    $temperature")
                self.addToConfig("useGroupPressure      yes")
                self.addToConfig("useFlexibleCell       no")
                self.addToConfig("useConstantArea       no")

            # Work out number of steps needed to exceed desired running time,
            # rounded up to the nearest 20.
            steps = ceil(self._protocol.runtime / 2e-6)
            steps = 20 * ceil(steps / 20)

            # Trajectory output frequency.
            self.addToConfig("DCDfreq               %d" % floor(steps / self._protocol.frames))

            # Run the simulation.
            self.addToConfig("run                   %d" % steps)

    def start(self):
        """Start the NAMD simulation."""

        # Store the current working directory.
        dir = getcwd()

        # Change to the working directory for the process.
        # This avoid problems with relative paths.
        chdir(self._work_dir)

        # Write the command-line process to a README.txt file.
        with open("README.txt", "w") as f:

            # Set the command-line string.
            self._command = "%s %s.namd" % (self._exe, self._name)

            # Write the command to file.
            f.write("# NAMD was run with the following command:\n")
            f.write("%s\n" % self._command)

        # Start the timer.
        self._timer = timer()

        # Start the simulation.
        self._process = Process.run(self._exe, "%s.namd" % self._name,
                                               "%s.out"  % self._name,
                                               "%s.err"  % self._name)

        # Change back to the original working directory.
        chdir(dir)

    def getSystem(self):
        """Get the latest molecular configuration as a Sire system."""

        # Read the PDB coordinate file and construct a parameterised molecular
        # system using the original PSF and param files.

        is_coor = False

        # First check for final configuration.
        if path.isfile("%s/%s_out.coor" % (self._work_dir, self._name)):
            coor_file = "%s/%s_out.coor" % (self._work_dir, self._name)
            is_coor = True

        # Otherwise check for a restart file.
        elif path.isfile("%s/%s_out.restart.coor" % (self._work_dir, self._name)):
            coor_file = "%s/%s_out.restart.coor" % (self._work_dir, self._name)
            is_coor = True

        # Try to find an XSC file.

        is_xsc = False

        # First check for final XSC file.
        if path.isfile("%s/%s_out.xsc" % (self._work_dir, self._name)):
            xsc_file = "%s/%s_out.xsc" % (self._work_dir, self._name)
            is_xsc = True

        # Otherwise check for a restart XSC file.
        elif path.isfile("%s/%s_out.restart.xsc" % (self._work_dir, self._name)):
            xsc_file = "%s/%s_out.restart.xsc" % (self._work_dir, self._name)
            is_xsc = True

        # We found a coordinate file.
        if is_coor:
            # List of files.
            files = [ coor_file, self._psf_file, self._param_file ]

            # Add the box information.
            if is_xsc:
                files.append(xsc_file)

            # Create and return the molecular system.
            return MoleculeParser.read(files)

        else:
            return None

    def getRecord(self, record, time_series=False):
        """Get a record from the stdout dictionary.

           Keyword arguments:

           record      -- The record keyword.
           time_series -- Whether to return a list of time series records.
        """
        self.stdout(0)
        return self._get_stdout_record(record, time_series)

    def getRecords(self):
        """Return the dictionary of stdout time-series records."""
        return self._stdout_dict

    def getTime(self, time_series=False):
        """Get the time (in nanoseconds)."""

        if self._protocol.type == ProtocolType.MINIMISATION:
            return None

        else:
            # Get the list of time steps.
            self.stdout(0)
            time_steps = self._get_stdout_record('TS', time_series)

            # Multiply by the integration time step (2fs).
            if time_steps is not None:
                if time_series:
                    return [x * 2e-6 for x in time_steps]
                else:
                    return 2e-6 * time_steps

    def getStep(self, time_series=False):
        """Get the number of integration steps."""
        self.stdout(0)
        return self._get_stdout_record('TS', time_series)

    def getBondEnergy(self, time_series=False):
        """Get the bond energy."""
        self.stdout(0)
        return self._get_stdout_record('BOND', time_series)

    def getAngleEnergy(self, time_series=False):
        """Get the angle energy."""
        self.stdout(0)
        return self._get_stdout_record('ANGLE', time_series)

    def getDihedralEnergy(self, time_series=False):
        """Get the dihedral energy."""
        self.stdout(0)
        return self._get_stdout_record('DIHED', time_series)

    def getImproperEnergy(self, time_series=False):
        """Get the improper energy."""
        self.stdout(0)
        return self._get_stdout_record('IMPRP', time_series)

    def getElectrostaticEnergy(self, time_series=False):
        """Get the electrostatic energy."""
        self.stdout(0)
        return self._get_stdout_record('ELECT', time_series)

    def getVanDerWaalsEnergy(self, time_series=False):
        """Get the Van der Vaals energy."""
        self.stdout(0)
        return self._get_stdout_record('VDW', time_series)

    def getBoundaryEnergy(self, time_series=False):
        """Get the boundary energy."""
        self.stdout(0)
        return self._get_stdout_record('BOUNDARY', time_series)

    def getMiscEnergy(self, time_series=False):
        """Get the external energy."""
        self.stdout(0)
        return self._get_stdout_record('MISC', time_series)

    def getKineticEnergy(self, time_series=False):
        """Get the kinetic energy."""
        self.stdout(0)
        return self._get_stdout_record('KINETIC', time_series)

    def getPotentialEnergy(self, time_series=False):
        """Get the potential energy."""
        self.stdout(0)
        return self._get_stdout_record('POTENTIAL', time_series)

    def getTotalEnergy(self, time_series=False):
        """Get the total energy."""
        self.stdout(0)
        return self._get_stdout_record('TOTAL', time_series)

    def getTotal2Energy(self, time_series=False):
        """Get the total energy. (Better KE conservation.)"""
        self.stdout(0)
        return self._get_stdout_record('TOTAL2', time_series)

    def getTotal3Energy(self, time_series=False):
        """Get the total energy. (Smaller short-time fluctuations.)"""
        self.stdout(0)
        return self._get_stdout_record('TOTAL3', time_series)

    def getTemperature(self, time_series=False):
        """Get the temperature."""
        self.stdout(0)
        return self._get_stdout_record('TEMP', time_series)

    def getTemperatureAverage(self, time_series=False):
        """Get the average temperature."""
        self.stdout(0)
        return self._get_stdout_record('TEMPAVG', time_series)

    def getPressure(self, time_series=False):
        """Get the pressure."""
        self.stdout(0)
        return self._get_stdout_record('PRESSURE', time_series)

    def getPressureAverage(self, time_series=False):
        """Get the average pressure."""
        self.stdout(0)
        return self._get_stdout_record('PRESSAVG', time_series)

    def getGPressure(self, time_series=False):
        """Get the pressure. (Hydrogens incorporated into bonded atoms.)"""
        self.stdout(0)
        return self._get_stdout_record('GPRESSURE', time_series)

    def getGPressureAverage(self, time_series=False):
        """Get the average pressure. (Hydrogens incorporated into bonded atoms.)"""
        self.stdout(0)
        return self._get_stdout_record('GPRESSAVG', time_series)

    def getVolume(self, time_series=False):
        """Get the volume."""
        self.stdout(0)
        return self._get_stdout_record('VOLUME', time_series)

    def eta(self):
        """Get the estimated time for the process to finish (in minutes)."""

        # Make sure the list of stdout records is up to date.
        # Print the last zero lines, i.e. no output.
        self.stdout(0)

        # Now search backwards through the list to find the last TIMING record.
        for x, record in reversed(list(enumerate(self._stdout))):

            # Split the record using whitespace.
            data = record.split()

            # We've found a TIMING record.
            if len(data) > 0:
                if data[0] == "TIMING:":

                    # Try to find the "hours" record.
                    # If found, return the entry preceeding it.
                    try:
                        return float(data[data.index("hours") - 1]) * 60

                    # No record found.
                    except ValueError:
                        return None

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

            # Split the record using whitespace.
            data = self._stdout[-1].split()

            # Make sure there is at least one record.
            if len(data) > 0:

                # Store the updated energy title.
                if data[0] == 'ETITLE:':
                    self._stdout_title = data[1:]

                # This is an energy record.
                elif data[0] == 'ENERGY:':
                    # Extract the data.
                    stdout_data = data[1:]

                    # Add the records to the dictionary.
                    if (len(stdout_data) == len(self._stdout_title)):
                        for title, data in zip(self._stdout_title, stdout_data):
                            self._stdout_dict[title] = data

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

    def _get_stdout_record(self, key, time_series=False):
        """Helper function to get a stdout record from the dictionary."""

        # No data!
        if len(self._stdout_dict) is 0:
            return None

        if type(time_series) is not bool:
            warn("Non-boolean time-series flag. Defaulting to False!")
            time_series = False

        # Return the list of dictionary values.
        if time_series:
            try:
                if key is 'TS':
                    return [int(x) for x in self._stdout_dict[key]]
                else:
                    return [float(x) for x in self._stdout_dict[key]]

            except KeyError:
                return None

        # Return the most recent dictionary value.
        else:
            try:
                if key is 'TS':
                    return int(self._stdout_dict[key][-1])
                else:
                    return float(self._stdout_dict[key][-1])

            except KeyError:
                return None
