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

class NamdProcess(process.Process):
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
        self._stdout_dict = process._MDict()
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
            self._generate_config_file()

        # Return the list of input files.
        return self._input_files

    def _generate_config_file(self):
        """Generate a NAMD configuration file."""

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

        # Open the configuration file for writing.
        f = open(self._config_file, "w")

        # Write generic configuration variables.

        # Topology.
        f.write("structure             %s\n" % path.basename(self._psf_file))
        f.write("coordinates           %s\n" % path.basename(self._pdb_file))

        # Velocities.
        if not self._velocity_file is None:
            f.write("velocities            %s\n" % path.basename(self._velocity_file))

        # Parameters.
        if self._is_charmm_params:
            f.write("paraTypeCharmm        on\n")
        f.write("parameters            %s\n" % path.basename(self._param_file))

        # Random number seed.
        if not self._seed is None:
            f.write("seed                  %s\n" % self._seed)

        # Exclusion policy.
        f.write("exclude               scaled1-4\n")

        # Non-bonded potential parameters.

        # Gas phase.
        if self._protocol.gas_phase:
            f.write("cutoff                999.\n")
            f.write("zeroMomentum          yes\n")
            f.write("switching             off\n")

        # Solvated.
        else:
            # Only use a cutoff if the box is large enough.
            if min(box_size) > 26:
                f.write("cutoff                12.\n")
                f.write("pairlistdist          14.\n")
                f.write("switching             on\n")
                f.write("switchdist            10.\n")

            # Load the XSC file.
            if has_box:
                f.write("extendedSystem        %s.xsc\n" % self._name)
                # We force the cell origin to be located at the system's centre
                # of geometry. This ensures a consistent periodic wrapping for
                # all NAMD output.
                f.write("cellOrigin            %.1f   %.1f   %.1f\n" % origin)

            # Set the cell using the axis-aligned bounding box.
            else:
                f.write("cellBasisVector1     %.1f   0.    0.\n" % box_size[0])
                f.write("cellBasisVector2      0.   %.1f   0.\n" % box_size[1])
                f.write("cellBasisVector3      0.    0.   %.1f\n" % box_size[2])
                f.write("cellOrigin            %.1f   %.1f   %.1f\n" % origin)

            # Wrap all molecular coordinates to the periodic box.
            f.write("wrapAll               on\n")

            # Periodic electrostatics.
            f.write("PME                   yes\n")
            f.write("PMEGridSpacing        1.\n")

        # Output file parameters.
        f.write("outputName            %s_out\n" % self._name)
        f.write("binaryOutput          no\n")
        f.write("binaryRestart         no\n")

        # Output frequency.
        f.write("restartfreq           500\n")
        f.write("xstFreq               500\n")

        # Printing frequency.
        f.write("outputEnergies        100\n")
        f.write("outputTiming          1000\n")

        # Add configuration variables for a minimisation simulation.
        if self._protocol.type() == ProtocolType.MINIMISATION:
            f.write("temperature           %s\n" % self._protocol.temperature)

            # Work out the number of steps. This must be a multiple of
            # stepspercycle, which is set the default of 20.
            steps = 20 * ceil(self._protocol.steps / 20)
            f.write("minimize              %s\n" % steps)

        # Add configuration variables for an equilibration simulation.
        elif self._protocol.type() == ProtocolType.EQUILIBRATION:
            # Set the Tcl temperature variable.
            if self._protocol.isConstantTemp():
                f.write("set temperature       %s\n" % self._protocol.temperature_start)
            else:
                f.write("set temperature       %s\n" % self._protocol.temperature_end)
            f.write("temperature           $temperature\n")

            # Integrator parameters.
            f.write("timestep              2.\n")
            f.write("rigidBonds            all\n")
            f.write("nonbondedFreq         1\n")
            f.write("fullElectFrequency    2\n")

            # Constant temperature control.
            f.write("langevin              on\n")
            f.write("langevinDamping       1.\n")
            f.write("langevinTemp          $temperature\n")
            f.write("langevinHydrogen      no\n")

            # Constant pressure control.
            if not self._protocol.isConstantTemp():
                f.write("langevinPiston        on\n")
                f.write("langevinPistonTarget  1.01325\n")
                f.write("langevinPistonPeriod  100.\n")
                f.write("langevinPistonDecay   50.\n")
                f.write("langevinPistonTemp    $temperature\n")
                f.write("useGroupPressure      yes\n")
                f.write("useFlexibleCell       no\n")
                f.write("useConstantArea       no\n")

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
                f.write("fixedAtoms            yes\n")
                f.write("fixedAtomsFile        %s.restrained\n" % self._name)

            # Work out number of steps needed to exceed desired running time,
            # rounded up to the nearest 20.
            steps = ceil(self._protocol.runtime / 2e-6)
            steps = 20 * ceil(steps / 20)

            # Heating/cooling simulation.
            if not self._protocol.isConstantTemp():
                # Work out temperature step size (assuming a unit increment).
                denom = abs(self._protocol.temperature_target - self._protocol.temperature_start)
                freq = floor(steps / denom)

                f.write("reassignFreq          %s\n" % freq)
                f.write("reassignTemp          %s\n" % self._protocol.temperature_start)
                f.write("reassignIncr          1.\n")
                f.write("reassignHold          %s\n" % self._protocol.temperature_target)

            # Run the simulation.
            f.write("run                   %s\n" % steps)

        # Add configuration variables for a production simulation.
        elif self._protocol.type() == ProtocolType.PRODUCTION:
            # Set the Tcl temperature variable.
            f.write("set temperature       %s\n" % self._protocol.temperature)
            f.write("temperature           $temperature\n")

            # Integrator parameters.
            f.write("timestep              2.\n")
            if self._protocol.first_step is not 0:
                f.write("firsttimestep         %s\n" % self._protocol.first_step)
            f.write("rigidBonds            all\n")
            f.write("nonbondedFreq         1\n")
            f.write("fullElectFrequency    2\n")

            # Constant temperature control.
            f.write("langevin              on\n")
            f.write("langevinDamping       1.\n")
            f.write("langevinTemp          $temperature\n")
            f.write("langevinHydrogen      no\n")

            # Constant pressure control.
            if self._protocol.ensemble is 'NPT':
                f.write("langevinPiston        on\n")
                f.write("langevinPistonTarget  1.01325\n")
                f.write("langevinPistonPeriod  100.\n")
                f.write("langevinPistonDecay   50.\n")
                f.write("langevinPistonTemp    $temperature\n")
                f.write("useGroupPressure      yes\n")
                f.write("useFlexibleCell       no\n")
                f.write("useConstantArea       no\n")

            # Work out number of steps needed to exceed desired running time,
            # rounded up to the nearest 20.
            steps = ceil(self._protocol.runtime / 2e-6)
            steps = 20 * ceil(steps / 20)

            # Trajectory output frequency.
            f.write("DCDfreq               %s\n" % floor(steps / self._protocol.frames))

            # Run the simulation.
            f.write("run                   %s\n" % steps)

        # Close the configuration file.
        f.close()

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
            f.write("# NamdProcess was run with the following command:\n")
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
                if data[0] == "ETITLE:":
                    self._stdout_title = data[1:]

                # This is an energy record.
                elif data[0] == "ENERGY:":
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
            time_seris = False

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
