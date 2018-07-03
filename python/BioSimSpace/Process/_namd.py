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
Functionality for running simulations using NAMD.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

import Sire as _Sire

from . import _process
from .._SireWrappers import System as _System
from ..Trajectory import Trajectory as _Trajectory

import BioSimSpace.Protocol as _Protocol

import math as _math
import os as _os
import pygtail as _pygtail
import timeit as _timeit
import warnings as _warnings

__all__ = ["Namd"]

class Namd(_process.Process):
    """A class for running simulations using NAMD."""

    def __init__(self, system, protocol, exe=None,
            name="namd", work_dir=None, seed=None, map={}):
        """Constructor.

           Positional arguments:

           system        -- The molecular system.
           protocol      -- The protocol for the NAMD process.

           Keyword arguments:

           exe           -- The full path to the NAMD executable.
           name          -- The name of the process.
           work_dir      -- The working directory for the process.
           seed          -- A random number seed.
           map           -- A dictionary that maps system "properties" to their user defined
                            values. This allows the user to refer to properties with their
                            own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(system, protocol, name, work_dir, seed, map)

        # Set the package name.
        self._package_name = "NAMD"

        # This process can generate trajectory data.
        self._has_trajectory = True

        # If the path to the executable wasn't specified, then search
        # for it in $PATH.
        if exe is None:
            self._exe = _Sire.Base.findExe("namd2").absoluteFilePath()

        else:
            # Make sure executable exists.
            if _os.path.isfile(exe):
                self._exe = exe
            else:
                raise IOError("NAMD executable doesn't exist: '%s'" % exe)

        # Initialise the stdout dictionary and title header.
        self._stdout_dict = _process._MultiDict()
        self._stdout_title = None

        # The names of the input files.
        self._psf_file = "%s/%s.psf" % (self._work_dir, name)
        self._top_file = "%s/%s.pdb" % (self._work_dir, name)
        self._param_file = "%s/%s.params" % (self._work_dir, name)
        self._velocity_file = None
        self._restraint_file = None

        # Set the path for the NAMD configuration file.
        self._config_file = "%s/%s.namd" % (self._work_dir, name)

        # Create the list of input files.
        self._input_files = [self._config_file, self._psf_file, self._top_file, self._param_file]

        # Now set up the working directory for the process.
        self._setup()

    def _setup(self):
        """Setup the input files and working directory ready for simulation."""

        # Create the input files...

        # PSF and parameter files.
        psf = _Sire.IO.CharmmPSF(self._system)
        psf.writeToFile(self._psf_file)

        # PDB file.
        pdb = _Sire.IO.PDB2(self._system)
        pdb.writeToFile(self._top_file)

        # Try to write a PDB "velocity" restart file.
        # The file will only be generated if all atoms in self._system have
        # a "velocity" property.

        # First generate a name for the velocity file.
        velocity_file = _os.path.splitext(self._top_file)[0] + ".vel"

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

        # When converting forcefields, improper terms may be represented as
        # dihedrals (cosine impropers). As such, there will be no improper
        # records in the PSF file.
        has_impropers = False

        # Open the PSF file for reading.
        with open(self._psf_file) as f:

            # Loop over all lines.
            for line in f:

                # There are improper records.
                if "!NIMPHI" in line:
                    has_impropers = True

                # There are donor records.
                if "!NDON" in line:
                    has_donors = True

                # There are acceptor records.
                elif "NACC" in line:
                    has_acceptors = True

                # There are non-bonded exclusion records.
                elif "NNB" in line:
                    has_non_bonded = True

        # Append empty improper record.
        if not has_impropers:
            f = open(self._psf_file, "a")
            f.write("\n%8d !NIMPHI: impropers\n" % 0)
            f.close()

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
        if type(self._protocol) is _Protocol.Custom:
            self.setConfig(self._protocol.getConfig())
        else:
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

        if "space" in self._map:
            prop = self._map["space"]
        else:
            prop = "space"

        # Check whether the system contains periodic box information.
        if prop in self._system.propertyKeys():
            # Flag that we have found a box.
            has_box = True

            # Get the box size.
            box_size = self._system.property(prop).dimensions()

            # Since the box is translationally invariant, we set the cell
            # origin to be the average of the atomic coordinates. This
            # ensures a consistent wrapping for coordinates in the  NAMD
            # output files.
            origin = tuple(_process._getAABox(self._system).center())

        # No box information. Assume this is a gas phase simulation.
        else:
            _warnings.warn("No simulation box found. Assuming gas phase simulation.")
            has_box = False

        if "param_format" in self._map:
            prop = self._map["param_format"]
        else:
            prop = "param_format"

        # Check whether the system contains parameter format information.
        if prop in self._system.propertyKeys():
            # Get the parameter format.
            if self._system.property(prop).toString() == "CHARMM":
                is_charmm = True
            else:
                is_charmm = False

        # No parameter format information. Assume these are CHARMM parameters.
        else:
            _warnings.warn("No parameter format information found. Assuming CHARMM parameters.")
            is_charmm = True

        # Append generic configuration variables.

        # Topology.
        self.addToConfig("structure             %s" % _os.path.basename(self._psf_file))
        self.addToConfig("coordinates           %s" % _os.path.basename(self._top_file))

        # Velocities.
        if self._velocity_file is not None:
            self.addToConfig("velocities            %s" % _os.path.basename(self._velocity_file))

        # Parameters.
        if is_charmm:
            self.addToConfig("paraTypeCharmm        on")
        self.addToConfig("parameters            %s" % _os.path.basename(self._param_file))

        # Random number seed.
        if self._is_seeded:
            self.addToConfig("seed                  %d" % self._seed)

        # Exclusion policy.
        self.addToConfig("exclude               scaled1-4")

        # Non-bonded potential parameters.

        # Gas phase.
        if not has_box:
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
            self.addToConfig("extendedSystem        %s.xsc" % self._name)
            # We force the cell origin to be located at the system's centre
            # of geometry. This ensures a consistent periodic wrapping for
            # all NAMD output.
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
        if type(self._protocol) is _Protocol.Minimisation:
            self.addToConfig("temperature           300")

            # Work out the number of steps. This must be a multiple of
            # stepspercycle, which is set the default of 20.
            steps = 20 * _math.ceil(self._protocol.getSteps() / 20)
            self.addToConfig("minimize              %d" % steps)

        # Add configuration variables for an equilibration simulation.
        if type(self._protocol) is _Protocol.Equilibration:
            # Set the Tcl temperature variable.
            if self._protocol.isConstantTemp():
                self.addToConfig("set temperature       %.2f"
                    % self._protocol.getStartTemperature().kelvin().magnitude())
            else:
                self.addToConfig("set temperature       %.2f"
                    % self._protocol.getEndTemperature().kelvin().magnitude())
            self.addToConfig("temperature           $temperature")

            # Integrator parameters.
            self.addToConfig("timestep              %.2f"
                % self._protocol.getTimeStep().femtoseconds().magnitude())
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
            if self._protocol.isRestrained():
                # Create a restrained system.
                restrained = _process._restrain_backbone(self._system)

                # Create a PDB object, mapping the "occupancy" property to "restrained".
                if "occupancy" in self._map:
                    prop = self._map["occupancy"]
                else:
                    prop = "occupancy"
                p = PDB2(restrained, {prop : "restrained"})

                # File name for the restraint file.
                self._restraint_file = "%s/%s.restrained" % (self._work_dir, self._name)

                # Write the PDB file.
                p.writeToFile(self._restraint_file)

                # Update the configuration file.
                self.addToConfig("fixedAtoms            yes")
                self.addToConfig("fixedAtomsFile        %s.restrained" % self._name)

            # Work out number of steps needed to exceed desired running time,
            # rounded up to the nearest 20.
            steps = _math.ceil(self._protocol.getRunTime().nanoseconds().magnitude() /
                               self._protocol.getTimeStep().nanoseconds().magnitude())
            steps = 20 * _math.ceil(steps / 20)

            # Heating/cooling simulation.
            if not self._protocol.isConstantTemp():
                # Work out temperature step size (assuming a unit increment).
                denom = abs(self._protocol.getEndTemperature().kelvin().magnitude() -
                            self._protocol.getStartTemperature().kelvin().magnitude())
                freq = _math.floor(steps / denom)

                self.addToConfig("reassignFreq          %d" % freq)
                self.addToConfig("reassignTemp          %.2f"
                    % self._protocol.getStartTemperature().kelvin().magnitude())
                self.addToConfig("reassignIncr          1.")
                self.addToConfig("reassignHold          %.2f"
                    % self._protocol.getEndTemperature().kelvin().magnitude())

            # Trajectory output frequency.
            self.addToConfig("DCDfreq               %d" % _math.floor(steps / self._protocol.getFrames()))

            # Run the simulation.
            self.addToConfig("run                   %d" % steps)

        # Add configuration variables for a production simulation.
        elif type(self._protocol) is _Protocol.Production:
            # Set the Tcl temperature variable.
            self.addToConfig("set temperature       %.2f"
                % self._protocol.getTemperature().kelvin().magnitude())
            self.addToConfig("temperature           $temperature")

            # Integrator parameters.
            self.addToConfig("timestep              %.2f"
                % self._protocol.getTimeStep().femtoseconds().magnitude())
            if self._protocol.getFirstStep() != 0:
                self.addToConfig("firsttimestep         %d" % self._protocol.getFirstStep())
            self.addToConfig("rigidBonds            all")
            self.addToConfig("nonbondedFreq         1")
            self.addToConfig("fullElectFrequency    2")

            # Constant temperature control.
            self.addToConfig("langevin              on")
            self.addToConfig("langevinDamping       1.")
            self.addToConfig("langevinTemp          $temperature")
            self.addToConfig("langevinHydrogen      no")

            # Constant pressure control.
            if self._protocol.getEnsemble() == "NPT":
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
            steps = _math.ceil(self._protocol.getRunTime().nanoseconds().magnitude() /
                               self._protocol.getTimeStep().nanoseconds().magnitude())
            steps = 20 * _math.ceil(steps / 20)

            # Trajectory output frequency.
            self.addToConfig("DCDfreq               %d" % _math.floor(steps / self._protocol.getFrames()))

            # Run the simulation.
            self.addToConfig("run                   %d" % steps)

    def start(self):
        """Start the NAMD process."""

        # Process is already running.
        if self._process is not None:
            if self._process.isRunning():
                return

        # Store the current working directory.
        dir = _os.getcwd()

        # Change to the working directory for the process.
        # This avoid problems with relative paths.
        _os.chdir(self._work_dir)

        # Write the command-line process to a README.txt file.
        with open("README.txt", "w") as f:

            # Set the command-line string.
            self._command = "%s %s.namd" % (self._exe, self._name)

            # Write the command to file.
            f.write("# NAMD was run with the following command:\n")
            f.write("%s\n" % self._command)

        # Start the timer.
        self._timer = _timeit.default_timer()

        # Start the simulation.
        self._process = _Sire.Base.Process.run(self._exe, "%s.namd" % self._name,
                                                          "%s.out"  % self._name,
                                                          "%s.err"  % self._name)

        # Change back to the original working directory.
        _os.chdir(dir)

        return self

    def getSystem(self, block="AUTO"):
        """Get the latest molecular configuration as a Sire system.

           Keyword arguments:

           block -- Whether to block until the process has finished running.
        """

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        # Read the PDB coordinate file and construct a parameterised molecular
        # system using the original PSF and param files.

        has_coor = False

        # First check for final configuration.
        if _os.path.isfile("%s/%s_out.coor" % (self._work_dir, self._name)):
            coor_file = "%s/%s_out.coor" % (self._work_dir, self._name)
            has_coor = True

        # Otherwise check for a restart file.
        elif _os.path.isfile("%s/%s_out.restart.coor" % (self._work_dir, self._name)):
            coor_file = "%s/%s_out.restart.coor" % (self._work_dir, self._name)
            has_coor = True

        # Try to find an XSC file.

        has_xsc = False

        # First check for final XSC file.
        if _os.path.isfile("%s/%s_out.xsc" % (self._work_dir, self._name)):
            xsc_file = "%s/%s_out.xsc" % (self._work_dir, self._name)
            has_xsc = True

        # Otherwise check for a restart XSC file.
        elif _os.path.isfile("%s/%s_out.restart.xsc" % (self._work_dir, self._name)):
            xsc_file = "%s/%s_out.restart.xsc" % (self._work_dir, self._name)
            has_xsc = True

        # We found a coordinate file.
        if has_coor:
            # List of files.
            files = [ coor_file, self._psf_file, self._param_file ]

            # Add the box information.
            if has_xsc:
                files.append(xsc_file)

            # Create and return the molecular system.
            return _System(_Sire.IO.MoleculeParser.read(files))

        else:
            return None

    def getCurrentSystem(self):
        """Get the latest molecular configuration as a Sire system."""
        return self.getSystem(block=False)

    def getTrajectory(self, block="AUTO"):
        """Return a trajectory object."""

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        return _Trajectory(process=self)

    def getRecord(self, record, time_series=False, block="AUTO"):
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
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        self.stdout(0)
        return self._get_stdout_record(record, time_series)

    def getCurrentRecord(self, record, time_series=False):
        """Get a current record from the stdout dictionary.

           Positional arguments:

           record      -- The record key.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        self.stdout(0)
        return self._get_stdout_record(record, time_series)

    def getRecords(self, block="AUTO"):
        """Return the dictionary of stdout time-series records.

           Keyword arguments:

           block       -- Whether to block until the process has finished running.
        """
        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        return self._stdout_dict.copy()

    def getCurrentRecords(self):
        """Return the current dictionary of stdout time-series records."""
        return getRecords(block=False)

    def getTime(self, time_series=False, block="AUTO"):
        """Get the time (in nanoseconds).

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """

        if type(self._protocol) is _Protocol.Minimisation:
            return None

        else:
            # Get the list of time steps.
            time_steps = self.getRecord("TS", time_series, block)

            # Convert the time step to nanoseconds.
            timestep = self._protocol.getTimeStep().nanoseconds().magnitude()

            # Multiply by the integration time step.
            if time_steps is not None:
                if time_series:
                    return [x * timestep for x in time_steps]
                else:
                    return timestep * time_steps

    def getCurrentTime(self, time_series=False):
        """Get the current time (in nanoseconds).

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getTime(time_series, block=False)

    def getStep(self, time_series=False, block="AUTO"):
        """Get the number of integration steps.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord("TS", time_series, block)

    def getCurrentStep(self, time_series=False):
        """Get the current number of integration steps.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getStep(time_series, block=False)

    def getBondEnergy(self, time_series=False, block="AUTO"):
        """Get the bond energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord("BOND", time_series, block)

    def getCurrentBondEnergy(self, time_series=False):
        """Get the current bond energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getBondEnergy(time_series, block=False)

    def getAngleEnergy(self, time_series=False, block="AUTO"):
        """Get the angle energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord("ANGLE", time_series, block)

    def getCurrentAngleEnergy(self, time_series=False):
        """Get the current angle energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getAngleEnergy(time_series, block=False)

    def getDihedralEnergy(self, time_series=False, block="AUTO"):
        """Get the dihedral energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord("DIHED", time_series, block)

    def getCurrentDihedralEnergy(self, time_series=False):
        """Get the current dihedral energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getDihedralEnergy(time_series, block=False)

    def getImproperEnergy(self, time_series=False, block="AUTO"):
        """Get the improper energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord("IMPRP", time_series, block)

    def getCurrentImproperEnergy(self, time_series=False):
        """Get the current improper energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getImproperEnergy(time_series, block=False)

    def getElectrostaticEnergy(self, time_series=False, block="AUTO"):
        """Get the electrostatic energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord("ELECT", time_series, block)

    def getCurrentElectrostaticEnergy(self, time_series=False):
        """Get the current electrostatic energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getElectrostaticEnergy(time_series, block=False)

    def getVanDerWaalsEnergy(self, time_series=False, block="AUTO"):
        """Get the Van der Vaals energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord("VDW", time_series)

    def getCurrentVanDerWaalsEnergy(self, time_series=False):
        """Get the current Van der Waals energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getVanDerWaalsEnergy(time_series, block=False)

    def getBoundaryEnergy(self, time_series=False, block="AUTO"):
        """Get the boundary energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord("BOUNDARY", time_series, block)

    def getCurrentBoundaryEnergy(self, time_series=False):
        """Get the current boundary energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getBoundaryEnergy(time_series, block=False)

    def getMiscEnergy(self, time_series=False, block="AUTO"):
        """Get the external energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord("MISC", time_series, block)

    def getCurrentMiscEnergy(self, time_series=False):
        """Get the current external energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getMiscEnergy(time_series, block=False)

    def getKineticEnergy(self, time_series=False, block="AUTO"):
        """Get the kinetic energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord("KINETIC", time_series, block)

    def getCurrentKineticEnergy(self, time_series=False):
        """Get the current kinetic energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getKineticEnergy(time_series, block=False)

    def getPotentialEnergy(self, time_series=False, block="AUTO"):
        """Get the potential energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord("POTENTIAL", time_series, block)

    def getCurrentPotentialEnergy(self, time_series=False):
        """Get the current potential energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getPotentialEnergy(time_series, block=False)

    def getTotalEnergy(self, time_series=False, block="AUTO"):
        """Get the total energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord("TOTAL", time_series, block)

    def getCurrentTotalEnergy(self, time_series=False):
        """Get the current potential energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getTotalEnergy(time_series, block=False)

    def getTotal2Energy(self, time_series=False, block="AUTO"):
        """Get the total energy. (Better KE conservation.)

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord("TOTAL2", time_series, block)

    def getCurrentTotal2Energy(self, time_series=False):
        """Get the current total energy. (Better KE conservation.)

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getTotal2Energy(time_series, block=False)

    def getTotal3Energy(self, time_series=False, block="AUTO"):
        """Get the total energy. (Smaller short-time fluctuations.)

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord("TOTAL3", time_series, block)

    def getCurrentTotal3Energy(self, time_series=False):
        """Get the total energy. (Smaller short-time fluctuations.)

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getTotal3Energy(time_series, block=False)

    def getTemperature(self, time_series=False, block="AUTO"):
        """Get the temperature.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord("TEMP", time_series, block)

    def getCurrentTemperature(self, time_series=False):
        """Get the temperature.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getTemperature(time_series, block=False)

    def getTemperatureAverage(self, time_series=False, block="AUTO"):
        """Get the average temperature.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord("TEMPAVG", time_series, block)

    def getCurrentTemperatureAverage(self, time_series=False):
        """Get the current average temperature.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getTemperatureAverage(time_series, block=False)

    def getPressure(self, time_series=False, block="AUTO"):
        """Get the pressure.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord("PRESSURE", time_series, block)

    def getCurrentPressure(self, time_series=False):
        """Get the current pressure.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getPressure(time_series, block=False)

    def getPressureAverage(self, time_series=False, block="AUTO"):
        """Get the average pressure.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord("PRESSAVG", time_series, block)

    def getCurrentPressureAverage(self, time_series=False):
        """Get the current average pressure.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getPressureAverage(time_series, block=False)

    def getGPressure(self, time_series=False, block="AUTO"):
        """Get the pressure. (Hydrogens incorporated into bonded atoms.)

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord("GPRESSURE", time_series, block)

    def getCurrentGPressure(self, time_series=False):
        """Get the current pressure. (Hydrogens incorporated into bonded atoms.)

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getGPressure(time_series, block=False)

    def getGPressureAverage(self, time_series=False, block="AUTO"):
        """Get the average pressure. (Hydrogens incorporated into bonded atoms.)

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord("GPRESSAVG", time_series, block)

    def getCurrentGPressureAverage(self, time_series=False):
        """Get the current average pressure. (Hydrogens incorporated into bonded atoms.)

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getGPressureAverage(time_series, block=False)

    def getVolume(self, time_series=False, block="AUTO"):
        """Get the volume.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord("VOLUME", time_series, block)

    def getCurrentVolume(self, time_series=False):
        """Get the current volume.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getVolume(time_series, block=False)

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
        for line in _pygtail.Pygtail(self._stdout_file):
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
        """Helper function to get a stdout record from the dictionary.

           Positional arguments:

           key         -- The record key.

           Keyword arguments:

           time_series -- Whether to return a time series of records.
        """

        # No data!
        if len(self._stdout_dict) is 0:
            return None

        if type(time_series) is not bool:
            _warnings.warn("Non-boolean time-series flag. Defaulting to False!")
            time_series = False

        # Return the list of dictionary values.
        if time_series:
            try:
                if key is "TS":
                    return [int(x) for x in self._stdout_dict[key]]
                else:
                    return [float(x) for x in self._stdout_dict[key]]

            except KeyError:
                return None

        # Return the most recent dictionary value.
        else:
            try:
                if key is "TS":
                    return int(self._stdout_dict[key][-1])
                else:
                    return float(self._stdout_dict[key][-1])

            except KeyError:
                return None

    def _get_trajectory_files(self):
        """Get all files associated with the molecular trajectory."""

        # Name of the trajectory file.
        traj_file = "%s/%s_out.dcd" % (self._work_dir, self._name)

        # Return the trajectory and topology file.
        if _os.path.isfile("%s/%s.nc" % (self._work_dir, self._name)):
            return (traj_file, self._top_file)

        # No trajectory file.
        else:
            return None
