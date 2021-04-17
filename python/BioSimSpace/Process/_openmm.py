######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2021
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BioSimSpace. If not, see <http://www.gnu.org/licenses/>.
#####################################################################

"""
Functionality for running simulations with OpenMM.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["OpenMM"]

import math as _math
import os as _os
import pygtail as _pygtail
import sys as _sys
import shutil as _shutil
import timeit as _timeit
import warnings as _warnings

from Sire import Base as _SireBase
from Sire import IO as _SireIO

from BioSimSpace import _isVerbose
from BioSimSpace._Exceptions import IncompatibleError as _IncompatibleError
from BioSimSpace._Exceptions import MissingSoftwareError as _MissingSoftwareError
from BioSimSpace._SireWrappers import System as _System
from BioSimSpace.Metadynamics import CollectiveVariable as _CollectiveVariable
from BioSimSpace.Types._type import Type as _Type

from BioSimSpace import IO as _IO
from BioSimSpace import Protocol as _Protocol
from BioSimSpace import Trajectory as _Trajectory
from BioSimSpace import Types as _Types
from BioSimSpace import Units as _Units
from BioSimSpace import _Utils as _Utils

from . import _process

class OpenMM(_process.Process):
    """A class for running simulations using OpenMM."""

    # Dictionary of platforms and their OpenMM keyword.
    _platforms = { "CPU"    : "CPU",
                   "CUDA"   : "CUDA",
                   "OPENCL" : "OpenCL" }

    def __init__(self, system, protocol, exe=None, name="openmm",
            platform="CPU", work_dir=None, seed=None, property_map={}):
        """Constructor.

           Parameters
           ----------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The molecular system.

           protocol : :class:`Protocol <BioSimSpace.Protocol>`
               The protocol for the OpenMM process.

           exe : str
               The full path to the Python interpreter used to run OpenMM.

           name : str
               The name of the process.

           platform : str
               The platform for the simulation: "CPU", "CUDA", or "OPENCL".
               For CUDA use the CUDA_VISIBLE_DEVICES environment variable to
               set the GPUs on which to run, e.g. to run on two GPUs indexed
               0 and 1 use: CUDA_VISIBLE_DEVICES=0,1. For OPENCL, instead use
               OPENCL_VISIBLE_DEVICES.

           work_dir :
               The working directory for the process.

           seed : int
               A random number seed.

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(system, protocol, name, work_dir, seed, property_map)

        # Set the package name.
        self._package_name = "OPENMM"

        # This process can generate trajectory data.
        self._has_trajectory = True

        # If the path a Python interpreter wasn't specified, then use the bundled sire_python.
        if exe is None:
            bin_dir = _SireBase.getBinDir()
            # Generate the name of the sire_python interpreter.
            if _sys.platform == "win32":
                self._exe = _os.path.join(_os.path.normpath(bin_dir), "sire_python.exe")
            else:
                self._exe = _os.path.join(_os.path.normpath(bin_dir), "sire_python")
        else:
            # Make sure executable exists.
            if _os.path.isfile(exe):
                self._exe = exe
            else:
                raise IOError("OpenMM Python interpreter doesn't exist: '%s'" % exe)

        if type(platform) is not str:
            raise TypeError("'platform' must be of type 'str'.")
        else:
            # Strip all whitespace and convert to upper case.
            platform = platform.replace(" ", "").upper()

            # Check for platform support.
            if platform not in self._platforms:
                raise ValueError("Supported platforms are: %s" % self._platforms.keys())
            else:
                self._platform = self._platforms[platform]

        # Initialise the stdout dictionary and title header.
        self._stdout_dict = _process._MultiDict()

        # Store the name of the OpenMM log file.
        self._log_file = "%s/%s.log" % (self._work_dir, name)

        # Initialise the log file separator.
        self._record_separator = None

        # Initialise a dictionary to map log file records to their column order.
        self._record_mapping = {}

        # The names of the input files. We choose to use AMBER files since they
        # are self-contained, but could equally work with GROMACS files.
        self._rst_file = "%s/%s.rst7" % (self._work_dir, name)
        self._top_file = "%s/%s.prm7" % (self._work_dir, name)

        # The name of the trajectory file.
        self._traj_file = "%s/%s.dcd" % (self._work_dir, name)

        # Set the path for the OpenMM Python script. (We use the concept of a
        # config file for consistency with other Process classes.)
        self._config_file = "%s/%s.py" % (self._work_dir, name)

        # Create the list of input files.
        self._input_files = [self._config_file, self._rst_file, self._top_file]

        # Now set up the working directory for the process.
        self._setup()

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.Process.%s: system=%s, protocol=%s, exe='%s', name='%s', platform='%s', work_dir='%s' seed=%s>" \
            % (self.__class__.__name__, str(self._system), self._protocol.__repr__(),
               self._exe + ("%s " % self._script if self._script else ""),
               self._name, self._platform, self._work_dir, self._seed)

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return "BioSimSpace.Process.%s(%s, %s, exe='%s', name='%s', platform='%s', work_dir='%s', seed=%s)" \
            % (self.__class__.__name__, str(self._system), self._protocol.__repr__(),
               self._exe + ("%s " % self._script if self._script else ""),
               self._name, self._platform, self._work_dir, self._seed)

    def _setup(self):
        """Setup the input files and working directory ready for simulation."""

        # Create the input files...

        # RST file (coordinates).
        try:
            rst = _SireIO.AmberRst7(self._system._sire_object, self._property_map)
            rst.writeToFile(self._rst_file)
        except Exception as e:
            msg = "Failed to write system to 'RST7' format."
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # PRM file (topology).
        try:
            prm = _SireIO.AmberPrm(self._system._sire_object, self._property_map)
            prm.writeToFile(self._top_file)

        except Exception as e:
            msg = "Failed to write system to 'PRM7' format."
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # Skip if the user has passed a custom config.
        if type(self._protocol) is _Protocol.Custom:
            self.setConfig(self._protocol.getConfig())
        else:
            self._generate_config()
        self.writeConfig(self._config_file)

        # Generate the dictionary of command-line arguments.
        self._generate_args()

        # Return the list of input files.
        return self._input_files

    def _generate_config(self):
        """Generate OpenMM Python script file strings."""

        # Clear the existing configuration list.
        self._config = []

        # Flag that this isn't a custom protocol.
        self._protocol._setCustomised(False)

        # Get the "space" property from the user mapping.
        prop = self._property_map.get("space", "space")

        # Check whether the system contains periodic box information.
        # For now, we'll not attempt to generate a box if the system property
        # is missing. If no box is present, we'll assume a non-periodic simulation.
        if prop in self._system._sire_object.propertyKeys():
            has_box = True
        else:
            _warnings.warn("No simulation box found. Assuming gas phase simulation.")
            has_box = False

        if type(self._protocol) is _Protocol.Minimisation:
            # Write the OpenMM import statements.
            self._add_config_imports()

            # Load the input files.
            self.addToConfig("\n# Load the topology and coordinate files.")
            self.addToConfig(f"prmtop = AmberPrmtopFile('{self._name}.prm7')")
            self.addToConfig(f"inpcrd = AmberInpcrdFile('{self._name}.rst7')")

            # Don't use a cut-off if this is a vacuum simulation or if box information
            # is missing.
            self.addToConfig("\n# Initialise the molecular system.")
            if not has_box or not self._has_water:
                self.addToConfig("system = prmtop.createSystem(nonbondedMethod=NoCutoff,")
            else:
                self.addToConfig("system = prmtop.createSystem(nonbondedMethod=PME,")
            self.addToConfig(    "                             nonbondedCutoff=1*nanometer,")
            self.addToConfig(    "                             constraints=HBonds)")

            # Set the integrator. (We use a zero-temperature Langevin integrator.)
            self.addToConfig("\n# Define the integrator.")
            self.addToConfig("integrator = LangevinIntegrator(0*kelvin,")
            self.addToConfig("                                1/picosecond,")
            self.addToConfig("                                0.002*picoseconds)")

            # Add the platform information.
            self._add_config_platform()

            # Set up the simulation object.
            self.addToConfig("\n# Initialise and configure the simulation object.")
            self.addToConfig("simulation = Simulation(prmtop.topology,")
            self.addToConfig("                        system,")
            self.addToConfig("                        integrator,")
            self.addToConfig("                        platform,")
            self.addToConfig("                        properties)")
            self.addToConfig("simulation.context.setPositions(inpcrd.positions)")
            self.addToConfig("simulation.minimizeEnergy()")

            # Add the reporters.
            self.addToConfig("\n# Add reporters.")
            if self._protocol.getSteps() < 100:
                self._add_config_reporters(state_interval=1, traj_interval=1)
            else:
                self._add_config_reporters(state_interval=100, traj_interval=100)

            # Now run the simulation.
            self.addToConfig("\n# Run the simulation.")
            self.addToConfig(f"simulation.step({self._protocol.getSteps()})")

        elif type(self._protocol) is _Protocol.Equilibration:
            # Write the OpenMM import statements and monkey-patches.
            self._add_config_imports()
            self._add_config_monkey_patches()

            # Load the input files.
            self.addToConfig("\n# Load the topology and coordinate files.")
            self.addToConfig(f"prmtop = AmberPrmtopFile('{self._name}.prm7')")
            self.addToConfig(f"inpcrd = AmberInpcrdFile('{self._name}.rst7')")

            # Don't use a cut-off if this is a vacuum simulation or if box information
            # is missing.
            is_periodic = True
            self.addToConfig("\n# Initialise the molecular system.")
            if not has_box or not self._has_water:
                is_periodic = False
                self.addToConfig("system = prmtop.createSystem(nonbondedMethod=NoCutoff,")
            else:
                self.addToConfig("system = prmtop.createSystem(nonbondedMethod=PME,")
            self.addToConfig(    "                             nonbondedCutoff=1*nanometer,")
            self.addToConfig(    "                             constraints=HBonds)")

            # Get the starting temperature and system pressure.
            temperature = self._protocol.getStartTemperature().kelvin().magnitude()
            pressure = self._protocol.getPressure()

            # Add a Monte Carlo barostat if the simulation is at constant pressure.
            is_const_pressure = False
            if pressure is not None:
                # Cannot use a barostat with a non-periodic system.
                if not is_periodic:
                    _warnings.warn("Cannot use a barostat for a vacuum or non-periodic simulation")
                else:
                    is_const_pressure = True

                    # Convert to bar and get the magnitude.
                    pressure = pressure.bar().magnitude()

                    # Create the barostat and add its force to the system.
                    self.addToConfig("\n# Add a barostat to run at constant pressure.")
                    self.addToConfig(f"barostat = MonteCarloBarostat({pressure}*bar, {temperature}*kelvin)")
                    if self._is_seeded:
                        self.addToConfig(f"barostat.setRandomNumberSeed({self._seed})")
                    self.addToConfig("system.addForce(barostat)")

            # Add backbone restraints. This uses the approach from:
            # https://github.com/openmm/openmm/issues/2262#issuecomment-464157489
            # Here zero-mass dummy atoms are bonded to the restrained atoms to avoid
            # issues with position rescaling during barostat updates.
            if self._protocol.isRestrained():
                # Get the list of backbone atom indices.
                restrained_atoms = self._system._getBackBoneAtoms()

                self.addToConfig("\n# Restrain the position of backbone atoms using zero-mass dummy atoms.")
                self.addToConfig("restraint = HarmonicBondForce()")
                self.addToConfig("restraint.setUsesPeriodicBoundaryConditions(True)")
                self.addToConfig("system.addForce(restraint)")
                self.addToConfig("nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]")
                self.addToConfig("dummy_indices = []")
                self.addToConfig("positions = inpcrd.positions")
                self.addToConfig(f"restrained_atoms = {restrained_atoms}")
                self.addToConfig("for i in restrained_atoms:")
                self.addToConfig("    j = system.addParticle(0)")
                self.addToConfig("    nonbonded.addParticle(0, 1, 0)")
                self.addToConfig("    nonbonded.addException(i, j, 0, 1, 0)")
                self.addToConfig("    restraint.addBond(i, j, 0*nanometers, 100*kilojoules_per_mole/nanometer**2)")
                self.addToConfig("    dummy_indices.append(j)")
                self.addToConfig("    positions.append(positions[i])")

            # Get the integration time step from the protocol.
            timestep = self._protocol.getTimeStep().picoseconds().magnitude()

            # Set the integrator.
            self.addToConfig( "\n# Define the integrator.")
            self.addToConfig(f"integrator = LangevinIntegrator({temperature}*kelvin,")
            self.addToConfig( "                                1/picosecond,")
            self.addToConfig(f"                                {timestep}*picoseconds)")
            if self._is_seeded:
                self.addToConfig(f"integrator.setRandomNumberSeed({self._seed})")

            # Add the platform information.
            self._add_config_platform()

            # Set up the simulation object.
            self.addToConfig("\n# Initialise and configure the simulation object.")
            self.addToConfig("simulation = Simulation(prmtop.topology,")
            self.addToConfig("                        system,")
            self.addToConfig("                        integrator,")
            self.addToConfig("                        platform,")
            self.addToConfig("                        properties)")
            if self._protocol.isRestrained():
                self.addToConfig("simulation.context.setPositions(positions)")
            else:
                self.addToConfig("simulation.context.setPositions(inpcrd.positions)")

            # Work out the number of integration steps.
            steps = _math.ceil(self._protocol.getRunTime() / self._protocol.getTimeStep())

            # Get the report and restart intervals.
            report_interval = self._protocol.getReportInterval()
            restart_interval = self._protocol.getRestartInterval()

            # Cap the intervals at the total number of steps.
            if report_interval > steps:
                report_interval = steps
            if restart_interval > steps:
                restart_interval = steps

            # Add the reporters.
            self.addToConfig("\n# Add reporters.")
            self._add_config_reporters(state_interval=report_interval, traj_interval=restart_interval)

            # Set initial velocities from temperature distribution.
            self.addToConfig("\n# Setting intial system velocities.")
            self.addToConfig(f"simulation.context.setVelocitiesToTemperature({temperature})")

            # Now run the simulation.
            self.addToConfig("\n# Run the simulation.")

            # Constant temperature equilibration.
            if self._protocol.isConstantTemp():
                self.addToConfig(f"simulation.step({steps})")

            # Heating / cooling cycle.
            else:
                # Adjust temperature every 100 cycles, assuming that there at
                # least that many cycles.
                if steps > 100:
                    # Work out the number of temperature cycles.
                    temp_cycles = _math.ceil(steps / 100)

                    # Work out the temperature change per cycle.
                    delta_temp = (self._protocol.getEndTemperature().kelvin().magnitude() -
                                  self._protocol.getStartTemperature().kelvin().magnitude()) / temp_cycles

                    self.addToConfig(f"start_temperature = {temperature}")
                    self.addToConfig(f"for x in range(0, {temp_cycles}):")
                    self.addToConfig(f"    temperature = {temperature} + x*{delta_temp}")
                    self.addToConfig(f"    integrator.setTemperature(temperature*kelvin)")
                    if is_const_pressure:
                        self.addToConfig(f"    barostat.setDefaultTemperature(temperature*kelvin)")
                    self.addToConfig( "    simulation.step(100)")
                else:
                    # Work out the temperature change per step.
                    delta_temp = (self._protocol.getEndTemperature().kelvin().magnitude() -
                                  self._protocol.getStartTemperature().kelvin().magnitude()) / steps

                    self.addToConfig(f"start_temperature = {temperature}")
                    self.addToConfig(f"for x in range(0, {temp_cycles}):")
                    self.addToConfig(f"    temperature = {temperature} + x*{delta_temp}")
                    self.addToConfig(f"    integrator.setTemperature(temperature*kelvin)")
                    if is_const_pressure:
                        self.addToConfig(f"    barostat.setDefaultTemperature(temperature*kelvin)")
                    self.addToConfig( "    simulation.step(1)")

        elif type(self._protocol) is _Protocol.Production:
            # Write the OpenMM import statements.
            self._add_config_imports()

            # Load the input files.
            self.addToConfig("\n# Load the topology and coordinate files.")
            self.addToConfig(f"prmtop = AmberPrmtopFile('{self._name}.prm7')")
            self.addToConfig(f"inpcrd = AmberInpcrdFile('{self._name}.rst7')")

            # Don't use a cut-off if this is a vacuum simulation or if box information
            # is missing.
            is_periodic = True
            self.addToConfig("\n# Initialise the molecular system.")
            if not has_box or not self._has_water:
                is_periodic = False
                self.addToConfig("system = prmtop.createSystem(nonbondedMethod=NoCutoff,")
            else:
                self.addToConfig("system = prmtop.createSystem(nonbondedMethod=PME,")
            self.addToConfig(    "                             nonbondedCutoff=1*nanometer,")
            self.addToConfig(    "                             constraints=HBonds)")

            # Get the starting temperature and system pressure.
            temperature = self._protocol.getTemperature().kelvin().magnitude()
            pressure = self._protocol.getPressure()

            # Add a Monte Carlo barostat if the simulation is at constant pressure.
            is_const_pressure = False
            if pressure is not None:
                # Cannot use a barostat with a non-periodic system.
                if not is_periodic:
                    _warnings.warn("Cannot use a barostat for a vacuum or non-periodic simulation")
                else:
                    is_const_pressure = True

                    # Convert to bar and get the magnitude.
                    pressure = pressure.bar().magnitude()

                    # Create the barostat and add its force to the system.
                    self.addToConfig("\n# Add a barostat to run at constant pressure.")
                    self.addToConfig(f"barostat = MonteCarloBarostat({pressure}*bar, {temperature}*kelvin)")
                    if self._is_seeded:
                        self.addToConfig(f"barostat.setRandomNumberSeed({self._seed})")
                    self.addToConfig("system.addForce(barostat)")

            # Get the integration time step from the protocol.
            timestep = self._protocol.getTimeStep().picoseconds().magnitude()

            # Set the integrator.
            self.addToConfig( "\n# Define the integrator.")
            self.addToConfig(f"integrator = LangevinIntegrator({temperature}*kelvin,")
            self.addToConfig( "                                1/picosecond,")
            self.addToConfig(f"                                {timestep}*picoseconds)")
            if self._is_seeded:
                self.addToConfig(f"integrator.setRandomNumberSeed({self._seed})")

            # Add the platform information.
            self._add_config_platform()

            # Set up the simulation object.
            self.addToConfig("\n# Initialise and configure the simulation object.")
            self.addToConfig("simulation = Simulation(prmtop.topology,")
            self.addToConfig("                        system,")
            self.addToConfig("                        integrator,")
            self.addToConfig("                        platform,")
            self.addToConfig("                        properties)")
            self.addToConfig("simulation.context.setPositions(inpcrd.positions)")

            # Work out the number of integration steps.
            steps = _math.ceil(self._protocol.getRunTime() / self._protocol.getTimeStep())

            # Get the report and restart intervals.
            report_interval = self._protocol.getReportInterval()
            restart_interval = self._protocol.getRestartInterval()

            # Cap the intervals at the total number of steps.
            if report_interval > steps:
                report_interval = steps
            if restart_interval > steps:
                restart_interval = steps

            # Add the reporters.
            self.addToConfig("\n# Add reporters.")
            self._add_config_reporters(state_interval=report_interval, traj_interval=restart_interval)

            # Set initial velocities from temperature distribution.
            self.addToConfig("\n# Setting intial system velocities.")
            self.addToConfig(f"simulation.context.setVelocitiesToTemperature({temperature})")

            # Now run the simulation.
            self.addToConfig("\n# Run the simulation.")
            self.addToConfig(f"simulation.step({steps})")

        elif type(self._protocol) is _Protocol.Metadynamics:
            colvar = self._protocol.getCollectiveVariable()
            if len(colvar) != 1 or (type(colvar[0]) is not _CollectiveVariable.Funnel):
                raise _IncompatibleError("We currently only support '%s' collective variables for '%s' protocols"
                        % (_CollectiveVariable.Funnel.__name__, self._protocol.__class__.__name__))

            # Create the path to the patched OpenMM metadynamics module.
            aux_file = "metadynamics.py"
            path = _os.path.dirname(_CollectiveVariable.__file__).replace("CollectiveVariable", "_aux") + "/" + aux_file

            # Copy the file into the working directory.
            _shutil.copyfile(path, self._work_dir + f"/{aux_file}")

            # The following OpenMM native implementation of the funnel metadynamics protocol
            # is adapted from funnel_maker.py by Dominykas Lukauskis.

            # Extract the only collective variable.
            colvar = colvar[0]

            # Write the OpenMM import statements.
            self._add_config_imports()
            self.addToConfig("from metadynamics import *")    # Use local patched metadynamics module.
            self.addToConfig("from glob import glob")
            self.addToConfig("import os")
            self.addToConfig("import shutil")

            # Load the input files.
            self.addToConfig("\n# Load the topology and coordinate files.")
            self.addToConfig(f"prmtop = AmberPrmtopFile('{self._name}.prm7')")
            self.addToConfig(f"inpcrd = AmberInpcrdFile('{self._name}.rst7')")

            # Don't use a cut-off if this is a vacuum simulation or if box information
            # is missing.
            is_periodic = True
            self.addToConfig("\n# Initialise the molecular system.")
            if not has_box or not self._has_water:
                is_periodic = False
                self.addToConfig("system = prmtop.createSystem(nonbondedMethod=NoCutoff,")
            else:
                self.addToConfig("system = prmtop.createSystem(nonbondedMethod=PME,")
            self.addToConfig(    "                             nonbondedCutoff=1*nanometer,")
            self.addToConfig(    "                             constraints=HBonds)")

            # Get the starting temperature and system pressure.
            temperature = self._protocol.getTemperature().kelvin().magnitude()
            pressure = self._protocol.getPressure()

            # Add a Monte Carlo barostat if the simulation is at constant pressure.
            is_const_pressure = False
            if pressure is not None:
                # Cannot use a barostat with a non-periodic system.
                if not is_periodic:
                    _warnings.warn("Cannot use a barostat for a vacuum or non-periodic simulation")
                else:
                    is_const_pressure = True

                    # Convert to bar and get the magnitude.
                    pressure = pressure.bar().magnitude()

                    # Create the barostat and add its force to the system.
                    self.addToConfig("\n# Add a barostat to run at constant pressure.")
                    self.addToConfig(f"barostat = MonteCarloBarostat({pressure}*bar, {temperature}*kelvin)")
                    if self._is_seeded:
                        self.addToConfig(f"barostat.setRandomNumberSeed({self._seed})")
                    self.addToConfig("system.addForce(barostat)")

            # Store the number of atoms in each molecule.
            atom_nums = []
            for mol in self._system:
                atom_nums.append(mol.nAtoms())

            # Sort the molecule indices by the number of atoms they contain.
            sorted_nums = sorted((value, idx) for idx, value in enumerate(atom_nums))

            # Set the ligand to the index of the second largest molecule.
            ligand = sorted_nums[-2][1]

            # Work out the start/end indices of the ligand within the system.
            idx_start = self._system.getIndex(self._system[ligand].getAtoms()[0])
            idx_end = self._system.getIndex(self._system[ligand+1].getAtoms()[0])

            # Add the funnel variables.
            p1_string = ",".join([str(x) for x in colvar.getAtoms0()])
            p2_string = ",".join([str(x) for x in colvar.getAtoms1()])
            self.addToConfig("\n# Set funnel variables.")
            self.addToConfig(f"p1 = [{p1_string}]")
            self.addToConfig(f"p2 = [{p2_string}]")
            self.addToConfig(f"lig = [x for x in range({idx_start}, {idx_end})]")

            sigma_proj = colvar.getHillWidth()[0].nanometers().magnitude()
            self.addToConfig( "\n# Create the bias variable for the funnel projection.")
            self.addToConfig( "projection = CustomCentroidBondForce(3, 'distance(g1,g2)*cos(angle(g1,g2,g3))')")
            self.addToConfig( "projection.addGroup(lig)")
            self.addToConfig( "projection.addGroup(p1)")
            self.addToConfig( "projection.addGroup(p2)")
            self.addToConfig( "projection.addBond([0,1,2])")
            self.addToConfig( "projection.setUsesPeriodicBoundaryConditions(True)")
            self.addToConfig(f"sigma_proj = {sigma_proj}")
            if colvar.getLowerBound() is None and colvar.getUpperBound() is None:
                # Sane defaults if no bounds are set.
                lower_wall = 0.5
                upper_wall = 4.5
            else:
                lower_wall = colvar.getLowerBound().getValue().nanometers().magnitude()
                upper_wall = colvar.getUpperBound().getValue().nanometers().magnitude()
            self.addToConfig(f"proj = BiasVariable(projection, {lower_wall-0.2}, {upper_wall+0.2}, {sigma_proj}, False, gridWidth=200)")

            sigma_ext = colvar.getHillWidth()[1].nanometers().magnitude()
            self.addToConfig("\n# Create the bias variable for the funnel extent.")
            self.addToConfig("extent = CustomCentroidBondForce(3, 'distance(g1,g2)*sin(angle(g1,g2,g3))')")
            self.addToConfig("extent.addGroup(lig)")
            self.addToConfig("extent.addGroup(p1)")
            self.addToConfig("extent.addGroup(p2)")
            self.addToConfig("extent.addBond([0,1,2])")
            self.addToConfig("extent.setUsesPeriodicBoundaryConditions(True)")
            extent_max = colvar.getWidth().nanometers().magnitude()  \
                       + colvar.getBuffer().nanometers().magnitude() \
                       + 0.2
            self.addToConfig(f"sigma_ext = {sigma_ext}")
            self.addToConfig(f"ext = BiasVariable(extent, 0.0, {extent_max}, {sigma_ext}, False, gridWidth=200)")

            # Add restraints.
            self.addToConfig("\n# Add restraints.")

            self.addToConfig("k1 = 10000*kilojoules_per_mole")
            self.addToConfig("k2 = 1000*kilojoules_per_mole")
            self.addToConfig(f"lower_wall = {lower_wall}*nanometer")
            self.addToConfig(f"upper_wall = {upper_wall}*nanometer")

            self.addToConfig("\n# Upper wall.")
            self.addToConfig("upper_wall_rest = CustomCentroidBondForce(3, '(k/2)*max(distance(g1,g2)*cos(angle(g1,g2,g3)) - upper_wall, 0)^2')")
            self.addToConfig("upper_wall_rest.addGroup(lig)")
            self.addToConfig("upper_wall_rest.addGroup(p1)")
            self.addToConfig("upper_wall_rest.addGroup(p2)")
            self.addToConfig("upper_wall_rest.addBond([0,1,2])")
            self.addToConfig("upper_wall_rest.addGlobalParameter('k', k1)")
            self.addToConfig("upper_wall_rest.addGlobalParameter('upper_wall', upper_wall)")
            self.addToConfig("upper_wall_rest.setUsesPeriodicBoundaryConditions(True)")
            self.addToConfig("system.addForce(upper_wall_rest)")

            self.addToConfig("\n# Sides of the funnel.")
            self.addToConfig(f"wall_width = {colvar.getWidth().nanometers().magnitude()}*nanometer")
            self.addToConfig(f"wall_buffer = {colvar.getBuffer().nanometers().magnitude()}*nanometer")
            self.addToConfig(f"beta_cent = {colvar.getSteepness()}")
            self.addToConfig(f"s_cent = {colvar.getInflection().nanometers().magnitude()}*nanometer")

            self.addToConfig("dist_restraint = CustomCentroidBondForce(3, '(k/2)*max(distance(g1,g2)*sin(angle(g1,g2,g3)) - (a/(1+exp(b*(distance(g1,g2)*cos(angle(g1,g2,g3))-c)))+d), 0)^2')")
            self.addToConfig("dist_restraint.addGroup(lig)")
            self.addToConfig("dist_restraint.addGroup(p1)")
            self.addToConfig("dist_restraint.addGroup(p2)")
            self.addToConfig("dist_restraint.addBond([0,1,2])")
            self.addToConfig("dist_restraint.addGlobalParameter('k', k2)")
            self.addToConfig("dist_restraint.addGlobalParameter('a', wall_width)")
            self.addToConfig("dist_restraint.addGlobalParameter('b', beta_cent)")
            self.addToConfig("dist_restraint.addGlobalParameter('c', s_cent)")
            self.addToConfig("dist_restraint.addGlobalParameter('d', wall_buffer)")
            self.addToConfig("dist_restraint.setUsesPeriodicBoundaryConditions(True)")
            self.addToConfig("system.addForce(dist_restraint)")

            self.addToConfig("\n# Lower wall.")
            self.addToConfig("lower_wall_rest = CustomCentroidBondForce(3, '(k/2)*min(distance(g1,g2)*cos(angle(g1,g2,g3)) - lower_wall, 0)^2')")
            self.addToConfig("lower_wall_rest.addGroup(lig)")
            self.addToConfig("lower_wall_rest.addGroup(p1)")
            self.addToConfig("lower_wall_rest.addGroup(p2)")
            self.addToConfig("lower_wall_rest.addBond([0,1,2])")
            self.addToConfig("lower_wall_rest.addGlobalParameter('k', k1)")
            self.addToConfig("lower_wall_rest.addGlobalParameter('lower_wall', lower_wall)")
            self.addToConfig("lower_wall_rest.setUsesPeriodicBoundaryConditions(True)")
            self.addToConfig("system.addForce(lower_wall_rest)")

            self.addToConfig("\n# Initialise the metadynamics object.")
            if self._protocol.getBiasFactor() is None:
                bias = 1.0
            else:
                bias = self._protocol.getBiasFactor()
            self.addToConfig(f"bias = {bias}")
            height = self._protocol.getHillHeight().kj_per_mol().magnitude()
            freq = self._protocol.getHillFrequency()

            # Work out the number of integration steps.
            steps = _math.ceil(self._protocol.getRunTime() / self._protocol.getTimeStep())

            # Get the report and restart intervals.
            report_interval = self._protocol.getReportInterval()
            restart_interval = self._protocol.getRestartInterval()

            # Cap the intervals at the total number of steps.
            if report_interval > steps:
                report_interval = steps
            if restart_interval > steps:
                restart_interval = steps

            # Work out the number of cycles.
            cycles = _math.ceil(steps / report_interval)

            self.addToConfig(f"meta = Metadynamics(system, [proj, ext], {temperature}*kelvin, {bias}, {height}*kilojoules_per_mole, {freq}, biasDir = '.', saveFrequency = {report_interval})")

            # Get the integration time step from the protocol.
            timestep = self._protocol.getTimeStep().picoseconds().magnitude()

            # Set the integrator.
            self.addToConfig( "\n# Define the integrator.")
            self.addToConfig(f"integrator = LangevinIntegrator({temperature}*kelvin,")
            self.addToConfig( "                                1/picosecond,")
            self.addToConfig(f"                                {timestep}*picoseconds)")
            if self._is_seeded:
                self.addToConfig(f"integrator.setRandomNumberSeed({self._seed})")

            # Add the platform information.
            self._add_config_platform()

            # Set up the simulation object.
            self.addToConfig("\n# Initialise and configure the simulation object.")
            self.addToConfig("simulation = Simulation(prmtop.topology,")
            self.addToConfig("                        system,")
            self.addToConfig("                        integrator,")
            self.addToConfig("                        platform,")
            self.addToConfig("                        properties)")
            self.addToConfig("simulation.context.setPositions(inpcrd.positions)")

            # Work out the number of integration.
            steps = _math.ceil(self._protocol.getRunTime() / self._protocol.getTimeStep())

            # Set initial velocities from temperature distribution.
            self.addToConfig("\n# Setting intial system velocities.")
            self.addToConfig(f"simulation.context.setVelocitiesToTemperature({temperature})")

            self.addToConfig("\n# Look for a restart file.")
            self.addToConfig(f"if os.path.isfile('{self._name}.chk'):")
            self.addToConfig(f"    simulation.loadCheckpoint('{self._name}.chk')")
            self.addToConfig(f"    shutil.copy('{self._name}.out','old_{self._name}.out')")
            self.addToConfig(f"    sim_log_file = [ line[:-2] for line in open('{self._name}.out').readlines()]")
            self.addToConfig( "    current_steps = int(sim_log_file[-1].split(',')[1])")
            self.addToConfig( "    steps -= current_steps")
            self.addToConfig( "    shutil.copy('COLVAR.npy','old_COLVAR.npy')")
            self.addToConfig(f"    shutil.copy('{self._name}.dcd','old_{self._name}.dcd')")

            # Add the reporters.
            self.addToConfig("\n# Add reporters.")
            self._add_config_reporters(state_interval=report_interval, traj_interval=restart_interval)
            self.addToConfig(f"simulation.reporters.append(CheckpointReporter('{self._name}.chk', {report_interval}))")

            # Create the HILLS file.
            self.addToConfig("\n# Create PLUMED compatible HILLS file.")
            self.addToConfig("file = open('HILLS','w')")
            self.addToConfig("file.write('#! FIELDS time pp.proj pp.ext sigma_pp.proj sigma_pp.ext height biasf\\n')")
            self.addToConfig("file.write('#! SET multivariate false\\n')")
            self.addToConfig("file.write('#! SET kerneltype gaussian\\n')")

            # Get the initial collective variables.
            self.addToConfig("\n# Initialise the collective variable array.")
            self.addToConfig("current_cvs = np.array(list(meta.getCollectiveVariables(simulation)) + [meta.getHillHeight(simulation)])")
            self.addToConfig("colvar_array = np.array([current_cvs])")

            # Run the metadynamics simulation.
            self.addToConfig("\n# Run the simulation.")
            self.addToConfig(f"steps = {steps}")
            self.addToConfig(f"cycles = {cycles}")
            self.addToConfig(f"steps_per_cycle = int({steps}/cycles)")
            self.addToConfig( "last_index = 0")
            self.addToConfig( "for x in range(0, cycles):")
            self.addToConfig( "    meta.step(simulation, steps_per_cycle)")
            self.addToConfig( "    current_cvs = np.array(list(meta.getCollectiveVariables(simulation)) + [meta.getHillHeight(simulation)])")
            self.addToConfig( "    colvar_array = np.append(colvar_array, [current_cvs], axis=0)")
            self.addToConfig( "    np.save('COLVAR.npy', colvar_array)")
            self.addToConfig( "    for index in range(last_index, np.shape(colvar_array)[0]):")
            self.addToConfig( "        line = colvar_array[index]")
            self.addToConfig( "        time = int(record_colvar_every.value_in_unit(picoseconds) * index)")
            self.addToConfig( "        write_line = f'{time:15} {line[0]:20.16f} {line[1]:20.16f}          {sigma_proj}           {sigma_ext} {line[2]:20.16f}            {bias}\\n'")
            self.addToConfig( "        file.write(write_line)")
            self.addToConfig( "    last_index = index")

        else:
            raise _IncompatibleError("Unsupported protocol: '%s'" % self._protocol.__class__.__name__)

    def getConfig(self):
        """Get the list of strings defining the OpenMM Python script.

           Returns
           -------

           config : [str]
               The list of configuration strings.
        """
        return super().getConfig()

    def setConfig(self, config):
        """Set the list of strings defining the OpenMM Python script.

           Parameters
           ----------

           config : str, [str]
               The list of configuration strings, or a path to a configuration
               file.
        """
        return super().setConfig()

    def addToConfig(self, config):
        """Add a string to the OpenMM Python script configuration.

           Parameters
           ----------

           config : str, [str]
               A configuration string, a list of configuration strings, or a
               path to a configuration file.
        """

        # Call the base class method.
        super().addToConfig(config)

    def resetConfig(self):
        """Reset the OpenMM Python script configuration."""
        self._generate_config()

    def setConfig(self, config):
        """Set the list of configuration file strings.

           Parameters
           ----------

           config : str, [str]
               The list of configuration strings, or a path to a configuration
               file.
        """

        # Call the base class method.
        super().setConfig(config)

    def start(self):
        """Start the OpenMM process.

           Returns
           -------

           process : :class:`Process.OpenMM <BioSimSpace.Process.OpenMM>`
               A handle to the OpenMM process.
        """

        # The process is currently queued.
        if self.isQueued():
            return

        # Process is already running.
        if self._process is not None:
            if self._process.isRunning():
                return

        # Clear any existing output.
        self._clear_output()

        # Run the process in the working directory.
        with _Utils.cd(self._work_dir):

            # Create the arguments string list.
            # The name of the Python script (config file) is the first argument.
            args = ["%s" % self._config_file]
            args.extend(self.getArgStringList())

            # Write the command-line process to a README.txt file.
            with open("README.txt", "w") as f:

                # Set the command-line string.
                self._command = "%s %s " % (self._exe, self._config_file) + self.getArgString()

                # Write the command to file.
                f.write("# OpenMM was run with the following command:\n")
                f.write("%s\n" % self._command)

            # Start the timer.
            self._timer = _timeit.default_timer()

            # Start the simulation.
            self._process = _SireBase.Process.run(self._exe, args,
                "%s.out" % self._name, "%s.err" % self._name)

        return self

    def getSystem(self, block="AUTO"):
        """Get the latest molecular system.

           Parameters
           ----------

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The latest molecular system.
        """

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        # Warn the user if the process has exited with an error.
        if self.isError():
            _warnings.warn("The process exited with an error!")

        # Try to get the most recent trajectory frame.
        try:
            # Handle minimisation protocols separately.
            if type(self._protocol) is _Protocol.Minimisation:
                traj = self.getTrajectory()

                # If there is no trajectory, simply return None.
                if traj is None:
                    return None

                # Get the last frame.
                new_system = traj.getFrames(-1)[0]

            else:
                # Work out the total number of trajectory frames.
                num_frames = int((self._protocol.getRunTime() / self._protocol.getTimeStep())
                    / self._protocol.getRestartInterval())

                # Work out the fraction of the simulation that has been completed.
                frac_complete = self._protocol.getRunTime() / self.getTime()

                # Work out the trajectory frame index, rounding down.
                index = int(frac_complete * num_frames)

                # Get the most recent frame.
                new_system = self.getFrame(index)

            # Copy the new coordinates back into the original system.
            old_system = self._system.copy()
            old_system._updateCoordinates(new_system,
                                          self._property_map,
                                          self._property_map)

            # Update the box information in the original system.
            if "space" in new_system._sire_object.propertyKeys():
                box = new_system._sire_object.property("space")
                old_system._sire_object.setProperty(self._property_map.get("space", "space"), box)

            return old_system

        except:
            return None

    def getCurrentSystem(self):
        """Get the latest molecular system.

           Returns
           -------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The latest molecular system.
        """
        return self.getSystem(block=False)

    def getTrajectory(self, block="AUTO"):
        """Return a trajectory object.

           Parameters
           ----------

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           trajectory : :class:`System <BioSimSpace.Trajectory.Trajectory>`
               The latest trajectory object.
        """

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        # Warn the user if the process has exited with an error.
        if self.isError():
            _warnings.warn("The process exited with an error!")

        if not _os.path.isfile(self._traj_file):
            return None
        else:
            return _Trajectory.Trajectory(process=self)

    def getFrame(self, index):
        """Return a specific trajectory frame.

           Parameters
           ----------

           index : int
               The index of the frame.

          Returns
          -------

          frame : :class:`System <BioSimSpace._SireWrappers.System>`
              The System object of the corresponding frame.
        """

        if type(index) is not int:
            raise TypeError("'index' must be of type 'int'")

        max_index = int((self._protocol.getRunTime() / self._protocol.getTimeStep())
                  / self._protocol.getRestartInterval())

        if index < 0 or index > max_index:
            raise ValueError(f"'index' must be in range [0, {max_index}].")

        try:
            new_system =  _Trajectory.getFrame(self._traj_file,
                                               self._top_file,
                                               index)

            # Copy the new coordinates back into the original system.
            old_system = self._system.copy()
            old_system._updateCoordinates(new_system,
                                          self._property_map,
                                          self._property_map)

            # Update the box information in the original system.
            if "space" in new_system._sire_object.propertyKeys():
                box = new_system._sire_object.property("space")
                old_system._sire_object.setProperty(self._property_map.get("space", "space"), box)

            return old_system

        except:
            return None

    def getRecord(self, record, time_series=False, unit=None, block="AUTO"):
        """Get a record from the stdout dictionary.

           Parameters
           ----------

           record : str
               The record key.

           time_series : bool
               Whether to return a list of time series records.

           unit : :class:`Unit <BioSimSpace.Units>`
               The unit to convert the record to.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           record : :class:`Type <BioSimSpace.Types>`
               The matching record.
        """

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        # Warn the user if the process has exited with an error.
        if self.isError():
            _warnings.warn("The process exited with an error!")

        self._update_stdout_dict()
        return self._get_stdout_record(record, time_series, unit)

    def getCurrentRecord(self, record, time_series=False, unit=None):
        """Get a current record from the stdout dictionary.

           Parameters
           ----------

           record : str
               The record key.

           time_series : bool
               Whether to return a list of time series records.

           unit : :class:`Unit <BioSimSpace.Units>`
               The unit to convert the record to.

           Returns
           -------

           record : :class:`Type <BioSimSpace.Types>`
               The matching record.
        """
        # Warn the user if the process has exited with an error.
        if self.isError():
            _warnings.warn("The process exited with an error!")

        self._update_stdout_dict()
        return self._get_stdout_record(record, time_series, unit)

    def getRecords(self, block="AUTO"):
        """Return the dictionary of stdout time-series records.

           Parameters
           ----------

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           records : :class:`MultiDict <BioSimSpace.Process._process._MultiDict>`
              The dictionary of time-series records.
        """
        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        # Warn the user if the process has exited with an error.
        if self.isError():
            _warnings.warn("The process exited with an error!")

        return self._stdout_dict.copy()

    def getCurrentRecords(self):
        """Return the current dictionary of stdout time-series records.

           Parameters
           ----------

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           records : :class:`MultiDict <BioSimSpace.Process._process._MultiDict>`
              The dictionary of time-series records.
        """
        return getRecords(block=False)

    def getTime(self, time_series=False, block="AUTO"):
        """Get the simulation time.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           time : :class:`Time <BioSimSpace.Types.Time>`
               The current simulation time in nanoseconds.
        """
        return self.getRecord("TIME(PS)", time_series, _Units.Time.picosecond, block)

    def getCurrentTime(self, time_series=False):
        """Get the current simulation time.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           time : :class:`Time <BioSimSpace.Types.Time>`
               The current simulation time in nanoseconds.
        """
        return self.getTime(time_series, block=False)

    def getStep(self, time_series=False, block="AUTO"):
        """Get the number of integration steps.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           step : int
               The current number of integration steps.
        """
        return self.getRecord("STEP", time_series, None, block)

    def getCurrentStep(self, time_series=False):
        """Get the current number of integration steps.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           step : int
               The current number of integration steps.
        """
        return self.getStep(time_series, block=False)

    def getPotentialEnergy(self, time_series=False, block="AUTO"):
        """Get the potential energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The potential energy.
        """
        return self.getRecord("POTENTIALENERGY(KJ/MOLE)", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentPotentialEnergy(self, time_series=False):
        """Get the current potential energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The potential energy.
        """
        return self.getPotentialEnergy(time_series, block=False)

    def getKineticEnergy(self, time_series=False, block="AUTO"):
        """Get the kinetic energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The kinetic energy.
        """
        return self.getRecord("KINETICENERGY(KJ/MOLE)", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentKineticEnergy(self, time_series=False):
        """Get the current kinetic energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The kinetic energy.
        """
        return self.getKineticEnergy(time_series, block=False)

    def getTotalEnergy(self, time_series=False, block="AUTO"):
        """Get the total energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The total energy.
        """
        return self.getRecord("TOTALENERGY(KJ/MOLE)", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentTotalEnergy(self, time_series=False):
        """Get the current total energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
               The total energy.
        """
        return self.getTotalEnergy(time_series, block=False)

    def getTemperature(self, time_series=False, block="AUTO"):
        """Get the temperature.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
               The temperature.
        """
        return self.getRecord("TEMPERATURE(K)", time_series, _Units.Temperature.kelvin, block)

    def getCurrentTemperature(self, time_series=False):
        """Get the current temperature.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
               The current temperature.
        """
        return self.getTemperature(time_series, block=False)

    def getVolume(self, time_series=False, block="AUTO"):
        """Get the volume.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           volume : :class:`Volume <BioSimSpace.Types.Volume>`
              The volume.
        """
        return self.getRecord("BOXVOLUME(NM^3)", time_series, _Units.Volume.nanometer3, block)

    def getCurrentVolume(self, time_series=False):
        """Get the current volume.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           volume : :class:`Volume <BioSimSpace.Types.Volume>`
              The volume.
        """
        return self.getVolume(time_series, block=False)

    def stdout(self, n=10):
        """Print the last n lines of the stdout buffer.

           Parameters
           ----------

           n : int
               The number of lines to print.
        """

        # Note that thermodynamic records, e.g. energy, pressure, temperture,
        # are redirected to a log file.

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

    def _add_config_imports(self):
        """Helper function to write the header (import statements) to the
           OpenMM Python script (config file).
        """
        self.addToConfig("from simtk.openmm.app import *")
        self.addToConfig("from simtk.openmm import *")
        self.addToConfig("from simtk.unit import *")

    def _add_config_platform(self):
        """Helper function to add platform information to the OpenMM
           Python script.
        """
        # Set the simulation platform.
        self.addToConfig("\n# Set the simulation platform.")
        self.addToConfig(f"platform = Platform.getPlatformByName('{self._platform}')")
        self.addToConfig("properties = {}")

        if self._platform == "CPU":
            self.addToConfig("properties = {}")
        elif self._platform == "CUDA":
            cuda_devices = _os.environ.get("CUDA_VISIBLE_DEVICES")
            if cuda_devices is None:
                raise EnvironmentError("'CUDA' platform selected but 'CUDA_VISIBLE_DEVICES' "
                                        "environment variable is unset.")
            else:
                self.addToConfig(f"properties = {{'CudaDeviceIndex': '{cuda_devices}'}}")
        elif self._platform == "OPENCL":
            opencl_devices = _os.environ.get("OPENCL_VISIBLE_DEVICES")
            if opencl_devices is None:
                raise EnvironmentError("'OpenCL' platform selected but 'OPENCL_VISIBLE_DEVICES' "
                                        "environment variable is unset.")
            else:
                self.addToConfig(f"properties = {{'OpenCLDeviceIndex': '{opencl_devices}'}}")

    def _add_config_monkey_patches(self):
        """Helper function to write any monkey-patches to the OpenMM Python
           script (config file).
        """
        # We monkey-patch the OpenMM DCDFile.writeModel method to avoid writing the
        # positions of any dummy atoms that are used as restraints to trajectory files.
        # This avoids the need to delete the dummies from the molecular system on read,
        # allowing us to make use of System._updateCoordinates which requires that the
        # number of atoms are consistent between systems. (Deleting the dummies from the
        # system is slower than not writing them in the first place.)
        self.addToConfig("\n# Monkey-patch the DCD.writeModel method to avoid writing dummy-atom positions.")

        # Store the original writeModel method.
        self.addToConfig("writeModel = DCDFile.writeModel")

        # Create a monkey-patch where we slice the positions list to match the
        # number of atoms in the topology, then pass this through to the original
        # writeModel method.
        self.addToConfig("def writeModelPatched(self, positions, unitCellDimensions=None, periodicBoxVectors=None):")
        self.addToConfig("    positions = positions[:len(list(self._topology.atoms()))]")
        self.addToConfig("    writeModel(self,")
        self.addToConfig("               positions,")
        self.addToConfig("               unitCellDimensions=unitCellDimensions,")
        self.addToConfig("               periodicBoxVectors=periodicBoxVectors)")

        # Replace the writeModel method with the monkey-patch.
        self.addToConfig("DCDFile.writeModel = writeModelPatched")

    def _add_config_reporters(self, state_interval=100, traj_interval=500):
        """Helper function to write the reporter (output statements) section
           to the OpenMM Python script (config file).

           Parameters
           ----------

           state_interval : int
               The frequency at which to write state information in
               integration steps.

           traj_interval : int
               The frequency at which to write trajectory frames in
               integration steps.
        """
        if type(state_interval) is not int:
            raise TypeError("'state_interval' must be of type 'int'.")
        if state_interval <= 0:
            raise ValueError("'state_interval' must be a positive integer.")
        if type(traj_interval) is not int:
            raise TypeError("'traj_interval' must be of type 'int'.")
        if traj_interval <= 0:
            raise ValueError("'traj_interval' must be a positive integer.")

		# Append to a trajectory file every 500 steps.
        self.addToConfig(f"simulation.reporters.append(DCDReporter('{self._name}.dcd', {traj_interval}))")

		# Write state information to file every 100 steps.
        self.addToConfig(f"simulation.reporters.append(StateDataReporter('{self._name}.log',")
        self.addToConfig(f"                                              {state_interval},")
        self.addToConfig( "                                              step=True,")
        self.addToConfig( "                                              time=True,")
        self.addToConfig( "                                              potentialEnergy=True,")
        self.addToConfig( "                                              kineticEnergy=True,")
        self.addToConfig( "                                              totalEnergy=True,")
        self.addToConfig( "                                              volume=True,")
        self.addToConfig( "                                              temperature=True,")
        self.addToConfig( "                                              totalSteps=True,")
        self.addToConfig( "                                              separator=' '))")

    def _update_stdout_dict(self):
        """Update the dictonary of thermodynamic records."""

        # Exit if log file hasn't been created.
        if not _os.path.isfile(self._log_file):
            return

        # A list of the new record lines.
        lines = []

        # Append any new lines.
        for line in _pygtail.Pygtail(self._log_file):
            lines.append(line)

        # Append any new records to the stdout dictionary.
        for line in lines:

            # Strip leading/trailing whitespace.
            line = line.strip()

            # This is the header record.
            if line[0] == "#":
                # Work out what records are in the file and the separator
                # that is used. While we use a standard format, this makes
                # sure that we can still parse the log file if the user
                # happens to have changed the formatting.

                # A tally counter for the number of quotes that we've seen
                # in the line so far.
                num_quotes = 0

                # Initalise the separator.
                sep = ""

                # Loop over the characters in the line.
                for c in line:
                    # Increment the number of quotes.
                    if c == '"':
                        num_quotes += 1
                    # This is the second quote we've seen, start adding
                    # characters to the separator.
                    if num_quotes == 2:
                        sep += c
                    # Break when we've reached the next quote.
                    elif num_quotes == 3:
                        break

                # The separator includes a leading " character, so delete it.
                self._record_separator = sep[1:]

                # Now split the line on the separator to work out the records.
                # We ignore the first character since it is a comment.
                # Here we use the full separator, i.e. including the both quotes,
                # so that we can correctly split record names with spaces in them.
                records = line[1:].split(sep + '"')

                # Store the number of records.
                num_records = len(records)

                # Map each record string to its position in the array (column order).
                for idx, record in enumerate(records):
                    # Strip the extra quotes from the record.
                    if idx == 0:
                        record = record[1:]
                    elif idx == num_records - 1:
                        record = record[:-1]

                    # Map the index to the record. Store records in upper
                    # case without whitespace to help catch typos from
                    # the user.
                    self._record_mapping[idx] = record.upper().replace(" ", "")

            # Extract the records and add them to the dictionary.
            else:
                # Split the line on the separator.
                records = line.split(self._record_separator)

                # Add each record to the appropriate key in the MultiDict.
                for idx, record in enumerate(records):
                    # Get the key for this record.
                    key = self._record_mapping[idx]

                    # Update the record dictionary for this key.
                    self._stdout_dict[key] = record

    def _get_stdout_record(self, key, time_series=False, unit=None):
        """Helper function to get a stdout record from the dictionary.

           Parameters
           ----------

           key : str
               The record key.

           time_series : bool
               Whether to return a time series of records.

           unit : BioSimSpace.Types._type.Type
               The unit to convert the record to.

           Returns
           -------

           record :
               The matching stdout record.
        """

        # No data!
        if len(self._stdout_dict) is 0:
            return None

        if type(time_series) is not bool:
            _warnings.warn("Non-boolean time-series flag. Defaulting to False!")
            time_series = False

        # Valdate the unit.
        if unit is not None:
            if not isinstance(unit, _Type):
                raise TypeError("'unit' must be of type 'BioSimSpace.Types'")

        # Return the list of dictionary values.
        if time_series:
            try:
                if unit is None:
                    return [float(x) for x in self._stdout_dict[key]]
                else:
                    return [(float(x) * unit)._default_unit() for x in self._stdout_dict[key]]

            except KeyError:
                return None

        # Return the most recent dictionary value.
        else:
            try:
                if unit is None:
                    return float(self._stdout_dict[key][-1])
                else:
                    return (float(self._stdout_dict[key][-1]) * unit)._default_unit()

            except KeyError:
                return None
