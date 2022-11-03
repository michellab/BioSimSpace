######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2022
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
Functionality for running simulations using AMBER.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Amber"]

from .._Utils import _try_import, _have_imported

_watchdog = _try_import("watchdog")

if _have_imported(_watchdog):
    from watchdog.events import PatternMatchingEventHandler as _PatternMatchingEventHandler
    from watchdog.observers import Observer as _Observer
else:
    _PatternMatchingEventHandler = _watchdog
    _Observer = _watchdog

import math as _math
import os as _os
import re as _re
import time as _time
import shutil as _shutil
import timeit as _timeit
import warnings as _warnings

from sire.legacy import Base as _SireBase
from sire.legacy import IO as _SireIO
from sire.legacy import Mol as _SireMol

from .. import _amber_home, _isVerbose
from .._Exceptions import IncompatibleError as _IncompatibleError
from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from .._SireWrappers import System as _System
from ..Types._type import Type as _Type

from .. import Protocol as _Protocol
from .. import Trajectory as _Trajectory
from .. import Units as _Units
from .. import _Utils

from . import _process

from ._plumed import Plumed as _Plumed

class _Watcher:
    """A class to watch for changes to the AMBER energy info file. An event handler
       is used trigger updates to the energy dictionary each time the file is modified.
    """
    def __init__(self, proc):
        """Constructor.

           Parameters
           ----------

           proc : :class:`Process.Amber <BioSimSpace.Process.Amber>`
               The Amber Process object.
        """

        self._process = proc
        self._observer = _Observer()

    def start(self):
        """Start the file watcher."""

        # Setup the event handler and observer.
        event_handler = _Handler(self._process)
        self._observer.schedule(event_handler, self._process._work_dir)
        self._observer.daemon = True
        self._observer.start()

if _have_imported(_watchdog):
    class _Handler(_PatternMatchingEventHandler):
        """An event handler to trigger updates to the energy dictionary each time
           the log file is changed.
        """

        def __init__(self, proc):
            """Constructor.

               Parameters
               ----------

               proc : :class:`Process.Amber <BioSimSpace.Process.Amber>`
                   The Amber Process object.
            """
            self._process = proc

            super(_Handler, self).__init__(case_sensitive=False,
                                           ignore_directories=True,
                                           ignore_patterns=[],
                                           patterns=["*.nrg"])

        def on_any_event(self, event):
            """Update the dictionary when the file is modified.

               Parameters
               ----------

               event : str
                   The file system event.
            """

            # N.B.
            #
            # Since the watchdog package is cross-platform it doesn't support
            # detection of "close-write" operations, so multiple "modified" events
            # can be triggered while the log file is being written. As such, we
            # check whether the file has been updated by seeing if the NSTEP record
            # is different to the most recent entry in the dictionary.
            # So far, no issues have been found with processing partially written
            # files, i.e. duplicate or missing records.

            if event.event_type == "modified":
                # If this is the first time the file has been modified since the
                # process started, then wipe the dictionary and flag that the file
                # is now being watched.
                if not self._process._is_watching:
                    self._process._stdout_dict = _process._MultiDict()
                    self._process._is_watching = True

                # Make sure the file exists.
                if _os.path.isfile(self._process._nrg_file):
                    # Now update the dictionary with any new records.
                    self._process._update_energy_dict()
else:
    _Handler = _watchdog

class Amber(_process.Process):
    """A class for running simulations using AMBER."""

    def __init__(self, system, protocol, exe=None, name="amber",
            work_dir=None, seed=None, property_map={}):
        """Constructor.

           Parameters
           ----------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The molecular system.

           protocol : :class:`Protocol <BioSimSpace.Protocol>`
               The protocol for the AMBER process.

           exe : str
               The full path to the AMBER executable.

           name : str
               The name of the process.

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
        self._package_name = "AMBER"

        # This process can generate trajectory data.
        self._has_trajectory = True

        # If the path to the executable wasn't specified, then search
        # for it in $PATH. For now, we'll just search for 'sander', which
        # is available free as part of AmberTools. In future, we will
        # look for all possible executables in order of preference: pmemd.cuda,
        # pmemd, sander, etc., as well as their variants, e.g. pmemd.MPI.
        if exe is None:
            # Search AMBERHOME, if set.
            if _amber_home is not None:
                exe = "%s/bin/sander" % _amber_home
                if _os.path.isfile(exe):
                    self._exe = exe
                else:
                    raise _MissingSoftwareError("'BioSimSpace.Process.Amber' is not supported. "
                                                "Please install AMBER (http://ambermd.org).")
        else:
            # Make sure executable exists.
            if _os.path.isfile(exe):
                self._exe = exe
            else:
                raise IOError("AMBER executable doesn't exist: '%s'" % exe)

        # Initialise the energy dictionary and header.
        self._stdout_dict = _process._MultiDict()

        # Create the name of the energy output file and wipe the
        # contents of any existing file.
        self._nrg_file = "%s/%s.nrg" % (self._work_dir, name)
        open(self._nrg_file, "w").close()

        # Initialise the energy watcher.
        self._watcher = None
        self._is_watching = False

        # The names of the input files.
        self._rst_file = "%s/%s.rst7" % (self._work_dir, name)
        self._top_file = "%s/%s.prm7" % (self._work_dir, name)

        # The name of the trajectory file.
        self._traj_file = "%s/%s.nc" % (self._work_dir, name)

        # Set the path for the AMBER configuration file.
        self._config_file = "%s/%s.cfg" % (self._work_dir, name)

        # Create the list of input files.
        self._input_files = [self._config_file, self._rst_file, self._top_file]

        # Now set up the working directory for the process.
        self._setup()

    def _setup(self):
        """Setup the input files and working directory ready for simulation."""

        # Create the input files...

        # Create a copy of the system.
        system = self._system.copy()

        # Convert the water model topology so that it matches the AMBER naming convention.
        system._set_water_topology("AMBER", self._property_map)

        # Check for perturbable molecules and convert to the chosen end state.
        system = self._checkPerturbable(system)

        # RST file (coordinates).
        try:
            rst = _SireIO.AmberRst7(system._sire_object, self._property_map)
            rst.writeToFile(self._rst_file)
        except Exception as e:
            msg = "Failed to write system to 'RST7' format."
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # PRM file (topology).
        try:
            prm = _SireIO.AmberPrm(system._sire_object, self._property_map)
            prm.writeToFile(self._top_file)

        except Exception as e:
            msg = "Failed to write system to 'PRM7' format."
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # Generate the AMBER configuration file.
        # Skip if the user has passed a custom config.
        if isinstance(self._protocol, _Protocol.Custom):
            self.setConfig(self._protocol.getConfig())
        else:
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

        prop = self._property_map.get("space", "space")

        # Check whether the system contains periodic box information.
        # For now, we'll not attempt to generate a box if the system property
        # is missing. If no box is present, we'll assume a non-periodic simulation.
        if prop in self._system._sire_object.propertyKeys():
            has_box = True
        else:
            _warnings.warn("No simulation box found. Assuming gas phase simulation.")
            has_box = False

        # Work out whether we're generating a config for PMEMD.
        if "pmemd" in self._exe.lower():
            is_pmemd = True
        else:
            is_pmemd = False

        # While the configuration parameters below share a lot of overlap,
        # we choose the keep them separate so that the user can modify options
        # for a given protocol in a single place.

        # Add configuration variables for a minimisation simulation.
        if isinstance(self._protocol, _Protocol.Minimisation):

            # Work out the number of steepest descent cycles.
            # This is 1000 or 10% of the number of steps, whichever is larger.
            num_steps = self._protocol.getSteps()
            if num_steps <= 1000:
                num_steep = num_steps
            else:
                num_steep = _math.ceil(num_steps/10)
                if num_steep < 1000:
                    num_steep = 1000

            self.addToConfig("Minimisation")
            self.addToConfig(" &cntrl")
            self.addToConfig("  imin=1,")                   # Minimisation simulation.
            self.addToConfig("  ntx=1,")                    # Only read coordinates from file.
            self.addToConfig("  ntxo=1,")                   # Output coordinates in ASCII.
            self.addToConfig("  ntpr=100,")                 # Output energies every 100 steps.
            self.addToConfig("  irest=0,")                  # Don't restart.
            self.addToConfig("  maxcyc=%d," % num_steps)    # Set the number of steps.
            self.addToConfig("  ncyc=%d," % num_steep)      # Set the number of steepest descent steps.
            if not has_box or not self._has_water:
                self.addToConfig("  ntb=0,")                # No periodic box.
                self.addToConfig("  cut=999.,")             # Non-bonded cut-off.
                if is_pmemd:
                    self.addToConfig("  igb=6,")            # Use vacuum generalised Born model.
            else:
                self.addToConfig("  cut=8.0,")              # Non-bonded cut-off.
            self.addToConfig(" /")

        # Add configuration variables for an equilibration simulation.
        elif isinstance(self._protocol, _Protocol.Equilibration):

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

            # Set the random number seed.
            if self._is_seeded:
                seed = self._seed
            else:
                seed = -1

            # Convert the timestep to picoseconds.
            timestep = self._protocol.getTimeStep().picoseconds().value()

            self.addToConfig("Equilibration.")
            self.addToConfig(" &cntrl")
            self.addToConfig("  ig=%d," % seed)                 # Random number seed.
            self.addToConfig("  ntx=1,")                        # Only read coordinates from file.
            self.addToConfig("  ntxo=1,")                       # Output coordinates in ASCII.
            self.addToConfig("  ntpr=%d," % report_interval)    # Interval between reporting energies.
            self.addToConfig("  ntwr=%d," % restart_interval)   # Interval between saving restart files.
            self.addToConfig("  ntwx=%d," % restart_interval)   # Trajectory sampling frequency.
            self.addToConfig("  irest=0,")                      # Don't restart.
            self.addToConfig("  dt=%.3f," % timestep)           # Time step.
            self.addToConfig("  nstlim=%d," % steps)            # Number of integration steps.
            self.addToConfig("  ntc=2,")                        # Enable SHAKE.
            self.addToConfig("  ntf=2,")                        # Don't calculate forces for constrained bonds.
            self.addToConfig("  ntt=3,")                        # Langevin dynamics.
            self.addToConfig("  gamma_ln=2,")                   # Collision frequency (ps).
            if not has_box or not self._has_water:
                self.addToConfig("  ntb=0,")                    # No periodic box.
                self.addToConfig("  cut=999.,")                 # Non-bonded cut-off.
                if is_pmemd:
                    self.addToConfig("  igb=6,")                # Use vacuum generalised Born model.
            else:
                self.addToConfig("  cut=8.0,")                  # Non-bonded cut-off.

            # Constant pressure control.
            if self._protocol.getPressure() is not None:
                # Don't use barostat for vacuum simulations.
                if has_box and self._has_water:
                    self.addToConfig("  ntp=1,")                # Isotropic pressure scaling.
                    self.addToConfig("  pres0=%.5f,"            # Pressure in bar.
                        % self._protocol.getPressure().bar().value())
                else:
                    _warnings.warn("Cannot use a barostat for a vacuum or non-periodic simulation")

            # Restrain the backbone.
            restraint = self._protocol.getRestraint()
            if restraint is not None:
                # Get the indices of the atoms that are restrained.
                if isinstance(restraint, str):
                    atom_idxs = self._system.getRestraintAtoms(restraint)
                else:
                    atom_idxs = restraint

                # Don't add restraints if there are no atoms to restrain.
                if len(atom_idxs) > 0:
                    # Generate the restraint mask based on atom indices.
                    restraint_mask = self._create_restraint_mask(atom_idxs)

                    # The restraintmask cannot be more than 256 characters.
                    if len(restraint_mask) > 256:

                        # AMBER has a limit on the length of the restraintmask
                        # so it's easy to overflow if we are matching by index
                        # on a large protein. As such, handle "backbone" and
                        # "heavy" restraints using a non-interoperable name mask.
                        if isinstance(restraint, str):
                            if restraint == "backbone":
                                restraint_mask = "@CA,C,O,N"
                            elif restraint == "heavy":
                                restraint_mask = "!:WAT & !@H"
                            elif restraint == "all":
                                restraint_mask = "!:WAT"

                        # We can't do anything about a custom restraint, since we don't
                        # know anything about the atoms.
                        else:
                            raise ValueError("AMBER atom 'restraintmask' exceeds 256 character limit!")

                    # Get the force constant value. The default is the same
                    # units as AMBER, i.e. kcal_per_mol/angstrom**2
                    force_constant = self._protocol.getForceConstant().value()

                    self.addToConfig( "  ntr=1,")
                    self.addToConfig(f"  restraint_wt={force_constant},")
                    self.addToConfig(f"  restraintmask='{restraint_mask}',")

            # Heating/cooling protocol.
            if not self._protocol.isConstantTemp():
                self.addToConfig("  tempi=%.2f," % self._protocol.getStartTemperature().kelvin().value())
                self.addToConfig("  temp0=%.2f," % self._protocol.getEndTemperature().kelvin().value())
                self.addToConfig("  nmropt=1,")
                self.addToConfig(" /")
                self.addToConfig("&wt TYPE='TEMP0', istep1=0, istep2=%d, value1=%.2f, value2=%.2f /"
                    % (steps, self._protocol.getStartTemperature().kelvin().value(),
                       self._protocol.getEndTemperature().kelvin().value()))
                self.addToConfig("&wt TYPE='END' /")

            # Constant temperature equilibration.
            else:
                self.addToConfig("  temp0=%.2f," % self._protocol.getStartTemperature().kelvin().value())
                self.addToConfig(" /")

        # Add configuration variables for a production simulation.
        elif isinstance(self._protocol, _Protocol.Production):

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

            # Set the random number seed.
            if self._seed is None:
                seed = -1
            else:
                seed = self._seed

            # Convert the timestep to picoseconds.
            timestep = self._protocol.getTimeStep().picoseconds().value()

            self.addToConfig("Production.")
            self.addToConfig(" &cntrl")
            self.addToConfig("  ig=%d," % seed)                 # Random number seed.
            if self._protocol.isRestart():
                self.addToConfig("  ntx=5,")                    # Read coordinates and velocities.
            else:
                self.addToConfig("  ntx=1,")                    # Only read coordinates.
            self.addToConfig("  ntxo=1,")                       # Output coordinates in ASCII.
            self.addToConfig("  ntpr=%d," % report_interval)    # Interval between reporting energies.
            self.addToConfig("  ntwr=%d," % restart_interval)   # Interval between saving restart files.
            self.addToConfig("  ntwx=%d," % restart_interval)   # Trajectory sampling frequency.
            if self._protocol.isRestart():
                self.addToConfig("  irest=1,")                  # Restart using previous velocities.
            else:
                self.addToConfig("  irest=0,")                  # Don't restart.
            self.addToConfig("  dt=%.3f," % timestep)           # Time step.
            self.addToConfig("  nstlim=%d," % steps)            # Number of integration steps.
            self.addToConfig("  ntc=2,")                        # Enable SHAKE.
            self.addToConfig("  ntf=2,")                        # Don't calculate forces for constrained bonds.
            self.addToConfig("  ntt=3,")                        # Langevin dynamics.
            self.addToConfig("  gamma_ln=2,")                   # Collision frequency (ps).
            if not has_box or not self._has_water:
                self.addToConfig("  ntb=0,")                    # No periodic box.
                self.addToConfig("  cut=999.,")                 # Non-bonded cut-off.
                if is_pmemd:
                    self.addToConfig("  igb=6,")                # Use vacuum generalised Born model.
            else:
                self.addToConfig("  cut=8.0,")                  # Non-bonded cut-off.
            if not self._protocol.isRestart():
                self.addToConfig("  tempi=%.2f,"                # Initial temperature.
                    % self._protocol.getTemperature().kelvin().value())
            self.addToConfig("  temp0=%.2f,"                    # Target temperature.
                % self._protocol.getTemperature().kelvin().value())

            # Constant pressure control.
            if self._protocol.getPressure() is not None:
                # Don't use barostat for vacuum simulations.
                if has_box and self._has_water:
                    self.addToConfig("  ntp=1,")                # Isotropic pressure scaling.
                    self.addToConfig("  pres0=%.5f,"            # Pressure in bar.
                        % self._protocol.getPressure().bar().value())
                else:
                    _warnings.warn("Cannot use a barostat for a vacuum or non-periodic simulation")

            self.addToConfig(" /")

        # Add configuration variables for a metadynamics simulation.
        elif isinstance(self._protocol, _Protocol.Metadynamics):

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

            # Set the random number seed.
            if self._seed is None:
                seed = -1
            else:
                seed = self._seed

            # Convert the timestep to picoseconds.
            timestep = self._protocol.getTimeStep().picoseconds().value()

            self.addToConfig("Production.")
            self.addToConfig(" &cntrl")
            self.addToConfig("  ig=%d," % seed)                 # Random number seed.
            self.addToConfig("  ntx=1,")                        # Only read coordinates.
            self.addToConfig("  ntxo=1,")                       # Output coordinates in ASCII.
            self.addToConfig("  ntpr=%d," % report_interval)    # Interval between reporting energies.
            self.addToConfig("  ntwr=%d," % restart_interval)   # Interval between saving restart files.
            self.addToConfig("  ntwx=%d," % restart_interval)   # Trajectory sampling frequency.
            self.addToConfig("  irest=0,")                      # Don't restart.
            self.addToConfig("  dt=%.3f," % timestep)           # Time step.
            self.addToConfig("  nstlim=%d," % steps)            # Number of integration steps.
            self.addToConfig("  ntc=2,")                        # Enable SHAKE.
            self.addToConfig("  ntf=2,")                        # Don't calculate forces for constrained bonds.
            self.addToConfig("  ntt=3,")                        # Langevin dynamics.
            self.addToConfig("  gamma_ln=2,")                   # Collision frequency (ps).
            if not has_box or not self._has_water:
                self.addToConfig("  ntb=0,")                    # No periodic box.
                self.addToConfig("  cut=999.,")                 # Non-bonded cut-off.
                if is_pmemd:
                    self.addToConfig("  igb=6,")                # Use vacuum generalised Born model.
            else:
                self.addToConfig("  cut=8.0,")                  # Non-bonded cut-off.
            self.addToConfig("  tempi=%.2f,"                    # Initial temperature.
                % self._protocol.getTemperature().kelvin().value())
            self.addToConfig("  temp0=%.2f,"                    # Target temperature.
                % self._protocol.getTemperature().kelvin().value())

            # Constant pressure control.
            if self._protocol.getPressure() is not None:
                # Don't use barostat for vacuum simulations.
                if has_box and self._has_water:
                    self.addToConfig("  ntp=1,")                # Isotropic pressure scaling.
                    self.addToConfig("  pres0=%.5f,"            # Pressure in bar.
                        % self._protocol.getPressure().bar().value())
                else:
                    _warnings.warn("Cannot use a barostat for a vacuum or non-periodic simulation")

            # Activate PLUMED and locate the plumed.dat file.
            self.addToConfig("  plumed=1,")
            self.addToConfig("  plumedfile='plumed.dat',")

            self.addToConfig(" /")

            # Create the PLUMED input file and copy auxiliary files to the working directory.
            self._plumed = _Plumed(self._work_dir)
            plumed_config, auxiliary_files = self._plumed.createConfig(self._system,
                                                                       self._protocol,
                                                                       self._property_map)
            self._setPlumedConfig(plumed_config)
            if auxiliary_files is not None:
                for file in auxiliary_files:
                    file_name = _os.path.basename(file)
                    _shutil.copyfile(file, self._work_dir + f"/{file_name}")
            self._input_files.append(self._plumed_config_file)

            # Expose the PLUMED specific member functions.
            setattr(self, "getPlumedConfig", self._getPlumedConfig)
            setattr(self, "getPlumedConfigFile", self._getPlumedConfigFile)
            setattr(self, "setPlumedConfig", self._setPlumedConfig)
            setattr(self, "getFreeEnergy", self._getFreeEnergy)
            setattr(self, "getCollectiveVariable", self._getCollectiveVariable)
            setattr(self, "sampleConfigurations", self._sampleConfigurations)
            setattr(self, "getTime", self._getTime)

        # Add configuration variables for a steered molecular dynamics simulation.
        elif isinstance(self._protocol, _Protocol.Steering):

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

            # Set the random number seed.
            if self._seed is None:
                seed = -1
            else:
                seed = self._seed

            # Convert the timestep to picoseconds.
            timestep = self._protocol.getTimeStep().picoseconds().value()

            self.addToConfig("Production.")
            self.addToConfig(" &cntrl")
            self.addToConfig("  ig=%d," % seed)                 # Random number seed.
            self.addToConfig("  ntx=1,")                        # Only read coordinates.
            self.addToConfig("  ntxo=1,")                       # Output coordinates in ASCII.
            self.addToConfig("  ntpr=%d," % report_interval)    # Interval between reporting energies.
            self.addToConfig("  ntwr=%d," % restart_interval)   # Interval between saving restart files.
            self.addToConfig("  ntwx=%d," % restart_interval)   # Trajectory sampling frequency.
            self.addToConfig("  irest=0,")                      # Don't restart.
            self.addToConfig("  dt=%.3f," % timestep)           # Time step.
            self.addToConfig("  nstlim=%d," % steps)            # Number of integration steps.
            self.addToConfig("  ntc=2,")                        # Enable SHAKE.
            self.addToConfig("  ntf=2,")                        # Don't calculate forces for constrained bonds.
            self.addToConfig("  ntt=3,")                        # Langevin dynamics.
            self.addToConfig("  gamma_ln=2,")                   # Collision frequency (ps).
            if not has_box or not self._has_water:
                self.addToConfig("  ntb=0,")                    # No periodic box.
                self.addToConfig("  cut=999.,")                 # Non-bonded cut-off.
                if is_pmemd:
                    self.addToConfig("  igb=6,")                # Use vacuum generalised Born model.
            else:
                self.addToConfig("  cut=8.0,")                  # Non-bonded cut-off.
            self.addToConfig("  tempi=%.2f,"                    # Initial temperature.
                % self._protocol.getTemperature().kelvin().value())
            self.addToConfig("  temp0=%.2f,"                    # Target temperature.
                % self._protocol.getTemperature().kelvin().value())

            # Constant pressure control.
            if self._protocol.getPressure() is not None:
                # Don't use barostat for vacuum simulations.
                if has_box and self._has_water:
                    self.addToConfig("  ntp=1,")                # Isotropic pressure scaling.
                    self.addToConfig("  pres0=%.5f,"            # Pressure in bar.
                        % self._protocol.getPressure().bar().value())
                else:
                    _warnings.warn("Cannot use a barostat for a vacuum or non-periodic simulation")

            # Activate PLUMED and locate the plumed.dat file.
            self.addToConfig("  plumed=1,")
            self.addToConfig("  plumedfile='plumed.dat',")

            self.addToConfig(" /")

            # Create the PLUMED input file and copy auxiliary files to the working directory.
            self._plumed = _Plumed(self._work_dir)
            plumed_config, auxiliary_files = self._plumed.createConfig(self._system,
                                                                       self._protocol,
                                                                       self._property_map)
            self._setPlumedConfig(plumed_config)
            if auxiliary_files is not None:
                for file in auxiliary_files:
                    file_name = _os.path.basename(file)
                    _shutil.copyfile(file, self._work_dir + f"/{file_name}")
            self._input_files.append(self._plumed_config_file)

            # Expose the PLUMED specific member functions.
            setattr(self, "getPlumedConfig", self._getPlumedConfig)
            setattr(self, "getPlumedConfigFile", self._getPlumedConfigFile)
            setattr(self, "setPlumedConfig", self._setPlumedConfig)
            setattr(self, "getCollectiveVariable", self._getCollectiveVariable)
            setattr(self, "getTime", self._getTime)

        else:
            raise _IncompatibleError("Unsupported protocol: '%s'" % self._protocol.__class__.__name__)

        # Flag that this isn't a custom protocol.
        self._protocol._setCustomised(False)

    def _generate_args(self):
        """Generate the dictionary of command-line arguments."""

        # Clear the existing arguments.
        self.clearArgs()

        # Add the default arguments.
        self.setArg("-O", True)                             # Overwrite.
        self.setArg("-i", "%s.cfg"   % self._name)          # Input file.
        self.setArg("-p", "%s.prm7"  % self._name)          # Topology file.
        self.setArg("-c", "%s.rst7"  % self._name)          # Coordinate file.
        self.setArg("-o", "%s.out"   % self._name)          # Redirect stdout to file.
        self.setArg("-r", "%s.crd"   % self._name)          # Restart file.
        self.setArg("-inf", "%s.nrg" % self._name)          # Energy info file.

        # Skip if the user has passed a custom protocol.
        if not isinstance(self._protocol, _Protocol.Custom):

            # Append a reference file if this a restrained equilibration.
            if isinstance(self._protocol, _Protocol.Equilibration):
                if self._protocol.getRestraint() is not None:
                    self.setArg("-ref", "%s.rst7" % self._name)

            # Append a trajectory file if this anything other than a minimisation.
            if not isinstance(self._protocol, _Protocol.Minimisation):
                self.setArg("-x", "%s.nc" % self._name)

    def start(self):
        """Start the AMBER process.

           Returns
           -------

           process : :class:`Process.Amber <BioSimSpace.Process.Amber>`
               The process object.
        """

        # The process is currently queued.
        if self.isQueued():
            return

        # Process is already running.
        if self._process is not None:
            if self._process.isRunning():
                return

        # Reset the watcher.
        self._is_watching = False

        # Run the process in the working directory.
        with _Utils.cd(self._work_dir):

            # Create the arguments string list.
            args = self.getArgStringList()

            # Write the command-line process to a README.txt file.
            with open("README.txt", "w") as file:

                # Set the command-line string.
                self._command = "%s " % self._exe + self.getArgString()

                # Write the command to file.
                file.write("# AMBER was run with the following command:\n")
                file.write("%s\n" % self._command)

            # Start the timer.
            self._timer = _timeit.default_timer()

            # Start the simulation. Pass a null string for the stdout file
            # since we've explicitly redirected AMBER output to file since
            # pmemd doesn't write to standard output.
            self._process = _SireBase.Process.run(self._exe, args,
                "", "%s.err"  % self._name)

        # Watch the energy info file for changes.
        self._watcher = _Watcher(self)
        self._watcher.start()

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

        # Create the name of the restart CRD file.
        restart = "%s/%s.crd" % (self._work_dir, self._name)

        # Check that the file exists.
        if _os.path.isfile(restart):
            # Do we need to get coordinates for the lambda=1 state.
            if "is_lambda1" in self._property_map:
                is_lambda1 = True
            else:
                is_lambda1 = False

            # Create a new molecular system from the restart file.
            new_system = _System(_SireIO.MoleculeParser.read([restart, self._top_file], self._property_map))

            # Create a copy of the existing system object.
            old_system = self._system.copy()

            # Update the coordinates and velocities and return a mapping between
            # the molecule indices in the two systems.
            sire_system, mapping = _SireIO.updateCoordinatesAndVelocities(
                    old_system._sire_object,
                    new_system._sire_object,
                    self._mapping,
                    is_lambda1,
                    self._property_map,
                    self._property_map)

            # Update the underlying Sire object.
            old_system._sire_object = sire_system

            # Store the mapping between the MolIdx in both systems so we don't
            # need to recompute it next time.
            self._mapping = mapping

            # Update the box information in the original system.
            if "space" in new_system._sire_object.propertyKeys():
                box = new_system._sire_object.property("space")
                old_system._sire_object.setProperty(self._property_map.get("space", "space"), box)

            return old_system

        else:
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

           trajectory : :class:`Trajectory <BioSimSpace.Trajectory.Trajectory>`
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

        try:
            return _Trajectory.Trajectory(process=self)

        except:
            return None

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

        if not type(index) is int:
            raise TypeError("'index' must be of type 'int'")

        max_index = int((self._protocol.getRunTime() / self._protocol.getTimeStep())
                  / self._protocol.getRestartInterval())

        if index < 0 or index > max_index:
            raise ValueError(f"'index' must be in range [0, {max_index}].")

        try:
            # Do we need to get coordinates for the lambda=1 state.
            if "is_lambda1" in self._property_map:
                is_lambda1 = True
            else:
                is_lambda1 = False

            # Get the latest trajectory frame.
            new_system =  _Trajectory.getFrame(self._traj_file,
                                               self._top_file,
                                               index)

            # Create a copy of the existing system object.
            old_system = self._system.copy()

            # Update the coordinates and velocities and return a mapping between
            # the molecule indices in the two systems.
            sire_system, mapping = _SireIO.updateCoordinatesAndVelocities(
                    old_system._sire_object,
                    new_system._sire_object,
                    self._mapping,
                    is_lambda1,
                    self._property_map,
                    self._property_map)

            # Update the underlying Sire object.
            old_system._sire_object = sire_system

            # Store the mapping between the MolIdx in both systems so we don't
            # need to recompute it next time.
            self._mapping = mapping

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

        return self._get_stdout_record(record.strip().upper(), time_series, unit)

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

        return self._get_stdout_record(record.strip().upper(), time_series, unit)

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

           Returns
           -------

           records : :class:`MultiDict <BioSimSpace.Process._process._MultiDict>`
              The dictionary of time-series records.
        """
        return self.getRecords(block=False)

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

        # No time records for minimisation protocols.
        if isinstance(self._protocol, _Protocol.Minimisation):
            return None

        # Get the list of time steps.
        time_steps = self.getRecord("TIME(PS)", time_series, None, block)

        # Convert from picoseconds to nanoseconds.
        if time_steps is not None:
            if time_series:
                return [(x * _Units.Time.picosecond)._to_default_unit() for x in time_steps]
            else:
                return (time_steps * _Units.Time.picosecond)._to_default_unit()

    def getCurrentTime(self, time_series=False):
        """Get the current simulation time.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

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
        return self.getRecord("NSTEP", time_series, None, block)

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

    def getBondEnergy(self, time_series=False, block="AUTO"):
        """Get the bond energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
              The bond energy.
        """
        return self.getRecord("BOND", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentBondEnergy(self, time_series=False):
        """Get the current bond energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
              The bond energy.
        """
        return self.getBondEnergy(time_series, block=False)

    def getAngleEnergy(self, time_series=False, block="AUTO"):
        """Get the angle energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
              The angle energy.
        """
        return self.getRecord("ANGLE", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentAngleEnergy(self, time_series=False):
        """Get the current angle energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
              The angle energy.
        """
        return self.getAngleEnergy(time_series, block=False)

    def getDihedralEnergy(self, time_series=False, block="AUTO"):
        """Get the total dihedral energy (proper + improper).

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
              The total dihedral energy.
        """
        return self.getRecord("DIHED", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentDihedralEnergy(self, time_series=False):
        """Get the current total dihedral energy (proper + improper).

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
              The total dihedral energy.
        """
        return self.getDihedralEnergy(time_series, block=False)

    def getElectrostaticEnergy(self, time_series=False, block="AUTO"):
        """Get the electrostatic energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
              The electrostatic energy.
        """
        return self.getRecord("EELECT", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentElectrostaticEnergy(self, time_series=False):
        """Get the current dihedral energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
              The electrostatic energy.
        """
        return self.getElectrostaticEnergy(time_series, block=False)

    def getElectrostaticEnergy14(self, time_series=False, block="AUTO"):
        """Get the electrostatic energy between atoms 1 and 4.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
              The electrostatic energy between atoms 1 and 4.
        """
        return self.getRecord("1-4 EEL", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentElectrostaticEnergy14(self, time_series=False):
        """Get the current electrostatic energy between atoms 1 and 4.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
              The electrostatic energy between atoms 1 and 4.
        """
        return self.getElectrostaticEnergy14(time_series, block=False)

    def getVanDerWaalsEnergy(self, time_series=False, block="AUTO"):
        """Get the Van der Vaals energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
              The Van der Vaals energy.
        """
        return self.getRecord("VDWAALS", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentVanDerWaalsEnergy(self, time_series=False):
        """Get the current Van der Vaals energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
              The Van der Vaals energy.
        """
        return self.getVanDerWaalsEnergy(time_series, block=False)

    def getHydrogenBondEnergy(self, time_series=False, block="AUTO"):
        """Get the hydrogen bond energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
              The hydrogen bond energy.
        """
        return self.getRecord("EHBOND", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentHydrogenBondEnergy(self, time_series=False):
        """Get the current hydrogen bond energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
              The hydrogen bond energy.
        """
        return self.getHydrogenBondEnergy(time_series, block=False)

    def getRestraintEnergy(self, time_series=False, block="AUTO"):
        """Get the restraint energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
              The restraint energy.
        """
        return self.getRecord("RESTRAINT", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentRestraintEnergy(self, time_series=False):
        """Get the current restraint energy.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
              The restraint energy.
        """
        return self.getRestraintEnergy(time_series, block=False)

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
        return self.getRecord("EPTOT", time_series, _Units.Energy.kcal_per_mol, block)

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
        return self.getRecord("EKTOT", time_series, _Units.Energy.kcal_per_mol, block)

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

    def getNonBondedEnergy14(self, time_series=False, block="AUTO"):
        """Get the non-bonded energy between atoms 1 and 4.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
              The non-bonded energy between atoms 1 and 4.
        """
        return self.getRecord("1-4 NB", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentNonBondedEnergy14(self, time_series=False):
        """Get the current non-bonded energy between atoms 1 and 4.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
              The non-bonded energy between atoms 1 and 4.
        """
        return self.getNonBondedEnergy14(time_series, block=False)

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
        if isinstance(self._protocol, _Protocol.Minimisation):
            return self.getRecord("ENERGY", time_series, _Units.Energy.kcal_per_mol, block)
        else:
            return self.getRecord("ETOT", time_series, _Units.Energy.kcal_per_mol, block)

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

    def getCentreOfMassKineticEnergy(self, time_series=False, block="AUTO"):
        """Get the kinetic energy of the centre of mass in translation.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
              The centre of mass kinetic energy.
        """
        return self.getRecord("EKCMT", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentCentreOfMassKineticEnergy(self, time_series=False):
        """Get the current kinetic energy of the centre of mass in translation.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           energy : :class:`Energy <BioSimSpace.Types.Energy>`
              The centre of mass kinetic energy.
        """
        return self.getCentreOfMassKineticEnergy(time_series, block=False)

    def getVirial(self, time_series=False, block="AUTO"):
        """Get the virial.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           virial : float
              The virial.
        """
        return self.getRecord("VIRIAL", time_series, block)

    def getCurrentVirial(self, time_series=False):
        """Get the current virial.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           virial : float
              The virial.
        """
        return self.getVirial(time_series, block=False)

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
        return self.getRecord("TEMP(K)", time_series, _Units.Temperature.kelvin, block)

    def getCurrentTemperature(self, time_series=False):
        """Get the current temperature.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
              The temperature.
        """
        return self.getTemperature(time_series, block=False)

    def getPressure(self, time_series=False, block="AUTO"):
        """Get the pressure.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
              The pressure.
        """
        return self.getRecord("PRESS", time_series, _Units.Pressure.bar, block)

    def getCurrentPressure(self, time_series=False):
        """Get the current pressure.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
              The pressure.
        """
        return self.getPressure(time_series, block=False)

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
        return self.getRecord("VOLUME", time_series, _Units.Volume.angstrom3, block)

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

    def getDensity(self, time_series=False, block="AUTO"):
        """Get the density.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.

           Returns
           -------

           density : float
              The density.
        """
        return self.getRecord("DENSITY", time_series, block)

    def getCurrentDensity(self, time_series=False):
        """Get the current density.

           Parameters
           ----------

           time_series : bool
               Whether to return a list of time series records.

           Returns
           -------

           density : float
              The density.
        """
        return self.getDensity(time_series, block=False)

    def _update_energy_dict(self):
        """Read the energy info file and update the dictionary."""

        # Return if the file doesn't exist, e.g. it could have been deleted.
        if not _os.path.isfile(self._nrg_file):
            return

        # Flag that this isn't a header line.
        is_header = False

        # Open the file for reading.
        with open(self._nrg_file, "r") as file:

            # Loop over all of the lines.
            for line in file:

                # Skip empty lines and summary reports.
                if len(line) > 0 and line[0] != "|":

                    # The output format is different for minimisation protocols.
                    if isinstance(self._protocol, _Protocol.Minimisation):

                        # No equals sign in the line.
                        if "=" not in line:

                            # Split the line using whitespace.
                            data = line.upper().split()

                            # If we find a header, jump to the top of the loop.
                            if len(data) > 0:
                                if data[0] == "NSTEP":
                                    is_header = True
                                    continue

                        # Process the header record.
                        if is_header:

                            # Split the line using whitespace.
                            data = line.upper().split()

                            # The file hasn't been updated.
                            if "NSTEP" in self._stdout_dict and data[0] == self._stdout_dict["NSTEP"][-1]:
                                return

                            else:
                                # Add the timestep and energy records to the dictionary.
                                self._stdout_dict["NSTEP"] = data[0]
                                self._stdout_dict["ENERGY"] = data[1]

                                # Turn off the header flag now that the data has been recorded.
                                is_header = False

                    # All other protocols have output that is formatted as RECORD = VALUE.

                    # Use a regex search to split the line into record names and values.
                    records = _re.findall(r"(\d*\-*\d*\s*[A-Z]+\(*[A-Z]*\)*)\s*=\s*(\-*\d+\.?\d*)", line.upper())

                    # Append each record to the dictionary.
                    for key, value in records:

                        # Strip whitespace from the record key.
                        key = key.strip()

                        # The file hasn't been updated.
                        if key == "NSTEP":
                            if "NSTEP" in self._stdout_dict and value == self._stdout_dict["NSTEP"][-1]:
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
        if not self._process is None and self._process.isRunning():
            self._process.kill()

    def wait(self, max_time=None):
        """Wait for the process to finish.

           Parameters
           ----------

           max_time: :class:`Time <BioSimSpace.Types.Time>`, int, float
               The maximum time to wait (in minutes).
        """

        # The process isn't running.
        if not self.isRunning():
            return

        if max_time is not None:
            # Convert int to float.
            if type(max_time) is int:
                max_time = float(max_time)

            # BioSimSpace.Types.Time
            if isinstance(max_time, _Type):
                max_time = max_time.minutes().value()

            # Float.
            elif isinstance(max_time, float):
                if max_time <= 0:
                    raise ValueError("'max_time' cannot be negative!")

            else:
                raise TypeError("'max_time' must be of type 'BioSimSpace.Types.Time' or 'float'.")

        # Loop until the process is finished.
        # For some reason we can't use Sire.Base.Process.wait() since it
        # doesn't work properly with the background threads used for the
        # watchdog observer.
        while self._process.isRunning():
            _time.sleep(1)

            # The maximum run time has been exceeded, kill the job.
            if max_time is not None:
                if self.runTime().value() > max_time:
                    self.kill()
                    return

        # Stop and join the watchdog observer.
        self._watcher._observer.stop()
        self._watcher._observer.join()

    def _get_stdout_record(self, key, time_series=False, unit=None):
        """Helper function to get a stdout record from the dictionary.

           Parameters
           ----------

           key : str
               The record key.

           time_series : bool
               Whether to return a time series of records.

           unit : :class:`Type <BioSimSpace.Types._type.Type>`
               The unit to convert the record to.

           Returns
           -------

           record :
               The matching stdout record.
        """

        # No data!
        if len(self._stdout_dict) == 0:
            return None

        if not isinstance(time_series, bool):
            _warnings.warn("Non-boolean time-series flag. Defaulting to False!")
            time_series = False

        # Validate the unit.
        if unit is not None:
            if not isinstance(unit, _Type):
                raise TypeError("'unit' must be of type 'BioSimSpace.Types'")

        # Return the list of dictionary values.
        if time_series:
            try:
                if key == "NSTEP":
                    return [int(x) for x in self._stdout_dict[key]]
                else:
                    if unit is None:
                        return [float(x) for x in self._stdout_dict[key]]
                    else:
                        return [(float(x) * unit)._to_default_unit() for x in self._stdout_dict[key]]

            except KeyError:
                return None

        # Return the most recent dictionary value.
        else:
            try:
                if key == "NSTEP":
                    return int(self._stdout_dict[key][-1])
                else:
                    if unit is None:
                        return float(self._stdout_dict[key][-1])
                    else:
                        return (float(self._stdout_dict[key][-1]) * unit)._to_default_unit()

            except KeyError:
                return None

    def _create_restraint_mask(self, atom_idxs):
        """Internal helper function to create an AMBER restraint mask from a
           list of atom indices.

           Parameters
           ----------

           atom_idxs : [int]
               A list of atom indices.

           Returns
           -------

           restraint_mask : str
               The AMBER restraint mask.
        """

        if not isinstance(atom_idxs, (list, tuple)):
            raise TypeError("'atom_idxs' must be a list of 'int' types.")

        if not all(type(x) is int for x in atom_idxs):
            raise TypeError("'atom_idxs' must be a list of 'int' types.")

        # AMBER has a restriction on the number of characters in the restraint
        # mask (not documented) so we can't just use comma-separated atom
        # indices. Instead we loop through the indices and use hyphens to
        # separate contiguous blocks of indices, e.g. 1-23,34-47,...

        # Create a set to sort and ensure no duplicates, then convert back to a list.
        # This should already by done, but do so again in case the user is accessing
        # the method directly.
        atom_idxs = list(set(atom_idxs))
        atom_idxs.sort()

        # Handle single atom restraints differently.
        if len(atom_idxs) == 1:
            restraint_mask =  f"@{atom_idxs[0]+1}"

        else:
            # Start the mask with the first atom index. (AMBER is 1 indexed.)
            restraint_mask = f"@{atom_idxs[0]+1}"

            # Store the current index.
            prev_idx = atom_idxs[0]

            # Store the lead index for this block.
            lead_idx = prev_idx

            # Loop over all other indices.
            for idx in atom_idxs[1:]:
                # There is a gap in the indices.
                if idx - prev_idx > 1:
                    if prev_idx != lead_idx:
                        restraint_mask += f"{prev_idx+1},{idx+1}"
                    else:
                        restraint_mask += f",{idx+1}"
                    lead_idx = idx
                else:
                    # This is the first index beyond the lead.
                    if idx - lead_idx == 1:
                        restraint_mask += "-"
                # Set the value of the previous index.
                prev_idx = idx

            # Add the final atom to the mask.
            if idx - atom_idxs[-2] == 1:
                restraint_mask += f"{idx+1}"
            else:
                if idx != lead_idx:
                    restraint_mask += f",{idx+1}"

        return restraint_mask
