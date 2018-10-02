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
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BioSimSpace. If not, see <http://www.gnu.org/licenses/>.
#####################################################################

"""
Functionality for running simulations with SOMD.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

import Sire as _Sire

from . import _process
from .._Exceptions import IncompatibleError as _IncompatibleError
from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from .._SireWrappers import System as _System
from ..Trajectory import Trajectory as _Trajectory

import BioSimSpace.Protocol as _Protocol
import BioSimSpace.Types._type as _Type
import BioSimSpace.Units as _Units

import math as _math
import os as _os
import timeit as _timeit
import warnings as _warnings

__all__ = ["Somd"]

class Somd(_process.Process):
    """A class for running simulations using SOMD."""

    def __init__(self, system, protocol, exe=None, name="somd",
            platform="GPU", work_dir=None, seed=None, map={}):
        """Constructor.

           Positional arguments
           --------------------

           system : BioSimSpace._SireWrappers.System
               The molecular system.

           protocol : BioSimSpace.Protocol
               The protocol for the SOMD process.


           Keyword arguments
           -----------------

           exe : str
               The full path to the SOMD executable.

           name : str
               The name of the process.

           platform : str
               The platform for the simulation: "GPU" or "CPU".

           work_dir :
               The working directory for the process.

           seed : int
               A random number seed.

           map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(system, protocol, name, work_dir, seed, map)

        # Set the package name.
        self._package_name = "SOMD"

        # This process can generate trajectory data.
        self._has_trajectory = True

        if type(platform) is not str:
            raise TypeError("'platform' must be of type 'str'.")
        else:
            # Strip all whitespace and convert to upper case.
            platform = platform.replace(" ", "").upper()

            # Check for platform support.
            if platform not in ["GPU", "CPU"]:
                raise ValueError("Supported platforms are: 'GPU' and 'CPU'")
            else:
                self._platform = platform

        # If the path to the executable wasn't specified, then use the bundled SOMD
        # executable.
        if exe is None:
            # Generate the name of the SOMD exe.
            somd_exe = _Sire.Base.getBinDir() + "/somd"
            if not _os.path.isfile(somd_exe):
                raise _MissingSoftwareError("'Cannot find SOMD executable in expected location: '%s'" % somd_exe)
            else:
                self._exe = somd_exe
        else:
            # Make sure executable exists.
            if _os.path.isfile(exe):
                self._exe = exe
            else:
                raise IOError("SOMD executable doesn't exist: '%s'" % exe)

        # The names of the input files.
        self._rst_file = "%s/%s.rst7" % (self._work_dir, name)
        self._top_file = "%s/%s.prm7" % (self._work_dir, name)

        # Set the path for the SOMD configuration file.
        self._config_file = "%s/%s.txt" % (self._work_dir, name)

        # Create the list of input files.
        self._input_files = [self._config_file, self._rst_file, self._top_file]

        # Now set up the working directory for the process.
        self._setup()

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.Process.%s: system=%s, protocol=%s, exe='%s', name='%s', platform='%s', work_dir='%s' seed=%s>" \
            % (self.__class__.__name__, str(_System(self._system)), self._protocol.__repr__(),
               self._exe, self._name, self._platform, self._work_dir, self._seed)

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return "BioSimSpace.Process.%s(%s, %s, exe='%s', name='%s', platform='%s', work_dir='%s', seed=%s)" \
            % (self.__class__.__name__, str(_System(self._system)), self._protocol.__repr__(),
               self._exe, self._name, self._platform, self._work_dir, self._seed)

    def _setup(self):
        """Setup the input files and working directory ready for simulation."""

        # Create the input files...

        # RST file (coordinates).
        try:
            rst = _Sire.IO.AmberRst7(self._system)
            rst.writeToFile(self._rst_file)
        except:
            raise IOError("Failed to write system to 'RST7' format.") from None

        # PRM file (topology).
        try:
            prm = _Sire.IO.AmberPrm(self._system)
            prm.writeToFile(self._top_file)
        except:
            raise IOError("Failed to write system to 'PRM7' format.") from None

        # Generate the SOMD configuration file.
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
        """Generate SOMD configuration file strings."""

        # Clear the existing configuration list.
        self._config = []

        # Check whether the system contains periodic box information.
        # For now, well not attempt to generate a box if the system property
        # is missing. If no box is present, we'll assume a non-periodic simulation.
        if "space" in self._system.propertyKeys():
            has_box = True
        else:
            _warnings.warn("No simulation box found. Assuming gas phase simulation.")
            has_box = False

        # Add configuration variables for a minimisation simulation.
        if type(self._protocol) is _Protocol.Minimisation:
            self.addToConfig("minimise = True")                         # Minimisation simulation.
            self.addToConfig("minimise maximum iterations = %d"         # Maximum number of steps.
                % self._protocol.getSteps())
            self.addToConfig("minimise tolerance = 1")                  # Convergence tolerance.
            self.addToConfig("ncycles = 1")                             # Perform a single SOMD cycle.
            self.addToConfig("nmoves = 1")                              # Perform a single MD move.
            self.addToConfig("save coordinates = True")                 # Save molecular coordinates.
            if not has_box:
                self.addToConfig("cutoff type = cutofftypenonperiodic") # No periodic box.
                self.addToConfig("cutoff distance = 1000 angstrom")     # Non-bonded cut-off.

        # Add configuration variables for an equilibration simulation.
        elif type(self._protocol) is _Protocol.Equilibration:
            # Only constant temperature equilibration simulations are supported.
            if not self._protocol.isConstantTemp():
                raise _IncompatibleError("SOMD only supports constant temperature equilibration.")

            # Work out the number of cycles. We save coordinates every cycle,
            # which is 100 MD steps (moves) in length (this is for consistency
            # with other MD drivers). Note that SOMD only save coordinates to
            # a DCD trajectory file, so it's impossible to decouple the
            # frequency of recording configurations and trajectory frames,
            # i.e. the number of trajectory frames specified in the protocol
            # is disregarded.
            ncycles = _math.ceil((self._protocol.getRunTime() / self._protocol.getTimeStep()) / 100)

            # Convert the timestep to picoseconds.
            timestep = self._protocol.getTimeStep().picoseconds().magnitude()

            self.addToConfig("ncycles = %d" % ncycles)                  # The number of SOMD cycles.
            self.addToConfig("nmoves = 100")                            # Perform 100 MD moves per cycle.
            self.addToConfig("save coordinates = True")                 # Save molecular coordinates.
            self.addToConfig("buffered coordinates frequency = 100")    # Save coordinates every 100 steps.
                                                                        # System temperature.
            self.addToConfig("temperature = %.2f kelvin" \
                % self._protocol.getStartTemperature().kelvin().magnitude())
            if not has_box:
                self.addToConfig("cutoff type = cutofftypenonperiodic") # No periodic box.
                self.addToConfig("cutoff distance = 1000 angstrom")     # Non-bonded cut-off.
            if self._is_seeded:
                self.addToConfig("random seed = %d" % self._seed)       # Random number seed.

        else:
            raise _IncompatibleError("Unsupported protocol: '%s'" % self._protocol.__class__.__name__)

    def _generate_args(self):
        """Generate the dictionary of command-line arguments."""

        # Clear the existing arguments.
        self.clearArgs()

        # Add the default arguments.
        self.setArg("-c", self._rst_file)       # Coordinate restart file.
        self.setArg("-t", self._top_file)       # Topology file.
        self.setArg("-C", self._config_file)    # Config file.
        if self._platform == "GPU":             # Platform.
            self.setArg("-p", "CUDA")
        else:
            self.setArg("-p", "CPU")

    def start(self):
        """Start the SOMD process."""

        # Process is already running.
        if self._process is not None:
            if self._process.isRunning():
                return

        # Clear any existing output.
        self._clear_output()

        # Store the current working directory.
        dir = _os.getcwd()

        # Change to the working directory for the process.
        # This avoids problems with relative paths.
        _os.chdir(self._work_dir)

        # Create the arguments string list.
        args = self.getArgStringList()

        # Write the command-line process to a README.txt file.
        with open("README.txt", "w") as f:

            # Set the command-line string.
            self._command = "%s " % self._exe + self.getArgString()

            # Write the command to file.
            f.write("# SOMD was run with the following command:\n")
            f.write("%s\n" % self._command)

        # Start the timer.
        self._timer = _timeit.default_timer()

        # Start the simulation.
        self._process = _Sire.Base.Process.run(self._exe, args,
            "%s.out"  % self._name, "%s.err"  % self._name)

        # Change back to the original working directory.
        _os.chdir(dir)

        return self

    def getSystem(self, block="AUTO"):
        """Get the latest molecular system.

           Keyword arguments
           -----------------

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           system : BioSimSpace._SireWrappers.System
               The latest molecular system.
        """

        # Get the trajectory object.
        traj = self.getTrajectory(block=block)

        # Try to get the latest frame from the trajectory.
        try:
            return traj.getFrames()[-1]

        except:
            return None

    def getCurrentSystem(self):
        """Get the latest molecular system.

           Returns
           -------

           system : BioSimSpace._SireWrappers.System
               The latest molecular system.
        """
        return self.getSystem(block=False)

    def getTrajectory(self, block="AUTO"):
        """Return a trajectory object.

           Keyword arguments
           -----------------

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           trajectory : BioSimSpace.Trajectory
        """

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        return _Trajectory(process=self)

    def _clear_output(self):
        """Reset stdout and stderr."""

        # Call the base class method.
        super()._clear_output()

        # Delete any simulation restart files in the working
        # directory.
        restart_file = "%s/sim_restart.s3"
        if _os.path.isfile(restart_file):
            _os.remove(restart_file)
