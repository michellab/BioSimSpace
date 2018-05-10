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
Functionality for running simulations with GROMACS.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

import Sire as _Sire

from . import _process

from .._System import System as _System
from ..Trajectory import Trajectory as _Trajectory

import BioSimSpace.Protocol as _Protocol

import math as _math
import os as _os
import subprocess as _subprocess
import timeit as _timeit
import warnings as _warnings

try:
    _pygtail = _Sire.try_import("pygtail")
except ImportError:
    raise ImportError("Pygtail is not installed. Please install pygtail in order to use BioSimSpace.")

__all__ = ["Gromacs"]

class Gromacs(_process.Process):
    """A class for running simulations using GROMACS."""

    def __init__(self, system, protocol, exe=None, name="gromacs",
            work_dir=None, seed=None):
        """Constructor.

           Positional arguments:

           system        -- The molecular system.
           protocol      -- The protocol for the GROMACS process.

           Keyword arguments:

           exe           -- The full path to the GROMACS executable.
           name          -- The name of the process.
           work_dir      -- The working directory for the process.
           seed          -- A random number seed.
        """

        # Call the base class constructor.
        super().__init__(system, protocol, name, work_dir, seed)

        # Set the package name.
        self._package_name = "GROMACS"

        # This process can generate trajectory data.
        self._has_trajectory = True

        if exe is None:
            # Search Sire bin directory.
            bin_dir = _Sire.Base.getBinDir()
            exe = "%s/gmx" % bin_dir

            if _os.path.isfile(exe):
                # Although the executable exists, it may not work because it was
                # precompiled on a system with different hardware instructions.
                # We test this by running the executable and checking the error
                # code.
                proc = _subprocess.run(exe, shell=True,
                    stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)

                # The executable runs.
                if proc.returncode == 0:
                    self._exe = exe
                # Search the system PATH.
                else:
                    self._exe = _Sire.Base.findExe("gmx").absoluteFilePath()

            # Search the system PATH.
            else:
                self._exe = _Sire.Base.findExe("gmx").absoluteFilePath()

        else:
            # Make sure executable exists.
            if _os.path.isfile(exe):
                self._exe = exe
            else:
                raise IOError("GROMACS executable doesn't exist: '%s'" % exe)

        # Now use the GROMACS exe to get the location of the data directory.

        # Generate the shell command.
        command = "%s -h 2>&1 | grep 'Data prefix' | awk -F ':' '{print $2}'" % self._exe

        # Run the command.
        proc = _subprocess.run(command, shell=True, stdout=_subprocess.PIPE)

        # Get the data prefix.
        if proc.returncode != 0:
            raise RuntimeError("Unable to determine GROMACS data prefix!")
        else:
            self._data_prefix = proc.stdout.decode("ascii").strip()

        # Initialise the stdout dictionary and title header.
        self._stdout_dict = _process._MultiDict()
        self._stdout_title = None

        # The names of the input files.
        self._gro_file = "%s/%s.gro" % (self._work_dir, name)
        self._top_file = "%s/%s.top" % (self._work_dir, name)

        # Set the path for the GROMACS configuration file.
        self._config_file = "%s/%s.mdp" % (self._work_dir, name)

        # Create the list of input files.
        self._run_files = [self._config_file, self._gro_file, self._top_file]

        # Now set up the working directory for the process.
        self._setup()

    def _setup(self):
        """Setup the input files and working directory ready for simulation."""

        # Create the input files...

        # GRO87 file.
        gro = _Sire.IO.Gro87(self._system)
        gro.writeToFile(self._gro_file)

        # TOP file.
        top = _Sire.IO.GroTop(self._system)
        top.writeToFile(self._top_file)

        # Create the binary input file name.
        self._tpr_file = "%s/%s.tpr" % (self._work_dir, self._name)
        self._run_files.append(self._tpr_file)

        # Generate the GROMACS configuration file.
        # Skip if the user has passed a custom config.
        if type(self._protocol) is _Protocol.Custom:
            self.setConfig(self._protocol.getConfig())
        else:
            self._generate_config()
        self.writeConfig(self._config_file)

        # Generate the dictionary of command-line arguments.
        self._generate_args()

        # Return the list of input files.
        return self._run_files

    def _generate_config(self):
        """Generate GROMACS configuration file strings."""

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
            self.addToConfig("integrator = cg")             # Use conjugate gradient.
            self.addToConfig("nsteps = %d"
                % self._protocol.getSteps())                # Set the number of steps.
            self.addToConfig("nstxout = 10")                # Write coordinates every 10 steps.
            self.addToConfig("nstlog = 10")                 # Write energies every 10 steps.
            self.addToConfig("cutoff-scheme = Verlet")      # Use Verlet pair lists.
            self.addToConfig("ns-type = grid")              # Use a grid to search for neighbours.
            if has_box:
                self.addToConfig("pbc = xyz")               # Simulate a fully periodic box.
            self.addToConfig("coulombtype = PME")           # Fast smooth Particle-Mesh Ewald.
            self.addToConfig("DispCorr = EnerPres")         # Dispersion corrections for energy and pressure.

    def start(self):
        """Start the GROMACS process."""

        # Process is already running.
        if self._process is not None:
            if self._process.isRunning():
                return

        # Reset the watcher.
        self._is_watching = False

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
            f.write("# GROMACS was run with the following command:\n")
            f.write("%s\n" % self._command)

        # Start the timer.
        self._timer = _timeit.default_timer()

        # Start the simulation.
        self._process = _Sire.Base.Process.run(self._exe, args,
            "%s.out" % self._name, "%s.out" % self._name)

        # For historical reasons (console message aggregation with MPI), Gromacs
        # writes the majority of its output to stderr. For user convenience, we
        # redirect all output to stdout, and place a message in the stderr file
        # to highlight this.
        with open(self._stderr_file, "w") as f:
            f.write("All output has been redirected to the stdout stream!")

        # Change back to the original working directory.
        _os.chdir(dir)

        return self

    def _generate_args(self):
        """Generate the dictionary of command-line arguments."""

        # Clear the existing arguments.
        self.clearArgs()

        # Add the default arguments.
        self.setArg("mdrun", True)          # Use mdrun.
        self.setArg("-v", True)             # Verbose output.
        self.setArg("-deffnm", self._name)  # Output file prefix.

    def _generate_binary_run_file(self):
        """Use grommp to generate the binary run input file."""

        # Use grompp to generate the portable binary run input file.
        command = "%s grompp -f %s -po %s.out.mdp -c %s -p %s -o %s" \
            % (self._exe, self._config_file, self._config_file.split(".")[0],
                self._gro_file, self._top_file, self._tpr_file)

        # Run the command.
        proc = _subprocess.run(command, shell=True,
            stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)

        # Get the data prefix.
        if proc.returncode != 0:
            raise RuntimeError("Unable to generate GROMACS binary run input file")

    def setConfig(self, config):
        """Set the list of configuration file strings."""

        # Call the base class method.
        super().setConfig(config)

        # Use grompp to generate the portable binary run input file.
        self._generate_binary_run_file()

    def addToConfig(self, config):
        """Add a string to the configuration list."""

        # Call the base class method.
        super().addToConfig(config)

        # Use grompp to generate the portable binary run input file.
        self._generate_binary_run_file()

    def resetConfig(self):
        """Reset the configuration parameters."""
        self._generate_config()

        # Use grompp to generate the portable binary run input file.
        self._generate_binary_run_file()
