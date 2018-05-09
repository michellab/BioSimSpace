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
                self._exe = exe

            # Search system PATH.
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
        command = "%s -h -h 2>&1 | grep 'Data prefix' | awk -F ':' '{print $2}'" % self._exe

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
        self._config_file = "%s/%s.gromacs" % (self._work_dir, name)

        # Create the list of input files.
        self._input_files = [self._config_file, self._gro_file, self._top_file]

        # Now set up the working directory for the process.
        self._setup()

    def _setup(self):
        """Setup the input files and working directory ready for simulation."""

        # Create the input files...

        # To read a full GROMACS system using Sire...
        # s = _Sire.IO.MoleculeParser.read(files, {"GROMACS_PATH":gromacs_path})

        # GRO87 file.
        gro = _Sire.IO.Gro87(self._system)
        gro.writeToFile(self._gro_file)

        # TOP file.
        top = _Sire.IO.GroTop(self._system)
        top.writeToFile(self._top_file)

        # Generate the GROMACS configuration file.
        # Skip if the user has passed a custom config.
        if type(self._protocol) is _Protocol.Custom:
            self.setConfig(self._protocol.getConfig())
        else:
            self._generate_config()
        self.writeConfig(self._config_file)

        # Return the list of input files.
        return self._input_files

    def _generate_config(self):
        """Generate GROMACS configuration file strings."""

        # Clear the existing configuration list.
        self._config = []
