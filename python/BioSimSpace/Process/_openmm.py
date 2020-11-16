######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2020
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
__email_ = "lester.hedges@gmail.com"

__all__ = ["OpenMM"]

import math as _math
import os as _os
import pygtail as _pygtail
import sys as _sys
import timeit as _timeit
import warnings as _warnings

from Sire import Base as _SireBase
from Sire import IO as _SireIO

from BioSimSpace._Exceptions import MissingSoftwareError as _MissingSoftwareError
from BioSimSpace._SireWrappers import System as _System
from BioSimSpace.Trajectory import Trajectory as _Trajectory
from BioSimSpace.Types._type import Type as _Type

from BioSimSpace import IO as _IO
from BioSimSpace import Protocol as _Protocol
from BioSimSpace import Types as _Types
from BioSimSpace import Units as _Units
from BioSimSpace import _Utils as _Utils

from . import _process

class OpenMM(_process.Process):
    """A class for running simulations using OpenMM."""

    def __init__(self, system, protocol, exe=None, name="openmm",
            work_dir=None, seed=None, property_map={}):
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

        # Initialise the stdout dictionary and title header.
        self._stdout_dict = _process._MultiDict()

        # Store the name of the OpenMM log file.
        self._log_file = "%s/%s.log" % (self._work_dir, name)

        # Initialise the log file separator.
        self._record_separator = None

        # Initialise a dictionary to map log file records to their column order.
        self._record_mapping = {}

        # The names of the input files. We choose to use GROMACS files but could
        # equally work with AMBER files.
        self._gro_file = "%s/%s.gro" % (self._work_dir, name)
        self._top_file = "%s/%s.top" % (self._work_dir, name)

        # The name of the trajectory file.
        self._traj_file = "%s/%s.dcd" % (self._work_dir, name)

        # Set the path for the OpenMM Python script. (We use the concept of a
        # config file for consistency with other Process classes.)
        self._config_file = "%s/%s.py" % (self._work_dir, name)

        # Create the list of input files.
        self._input_files = [self._config_file, self._gro_file, self._top_file]

        # Now set up the working directory for the process.
        self._setup()

    def _setup(self):
        """Setup the input files and working directory ready for simulation."""

        # Create the input files...

        # GRO87 file.
        gro = _SireIO.Gro87(self._system._sire_object, self._property_map)
        gro.writeToFile(self._gro_file)

        # TOP file.
        top = _SireIO.GroTop(self._system._sire_object, self._property_map)
        top.writeToFile(self._top_file)

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

    def _generate_args(self):
        """Generate the dictionary of command-line arguments."""

        # Clear the existing arguments.
        self.clearArgs()

    def addToConfig(self, config):
        """Add a string to the configuration list.

           Parameters
           ----------

           config : str, [str]
               A configuration string, a list of configuration strings, or a
               path to a configuration file.
        """

        # Call the base class method.
        super().addToConfig(config)

    def resetConfig(self):
        """Reset the configuration parameters."""
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
                f.write("# GROMACS was run with the following command:\n")
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

        if type(self._protocol) is _Protocol.Minimisation or \
          (type(self._protocol) is _Protocol.Custom and _is_minimisation(self.getConfig())):
            # Create the name of the restart GRO file.
            restart = "%s/%s.gro" % (self._work_dir, self._name)

            # Check that the file exists.
            if _os.path.isfile(restart):
                # Read the molecular system.
                new_system = _System(_SireIO.MoleculeParser.read([restart, self._top_file], self._property_map))

                # Copy the new coordinates back into the original system.
                old_system = self._system.copy()
                old_system._updateCoordinates(new_system,
                                              self._property_map,
                                              self._property_map)

                # Update the periodic box information in the original system.
                if "space" in new_system._sire_object.propertyKeys():
                    box = new_system._sire_object.property("space")
                    old_system._sire_object.setProperty(self._property_map.get("space", "space"), box)

                return old_system

            else:
                return None

        else:
            # Grab the most recent frame from the trajectory file.
            return self._getFrame(self.getTime())

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

        try:
            # Locate the trajectory file.
            traj_file = self._find_trajectory_file()

            if traj_file is None:
                return None
            else:
                self._traj_file = traj_file

            return _Trajectory(process=self)

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
        """Get the time (in nanoseconds).

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
        """Get the current time (in nanoseconds).

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
                records = line[1:].split(self._record_separator)

                # Map each record string to its position in the array (column order).
                for idx, record in enumerate(records):
                    # Strip the extra quotes from the record.
                    record = record[1:-1]

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
                    return [float(x) * unit for x in self._stdout_dict[key]]

            except KeyError:
                return None

        # Return the most recent dictionary value.
        else:
            try:
                if unit is None:
                    return float(self._stdout_dict[key][-1])
                else:
                    return float(self._stdout_dict[key][-1]) * unit

            except KeyError:
                return None
