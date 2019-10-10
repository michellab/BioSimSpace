#####################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2019
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
Functionality for metadynamics protocols.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["Metadynamics"]

import math as _math
import os as _os

from BioSimSpace import Types as _Types
from BioSimSpace.Metadynamics import CollectiveVariable as _CollectiveVariable

from ._protocol import Protocol as _Protocol

# Store the collective variable base type.
_colvar_type = _CollectiveVariable._collective_variable.CollectiveVariable

class Metadynamics(_Protocol):
    """A class for storing metadynamics protocols."""

    def __init__(self,
                 collective_variable,
                 timestep=_Types.Time(2, "femtosecond"),
                 runtime=_Types.Time(1, "nanosecond"),
                 temperature=_Types.Temperature(300, "kelvin"),
                 pressure=_Types.Pressure(1, "atmosphere"),
                 hill_height=_Types.Energy(1, "kj per mol"),
                 hill_frequency=1000,
                 bias_factor=None,
                 hills_file=None,
                 grid_file=None,
                 colvar_file=None
                ):
        """Constructor.

           Parameters
           ----------

           collective_variable : :class:`CollectiveVariable <BioSimSpace.Metadynamics.CollectiveVariable>`, \
                                [:class:`CollectiveVariable <BioSimSpace.Metadynamics.CollectiveVariable>`]
               The collective variable (or variables) for the simulation.

           timestep : :class:`Time <BioSimSpace.Types.Time>`
               The integration timestep.

           runtime : :class:`Time <BioSimSpace.Types.Time>`
               The running time.

           temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
               The temperature.

           pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
               The pressure. Pass pressure=None to use the NVT ensemble.

           hill_height : :class:`Energy <BioSimSpace.Types.Energy>`
               The height of the Gaussian hills.

           hill_frequency : int
               The frequency at which hills are deposited.

           bias_factor : float
               The bias factor for well tempered metadynamics.

           hills_file : str
               The path to a HILLS file from a previous simulation. This can
               be used to restart in order to contiune sampling. The information
               in the file must be consistent with the 'collective_variable'
               argument.

           grid_file : str
               The path to a GRID file from a previous simulation. This can
               be used to restart in order to continue sampling. The information
               in the file must be consistent with the 'collective_variable'
               argument.

           colvar_file : str
               The path to a COLVAR file from a previous simulation. The
               information in the file must be consistent with the
               'collective_variable' argument.
        """

        # Call the base class constructor.
        super().__init__()

        # Whether this is a newly created object.
        self._is_new_object = True

        # Set the collective variable.
        self.setCollectiveVariable(collective_variable)

        # Set the time step.
        self.setTimeStep(timestep)

        # Set the runtime.
        self.setRunTime(runtime)

        # Set the system temperature.
        self.setTemperature(temperature)

        # Set the system pressure.
        if pressure is not None:
            self.setPressure(pressure)
        else:
            self._pressure = None

        # Set the hill parameters: height, frequency.
        self.setHillHeight(hill_height)
        self.setHillFrequency(hill_frequency)

        # Set the bias factor for well tempered metadynamics.
        if bias_factor is not None:
            self.setBiasFactor(bias_factor)
        else:
            self._bias_factor = None

        # Set the restart files.
        if hills_file is not None:
            self.setHillsFile(hills_file)
        else:
            self._hills_file = None
        if grid_file is not None:
            self.setGridFile(grid_file)
        else:
            self._grid_file = None
        if colvar_file is not None:
            self.setColvarFile(colvar_file)
        else:
            self._colvar_file = None

        # Flag that the object has been created.
        self._is_new_object = False

    def __str__(self):
        """Return a human readable string representation of the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            string = "<BioSimSpace.Protocol.Metadynamics: "
            string += "collective_variable=%s" % self._collective_variable
            string += ", timestep=%s" % self._timestep
            string += ", runtime=%s" % self._runtime
            string += ", temperature=%s" % self._temperature
            if self._pressure is not None:
                string += ", pressure=%s" % self._pressure
            string += ", hill_height=%s" % self._hill_height
            string += ", hill_frequency=%s" % self._hill_frequency
            if self._bias_factor is not None:
                string += ", bias_factor=%s" % self._bias_factor
            if self._hills_file is not None:
                string += ", hills_file=%r" % self._hills_file
            if self._grid_file is not None:
                string += ", grid_file=%r" % self._grid_file
            if self._colvar_file is not None:
                string += ", colvar_file=%r" % self._colvar_file
            string += ">"

            return string

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return self.__str__()

    def getCollectiveVariable(self):
        """Return the collective variable (or variables).

           Returns
           -------

           collective_variable : [:class:`CollectiveVariable <BioSimSpace.Metadynamics.CollectiveVariable>`]
               The collective variable (or variables) for the simulation.
        """
        return self._collective_variable.copy()

    def setCollectiveVariable(self, collective_variable):
        """Set the collective variable (or variables).

           Parameters
           ----------

           collective_variable : :class:`CollectiveVariable <BioSimSpace.Metadynamics.CollectiveVariable>`, \
                                [:class:`CollectiveVariable <BioSimSpace.Metadynamics.CollectiveVariable>`]
               The collective variable (or variables) for the simulation.
        """

        # A single collective variable.
        if isinstance(collective_variable, _colvar_type):
            self._collective_variable = [collective_variable]
            return

        # Convert tuple to list.
        if type(collective_variable) is tuple:
            collective_variable = tuple(collective_variable)

        if type(collective_variable) is list:
            if not all(isinstance(x, _colvar_type) for x in collective_variable):
                raise TypeError("'collective_variable' must all be of type "
                                "'BioSimSpace.Metadynamics.CollectiveVariable'")
        else:
            raise TypeError("'collective_variable' must be of type "
                            "'BioSimSpace.Metadynamics.CollectiveVariable' "
                            "or a list of 'BioSimSpace.Metadynamics.CollectiveVariable' types.")

        # Make sure all of the collective variables are consistent. If any have
        # a grid set, then so must all other variables.

        num_grid = 0

        for colvar in collective_variable:
            if colvar.getGrid() is not None:
                num_grid += 1

        if num_grid > 0 and num_grid != len(collective_variable):
            raise ValueError("If a 'grid' is desired, then all collective "
                             "variables must define one.")

        self._collective_variable = collective_variable

        # If the object has already been created, then check that other member
        # data is consistent.
        if not self._is_new_object:
            self.setHillHeight(self._hill_height)

    def getTimeStep(self):
        """Return the time step.

           Returns
           -------

           timestep : :class:`Time <BioSimSpace.Types.Time>`
               The integration time step.
        """
        return self._timestep

    def setTimeStep(self, timestep):
        """Set the time step.

           Parameters
           ----------

           timestep : :class:`Time <BioSimSpace.Types.Time>`
               The integration time step.
        """
        if type(timestep) is _Types.Time:
            self._timestep = timestep
        else:
            raise TypeError("'timestep' must be of type 'BioSimSpace.Types.Time'")

    def getRunTime(self):
        """Return the running time.

           Returns
           -------

           runtime : :class:`Time <BioSimSpace.Types.Time>`
               The simulation run time.
        """
        return self._runtime

    def setRunTime(self, runtime):
        """Set the running time.

           Parameters
           ----------

           runtime : :class:`Time <BioSimSpace.Types.Time>`
               The simulation run time.
        """
        if type(runtime) is _Types.Time:
            self._runtime = runtime
        else:
            raise TypeError("'runtime' must be of type 'BioSimSpace.Types.Time'")

    def getTemperature(self):
        """Return temperature.

           Returns
           -------

           temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
               The simulation temperature.
        """
        return self._temperature

    def setTemperature(self, temperature):
        """Set the temperature.

           Parameters
           ----------

           temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
               The simulation temperature.
        """
        if type(temperature) is _Types.Temperature:
            self._temperature = temperature
        else:
            raise TypeError("'temperature' must be of type 'BioSimSpace.Types.Temperature'")

    def getPressure(self):
        """Return the pressure.

           Returns
           -------

           pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
               The pressure.
        """
        return self._pressure

    def setPressure(self, pressure):
        """Set the pressure.

           Parameters
           ----------

           pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
               The pressure.
        """
        if type(pressure) is _Types.Pressure:
            self._pressure = pressure
        else:
            raise TypeError("'pressure' must be of type 'BioSimSpace.Types.Pressure'")

    def getHillHeight(self):
        """Return the height of the Gaussian hills.

           Returns
           -------

           hill_width : :class:`Energy <BioSimSpace.Types.Energy>`
               The height of the Gaussian hills.
        """
        return self._hill_height

    def setHillHeight(self, hill_height):
        """Set the height of the Gaussian hills.

           Parameters
           ----------

           hill_height : :class:`Energy <BioSimSpace.Types.Energy>`
               The hill height.
        """

        if type(hill_height) is not _Types.Energy:
            raise TypeError("'hill_height' must be of type 'BioSimSpace.Types.Energy'")

        # Check that heights is greater than zero.
        if hill_height.magnitude() <= 0:
            raise ValueError("Cannot have 'hill_height' with magnitude <= 0")

        self._hill_height = hill_height.kj_per_mol()

    def getHillFrequency(self):
        """Return the frequency at which Gaussian hills are deposited.

           Returns
           -------

           hill_frequency : int
               The frequency at which hills are deposited.
        """
        return self._hill_frequency

    def setHillFrequency(self, hill_frequency):
        """Set the frequency at which Gaussian hills are deposited.

           Parameters
           ----------

           hill_frequency : int
               The frequency at which hills are deposited.
        """

        try:
            hill_frequency = int(hill_frequency)
        except:
            raise TypeError("'hill_frequency' must be of type 'int'")

        if hill_frequency < 1:
            raise ValueError("'hill_frequency' must be >= 1")

        self._hill_frequency = hill_frequency

    def getBiasFactor(self):
        """Return the bias factor for well tempered metadynamics.

           Returns
           -------

           bias_factor : float
               The bias factor for well tempered metadynamics.
        """
        return self._bias_factor

    def setBiasFactor(self, bias_factor=None):
        """Set the bias factor for well tempered metadynamics.
           Call with no arguments to clear the bias factor.

           Parameters
           ----------

           bias_factor : float
               The bias factor for well tempered metadynamics.
        """

        if bias_factor is None:
            self._bias_factor = None
            return

        try:
            bias_factor = float(bias_factor)
        except:
            raise TypeError("'bias_factor' must be of type 'float'")

        if bias_factor <= 0:
            raise ValueError("'bias_factor' must be > 0")

        self._bias_factor = bias_factor

    def getHillsFile(self):
        """Return the path to the HILLS file.

           Returns
           -------

           hills_file : str
               The path to the HILLS file.
        """
        return self._hills_file

    def setHillsFile(self, hills_file):
        """Set the location of an existing HILLS file.

           Parameters
           ----------

           hills_file : str
               The path to an existing HILLS file.
        """
        if type(hills_file) is not str:
            raise ValueError("'hills_file' must be of type 'str'")

        if not _os.path.isfile(hills_file):
            raise ValueError("'hills_file' doesn't exist: %s" % hills_file)

        # Read the header and make sure that it is consistent with the
        # collective variables.
        with open(hills_file, "r") as file:
            # Read the header record.
            header = file.readline()

            # Split on whitespace to get the records.
            records = header.split()

            # Make sure this is a valid header.
            if records[0] != "#!" or records[1] != "FIELDS" or records[2] != "time":
                raise ValueError("'hills_file' doesn't contain valid header information!")

            # Make sure there are the right number of collective variables.
            num_colvar = 0
            for x in range(3, len(records)):
                # We've incremented past the collective variable records.
                if "sigma" in records[x]:
                    break
                else:
                    num_colvar += 1

            try:
                if num_colvar != len(self._collective_variable):
                    raise ValueError("'hills_file' contains %d collective variable records, "
                                     "should have %d" % (num_colvar, len(collective_variable)))
            except:
                if num_colvar != 1:
                    raise ValueError("'hills_file' contains %d collective variable records, "
                                     "should have 1" % num_colvar)

        self._hills_file = hills_file

    def getGridFile(self):
        """Return the path to the GRID file.

           Returns
           -------

           grid_file : str
               The path to the GRID file.
        """
        return self._grid_file

    def setGridFile(self, grid_file):
        """Set the location of an existing GRID file.

           Parameters
           ----------

           grid_file : str
               The path to an existing GRID file.
        """
        if type(grid_file) is not str:
            raise ValueError("'grid_file' must be of type 'str'")

        if not _os.path.isfile(grid_file):
            raise ValueError("'grid_file' doesn't exist: %s" % grid_file)

        # Read the header and make sure that it is consistent with the
        # collective variables.
        with open(grid_file, "r") as file:
            # Read the header record.
            header = file.readline()

            # Split on whitespace to get the records.
            records = header.split()

            # Make sure this is a valid header.
            if records[0] != "#!" or records[1] != "FIELDS":
                raise ValueError("'grid_file' doesn't contain valid header information!")

            # Make sure there are the right number of collective variables.
            num_colvar = 0
            for x in range(2, len(records)):
                # We've incremented past the collective variable records.
                if "bias" in records[x]:
                    break
                else:
                    num_colvar += 1

            try:
                if num_colvar != len(self._collective_variable):
                    raise ValueError("'grid_file' contains %d collective variable records, "
                                     "should have %d" % (num_colvar, len(collective_variable)))
            except:
                if num_colvar != 1:
                    raise ValueError("'grid_file' contains %d collective variable records, "
                                     "should have 1" % num_colvar)

        self._grid_file = grid_file

    def getColvarFile(self):
        """Return the path to the COLVAR file.

           Returns
           -------

           colvar_file : str
               The path to the COLVAR file.
        """
        return self._colvar_file

    def setColvarFile(self, colvar_file):
        """Set the location of an existing COLVAR file.

           Parameters
           ----------

           colvar_file : str
               The path to an existing COLVAR file.
        """
        if type(colvar_file) is not str:
            raise ValueError("'colvar_file' must be of type 'str'")

        if not _os.path.isfile(colvar_file):
            raise ValueError("'colvar_file' doesn't exist: %s" % colvar_file)

        self._colvar_file = colvar_file
