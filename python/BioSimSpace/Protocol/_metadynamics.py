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

from ._protocol import Protocol as _Protocol

import BioSimSpace.Types as _Types
import BioSimSpace.Metadynamics.CollectiveVariable as _CollectiveVariable

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
                 hill_width=0.1,
                 hill_height=1,
                 hill_frequency=1000,
                 bias_factor=None,
                 restart=False
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

           hill_width : float, [float]
               The width of the Guassian hills.

           hill_height : float
               The height of the Guassian hills.

           hill_frequency : int
               The frequency at which hills are deposited.

           bias_factor : float
               The bias factor for well tempered metadynamics.

           restart : bool
               Whether this is a continuation of a previous simulation.
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

        # Set the hill parameters: width, height, frequency.
        self.setHillWidth(hill_width)
        self.setHillHeight(hill_height)
        self.setHillFrequency(hill_frequency)

        # Set the bias factor for well tempered metadynamics.
        if bias_factor is not None:
            self.setBiasFactor(bias_factor)
        else:
            self._bias_factor = None

        # Set the restart flag.
        self.setRestart(restart)

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
            string += ", temperature=%s" % self._runtime
            if self._pressure is not None:
                string += ", pressure=%s" % self._pressure
            string += ", hill_width=%s" % self._hill_width
            string += ", hill_height=%s" % self._hill_height
            string += ", hill_frequency=%s" % self._hill_frequency
            if self._bias_factor is not None:
                string += ", bias_factor=%s" % self._bias_factor
            string += ", restart=%s>" % self._restart

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
            self.setHillWidth(self._hill_width)
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

    def getHillWidth(self):
        """Return the width of the Guassian hills.

           Returns
           -------

           hill_width : [float]
               The hill width for each collective variable.
        """
        return self._hill_width

    def setHillWidth(self, hill_width):
        """Set the width of the Guassian hills.

           Parameters
           ----------

           hill_width : [float]
               The hill width for each collective variable.
        """

        # Convert tuple to list.
        if type(hill_width) is tuple:
            hill_width = list(hill_width)

        # A single value.
        if type(hill_width) is float or type(hill_width) is int:
            hill_width = [float(hill_width) for x in range(0, len(self._collective_variable))]

        # A list of values.
        elif type(hill_width) is list:
            widths = []
            for width in hill_width:
                try:
                    widths.append(float(width))
                except:
                    raise TypeError("'hill_width' must be a list of 'float' types.")
            hill_width = widths
            # Must be one for each collective variable.
            if len(hill_width) != len(self._collective_variable):
                raise ValueError(("Number of 'hill_width' parameters (%d) doesn't "
                                  "match the number of collective variables (%d)")
                                  % (len(hill_width), len(self._collective_variable)))
        else:
            raise TypeError("'hill_width' must be of type 'float', or a list of 'float' types.")

        # Check that all widths are greater than zero.
        for width in hill_width:
            if width <= 0:
                raise ValueError("Cannot have 'hill_width' <= 0")

        self._hill_width = hill_width

    def getHillHeight(self):
        """Return the height of the Guassian hills.

           Returns
           -------

           hill_height : float
               The hill height.
        """
        return self._hill_height

    def setHillHeight(self, hill_height):
        """Set the height of the Guassian hills.

           Parameters
           ----------

           hill_height : float
               The hill height.
        """

        if type(hill_height) is float or type(hill_height) is int:
            hill_height = float(hill_height)
        else:
            raise TypeError("'hill_height' must be of type 'float'")

        # Check that heights is greater than zero.
        if hill_height <= 0:
            raise ValueError("Cannot have 'hill_height' <= 0")

        self._hill_height = hill_height

    def getHillFrequency(self):
        """Return the frequency at which Guassian hills are deposited.

           Returns
           -------

           hill_frequency : int
               The frequency at which hills are deposited.
        """
        return self._hill_frequency

    def setHillFrequency(self, hill_frequency):
        """Set the frequency at which Guassian hills are deposited.

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

    def isRestart(self):
        """Return whether this restart simulation.

           Returns
           -------

           is_restart : bool
               Whether this is a restart simulation.
        """
        return self._restart

    def setRestart(self, restart):
        """Set the restart flag.

           Parameters
           ----------

           restart : bool
               Whether this is a restart simulation.
        """
        if type(restart) is bool:
            self._restart = restart
        else:
            warn("Non-boolean restart flag. Defaulting to False!")
            self._restart = False
