######################################################################
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
Functionality for production protocols.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

from ._protocol import Protocol as _Protocol

import BioSimSpace.Types as _Types

import math as _math

__all__ = ["Production"]

class Production(_Protocol):
    """A class for storing production protocols."""

    # A list of allowed thermodynamic ensembles.
    # Update as support is added.
    _ensembles = ["NVT", "NPT"]

    def __init__(self,
                 timestep=_Types.Time(2, "femtosecond"),
                 runtime=_Types.Time(1, "nanosecond"),
                 temperature=_Types.Temperature(300, "kelvin"),
                 frames=20,
                 ensemble="NPT",
                 first_step=0,
                 restart=False
                ):
        """Constructor.

           Parameters
           ----------

           timestep : BioSimSpace.Types.Time
               The integration timestep.

           runtime : BioSimSpace.Types.Time
               The running time.

           temperature : BioSimSpace.Types.Temperature
               The temperature.

           frames : int
               The number of trajectory frames to record.

           ensemble : str
               The thermodynamic ensemble.

           first_step : int
               The initial time step (for restart simulations).

           restart : bool
               Whether this is a continuation of a previous simulation.
        """

        # Set the time step.
        self.setTimeStep(timestep)

        # Set the runtime.
        self.setRunTime(runtime)

        # Set the system temperature.
        self.setTemperature(temperature)

        # Set the number of trajectory frames.
        self.setFrames(frames)

        # Set the thermodynamic ensemble.
        self.setEnsemble(ensemble)

        # Set the restart flag.
        self.setRestart(restart)

        # Set the first time step.
        self.setFirstStep(first_step)

    def __str__(self):
        """Return a human readable string representation of the object."""
        return ("<BioSimSpace.Protocol.Production: timestep=%s, runtime=%s, "
                "temperature=%s, frames=%d, ensemble=%r, first_step=%d, restart=%r>"
               ) % (self._timestep, self._runtime, self._temperature, self._frames,
                       self._ensemble, self._first_step, self._restart)

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return ("BioSimSpace.Protocol.Production(timestep=%s, runtime=%s, "
                "temperature=%s, frames=%d, ensemble=%r, first_step=%d, restart=%r)"
               ) % (self._timestep, self._runtime, self._temperature, self._frames,
                       self._ensemble, self._first_step, self._restart)

    def getTimeStep(self):
        """Return the time step.

           Returns
           -------

           time : BioSimSpace.Types.Time
               The integration time step.
        """
        return self._timestep

    def setTimeStep(self, timestep):
        """Set the time step.

           Parameters
           ----------

           timestep : BioSimSpace.Types.Time
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

           time : BioSimSpace.Types.Time
               The simulation run time.
        """
        return self._runtime

    def setRunTime(self, runtime):
        """Set the running time.

           Parameters
           ----------

           runtime : BioSimSpace.Types.Time
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

           temperature : BioSimSpace.Types.Temperature
               The simulation temperature.
        """
        return self._temperature

    def setTemperature(self, temperature):
        """Set the temperature.

           Parameters
           ----------

           temperature : BioSimSpace.Types.Time
               The simulation temperature.
        """
        if type(temperature) is _Types.Temperature:
            self._temperature = temperature
        else:
            raise TypeError("'temperature' must be of type 'BioSimSpace.Types.Temperature'")

    def getFrames(self):
        """Return the number of frames.

           Returns
           -------

           frames : int
               The number of trajectory frames.
        """
        return self._frames

    def setFrames(self, frames):
        """Set the number of frames.

           Parameters
           ----------

           frames : int
               The number of trajectory frames.
        """
        if type(frames) is not int:
            raise TypeError("'frames' must be of type 'int'")

        if frames <= 0:
            warn("The number of frames must be positive. Using default (20).")
            self._frames = 20
        else:
            self._frames = _math.ceil(frames)

    def getEnsemble(self):
        """Return the thermodynamic ensemble.

           Returns
           -------

           ensemble : str
               The thermodynamic ensemble.
        """
        return self._ensemble

    def setEnsemble(self, ensemble):
        """Set the thermodynamic ensemble.

           Parameters
           ----------

           ensemble : str
               The thermodynamic ensemble.
        """
        if ensemble.replace(" ", "").upper() not in self._ensembles:
            warn("Unsupported thermodynamic ensemble. Using default ('NPT').")
            self._ensemble = "NPT"
        else:
            self._ensemble = ensemble.replace(" ", "").upper()

    def getFirstStep(self):
        """Return the first time step.

           Returns
           -------

           step : int
               The first time step.
        """
        return self._first_step

    def setFirstStep(self, first_step):
        """Set the initial time step.

           Parameters
           ----------

           step : int
               The first time step.
        """
        if type(first_step) is not int:
            raise TypeError("'first_step' must be of type 'int'")

        if first_step < 0:
            warn("The initial time step must be positive. Using default (0).")
            self._first_step = 0
        else:
            self._first_step = _math.ceil(first_step)

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
