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
Functionality for equilibration protocols.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

from ._protocol import Protocol as _Protocol

import BioSimSpace.Types as _Types

import math as _math
import pytest as _pytest
import warnings as _warnings

__all__ = ["Equilibration"]

class Equilibration(_Protocol):
    """A class for storing equilibration protocols."""

    def __init__(self,
                 timestep=_Types.Time(2, "femtosecond"),
                 runtime=_Types.Time(0.2, "nanoseconds"),
                 temperature_start=_Types.Temperature(300, "kelvin"),
                 temperature_end=_Types.Temperature(300, "kelvin"),
                 temperature=None,
                 frames=20,
                 restrain_backbone=False
                ):
        """Constructor.

           Keyword arguments
           -----------------

           timestep : BioSimSpace.Types.Time
               The integration timestep (in femtoseconds).

           runtime : BioSimSpace.Types.Time
               The running time (in nanoseconds).

           temperature_start : BioSimSpace.Types.Temperature
               The starting temperature (in Kelvin).

           temperature_end : BioSimSpace.Types.Temperature
               The final temperature (in Kelvin).

           temperature : BioSimSpace.Types.Temperature
               The equilibration temperature (in Kelvin). This takes
               precedence of over the other temperatures, i.e. to run
               at fixed temperature.

           frames : int
               The number of trajectory frames to record.

           restrain_backbone : bool
               Whether the atoms in the backbone are fixed.
        """

        # Set the time step.
        self.setTimeStep(timestep)

        # Set the running time.
        self.setRunTime(runtime)

        # Constant temperature equilibration.
        if temperature is not None:
            self.setStartTemperature(temperature)
            self.setEndTemperature(temperature)
            self._is_const_temp = True

        # Heating / cooling simulation.
        else:
            self._is_const_temp = False

            # Set the start temperature.
            self.setStartTemperature(temperature_start)

            # Set the final temperature.
            self.setEndTemperature(temperature_end)

            # Constant temperature simulation.
            if self._temperature_start == self._temperature_end:
                self._is_const_temp = True

        # Set the number of trajectory frames.
        self.setFrames(frames)

        # Set the backbone restraint.
        self.setRestraint(restrain_backbone)

    def __str__(self):
        """Return a human readable string representation of the object."""
        return ("<BioSimSpace.Protocol.Equilibration: timestep=%s, runtime=%s, "
                "temperature_start=%s, temperature_end=%s, frames=%d, restrain_backbone=%r>"
               ) % (self._timestep, self._runtime, self._temperature_start,
                       self._temperature_end, self._frames, self._is_restrained)

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return ("BioSimSpace.Protocol.Equilibration(timestep=%s, runtime=%s, "
                "temperature_start=%s, temperature_end=%s, frames=%d, restrain_backbone=%r)"
               ) % (self._timestep, self._runtime, self._temperature_start,
                       self._temperature_end, self._frames, self._is_restrained)

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


           Positional arguments
           --------------------

           time : BioSimSpace.Types.Time
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


           Positional arguments
           --------------------

           runtime : BioSimSpace.Types.Time
               The simulation run time.
        """
        if type(runtime) is _Types.Time:
            self._runtime = runtime
        else:
            raise TypeError("'runtime' must be of type 'BioSimSpace.Types.Time'")

    def getStartTemperature(self):
        """Return the starting temperature.

           Returns
           -------

           temperature : BioSimSpace.Types.Temperature
               The starting temperature.
        """
        return self._temperature_start

    def setStartTemperature(self, temperature):
        """Set the starting temperature.


           Positional arguments
           --------------------

           temperature : BioSimSpace.Types.Time
               The starting temperature.
        """

        if type(temperature) is _Types.Temperature:
            if temperature.kelvin().magnitude() == _pytest.approx(0):
                temperature._magnitude = 0.01
            self._temperature_start = temperature
        else:
            raise TypeError("'temperature_start' must be of type 'BioSimSpace.Types.Temperature'")

    def getEndTemperature(self):
        """Return the final temperature.

           Returns
           -------

           temperature : BioSimSpace.Types.Temperature
               The final temperature.
        """
        return self._temperature_end

    def setEndTemperature(self, temperature):
        """Set the final temperature.


           Positional arguments
           --------------------

           temperature : BioSimSpace.Types.Time
               The final temperature.
        """
        if type(temperature) is _Types.Temperature:
            if temperature.kelvin().magnitude() == _pytest.approx(0):
                temperature._magnitude = 0.01
            self._temperature_end = temperature
        else:
            raise TypeError("'temperature_start' must be of type 'BioSimSpace.Types.Temperature'")

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


           Positional arguments
           --------------------

           frames : int
               The number of trajectory frames.
        """
        if type(frames) is not int:
            raise TypeError("'frames' must be of type 'int'")

        if frames <= 0:
            _warnings.warn("The number of frames must be positive. Using default (20).")
            self._frames = 20
        else:
            self._frames = _math.ceil(frames)

    def isRestrained(self):
        """Return whether the backbone is restrained.

           Returns
           -------

           is_restrained : bool
               Whether the backbone is restrained.
        """
        return self._is_restrained

    def setRestraint(self, restrain_backbone):
        """Set the backbone restraint.


           Positional arguments
           --------------------

           is_restrained : bool
               Whether the backbone is restrained.
        """
        if type(restrain_backbone) is bool:
            self._is_restrained = restrain_backbone
        else:
            _warnings.warn("Non-boolean backbone restraint flag. Defaulting to no restraint!")
            self._is_restrained = False

    def isConstantTemp(self):
        """Return whether the protocol has a constant temperature.

           Returns
           -------

           is_const_temp : bool
               Whether the temperature is fixed.
        """
        return self._temperature_start == self._temperature_end
