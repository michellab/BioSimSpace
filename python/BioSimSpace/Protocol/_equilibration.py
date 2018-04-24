######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright 2017-2018
#
# Authors: Lester Hedges

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
Functionality for equilibration protocols.
Author: Lester Hedges
"""

from ._protocol import Protocol

from math import ceil
from pytest import approx
from warnings import warn

__all__ = ["Equilibration"]

class Equilibration(Protocol):
    """A class for storing equilibration protocols."""

    def __init__(self, timestep=2, runtime=0.2, temperature_start=300,
            temperature_end=None, frames=20, restrain_backbone=False):
        """Constructor.

           Keyword arguments:

           timestep          -- The integration timestep (in femtoseconds).
           runtime           -- The running time (in nanoseconds).
           temperature_start -- The starting temperature (in Kelvin).
           temperature_end   -- The final temperature (in Kelvin).
           frames            -- The number of trajectory frames to record.
           restrain_backbone -- Whether the atoms in the backbone are fixed.
        """

        # Set the time step.
        self.setTimeStep(timestep)

        # Set the running time.
        self.setRunTime(runtime)

        # Set the start temperature.
        self.setStartTemperature(temperature_start)

        # Set the final temperature.
        if temperature_end is not None:
            self.setEndTemperature(temperature_end)
            self._is_const_temp = False

            # Start and end temperature is the same.
            if (self._temperature_start == self._temperature_end):
                warn("Start and end temperatures are the same!")
                self._is_const_temp = True

        # Constant temperature simulation.
        else:
            self._temperature_end = self._temperature_start
            self._is_const_temp = True

        # Set the number of trajectory frames.
        self.setFrames(frames)

        # Set the backbone restraint.
        self.setRestraint(restrain_backbone)

    def __str__(self):
        """Return a human readable string representation of the object."""
        return ("<BioSimSpace.Protocol.Equilibration: timestep=%.2f, runtime=%.2f, "
                "temperature_start=%.2f, temperature_end=%.2f, frames=%d, restrain_backbone=%r>"
               ) % (self._timestep, self._runtime, self._temperature_start,
                       self._temperature_end, self._frames, self._is_restrained)

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return ("BioSimSpace.Protocol.Equilibration(timestep=%.2f, runtime=%.2f, "
                "temperature_start=%.2f, temperature_end=%.2f, frames=%d, restrain_backbone=%r)"
               ) % (self._timestep, self._runtime, self._temperature_start,
                       self._temperature_end, self._frames, self._is_restrained)

    def getTimeStep(self):
        """Return the time step."""
        return self._timestep

    def setTimeStep(self, timestep):
        """Set the time step."""

        if type(timestep) is int:
            timestep = float(timestep)

        if type(timestep) is not float:
            raise TypeError("'timestep' must be of type 'float'")

        if timestep <= 0:
            warn("The time step must be positive. Using default (2 fs).")
            self._timestep = 2

        else:
            self._timestep = timestep

    def getRunTime(self):
        """Return the running time."""
        return self._runtime

    def setRunTime(self, runtime):
        """Set the running time."""

        if type(runtime) is int:
            runtime = float(runtime)

        if type(runtime) is not float:
            raise TypeError("'runtime' must be of type 'float'")

        if runtime <= 0:
            warn("The running time must be positive. Using default (0.2 ns).")
            self._runtime = 0.2

        else:
            self._runtime = runtime

    def getStartTemperature(self):
        """Return the starting temperature."""
        return self._temperature_start

    def setStartTemperature(self, temperature):
        """Set the starting temperature."""

        if type(temperature) is int:
            temperature = float(temperature)

        if type(temperature) is not float:
            raise TypeError("'temperature' must be of type 'float'")

        if temperature < 0:
            warn("Starting temperature must be positive. Using default (300 K).")
            self._temperature_start = 300

        elif temperature == approx(0):
            self._temperature_start = 0.01

        else:
            self._temperature_start = temperature

    def getEndTemperature(self):
        """Return the final temperature."""
        return self._temperature_end

    def setEndTemperature(self, temperature):
        """Set the final temperature."""

        if type(temperature) is int:
            temperature = float(temperature)

        if type(temperature) is not float:
            raise TypeError("'temperature' must be of type 'float'")

        if temperature < 0:
            warn("Final temperature must be positive. Using default (300 K).")
            self._temperature_end = 300

        elif temperature == approx(0):
            self._temperature_end = 0.01

        else:
            self._temperature_end = temperature

    def getFrames(self):
        """Return the number of frames."""
        return self._frames

    def setFrames(self, frames):
        """Set the number of frames."""

        if type(frames) is not int:
            raise TypeError("'frames' must be of type 'int'")

        if frames <= 0:
            warn("The number of frames must be positive. Using default (20).")
            self._frames = 20

        else:
            self._frames = ceil(frames)

    def isRestrained(self):
        """Return whether the backbone is restrained."""
        return self._is_restrained

    def setRestraint(self, restrain_backbone):
        """Set the backbone restraint."""

        if type(restrain_backbone) is bool:
            self._is_restrained = restrain_backbone

        else:
            warn("Non-boolean backbone restraint flag. Defaulting to no restraint!")
            self._is_restrained = False

    def isConstantTemp(self):
        """Return whether the protocol has a constant temperature."""
        return self._is_const_temp
