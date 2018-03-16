"""
@package biosimspace
@author  Lester Hedges
@brief   A class for storing production protocols.
"""

from .protocol import Protocol, ProtocolType

from math import ceil

# A list of allowed thermodynamic ensembles.
# Update as support is added.
_ensembles = ['NVT', 'NPT']

class Production(Protocol):
    """A class for storing production protocols."""

    def __init__(self, timestep=2, runtime=1, temperature=300, frames=20,
            ensemble="NPT", first_step=0, restart=False, gas_phase=False):
        """Constructor.

           Keyword arguments:

           timestep    -- The integration timestep (in femtoseconds).
           runtime     -- The running time (in nanoseconds).
           temperature -- The temperature (in Kelvin).
           frames      -- The number of trajectory frames to record.
           ensemble    -- The thermodynamic ensemble.
           first_step  -- The initial time step (for restart simulations).
           restart     -- Whether this is a continuation of a previous simulation.
           gas_phase   -- Whether this a gas phase simulation.
        """

        # Call the base class constructor.
        super().__init__(ProtocolType.PRODUCTION, gas_phase)

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
        step = self.setFirstStep(first_step)

    def getTimeStep(self):
        """Return the time step."""
        return self._timestep

    def setTimeStep(self, timestep):
        """Set the time step."""

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

        if runtime <= 0:
            warn("The running time must be positive. Using default (1 ns).")
            self._runtime = 1

        else:
            self._runtime = runtime

    def getTemperature(self):
        """Return temperature."""
        return self._temperature

    def setTemperature(self, temperature):
        """Set the temperature."""

        if temperature <= 0:
            warn("Temperature must be positive. Using default (300 K).")
            self._temperature = 300

        else:
            self._temperature = temperature

    def getFrames(self):
        """Return the number of frames."""
        return self._frames

    def setFrames(self, frames):
        """Set the number of frames."""

        if frames <= 0:
            warn("The number of frames must be positive. Using default (20).")
            self._frames = 20

        else:
            self._frames = ceil(frames)

    def getEnsemble(self):
        """Return the thermodynamic ensemble."""
        return self._ensemble

    def setEnsemble(self, ensemble):
        """Set the thermodynamic ensemble."""

        if ensemble.strip().upper() not in _ensembles:
            warn("Unsupported thermodynamic ensemble. Using default ('NPT').")
            self._ensemble = 'NPT'

        else:
            self._ensemble = ensemble.strip().upper()

    def getFirstStep(self):
        """Return the first time step."""
        return self._first_step

    def setFirstStep(self, first_step):
        """Set the initial time step."""

        if first_step < 0:
            warn("The initial time step must be positive. Using default (0).")
            self._first_step = 0

        else:
            self._first_step = ceil(first_step)

    def isRestart(self):
        """Return whether this restart simulation."""
        return self._restart

    def setRestart(self, restart):
        """Set the restart flag."""

        if type(restart) is bool:
            self._restart = restart

        else:
            warn("Non-boolean restart flag. Defaulting to False!")
            self._restart = False
