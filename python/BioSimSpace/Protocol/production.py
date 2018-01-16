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

    def __init__(self, runtime=1, temperature=300, frames=20, ensemble="NPT", first_step=0):
        """Constructor.

           Keyword arguments:

           runtime     -- The running time (in nanoseconds).
           temperature -- The temperature (in Kelvin).
           frames      -- The number of trajectory frames to record.
           ensemble    -- The thermodynamic ensemble.
           first_step  -- The initial time step (for restart simulations).
        """

        # Call the base class constructor.
        super().__init__(ProtocolType.PRODUCTION)

        # Set the runtime.
        self._runtime = runtime

        # Set the system temperature.
        self._temperature = temperature

        # Set the number of trajectory frames.
        self._frames = frames

        # Set the thermodynamic ensemble.
        self._ensemble = ensemble

        # Set the first time step.
        self._first_step = first_step

    def type(self):
        """Return the protocol type."""
        return self._type
    @property
    def runtime(self):
        """Return the running time."""
        return self._runtime

    @runtime.setter
    def runtime(self, runtime):
        """Set the running time."""

        if runtime <= 0:
            warn("The running time must be positive. Using default (1 ns).")
            self._runtime = 1

        else:
            self._runtime = runtime

    @property
    def temperature(self):
        """Return temperature."""
        return self._temperature

    @temperature.setter
    def temperature(self, temperature):
        """Set the temperature."""

        if temperature <= 0:
            warn("Temperature must be positive. Using default (300 K).")
            self._temperature = 300

        else:
            self._temperature = temperature

    @property
    def frames(self):
        """Return the number of frames."""
        return self._frames

    @frames.setter
    def frames(self, frames):
        """Set the number of frames."""

        if frames <= 0:
            warn("The number of frames must be positive. Using default (20).")
            self._frames = 20

        else:
            self._frames = ceil(frames)

    @property
    def ensemble(self):
        """Return the thermodynamic ensemble."""
        return self._ensemble

    @ensemble.setter
    def ensemble(self, ensemble):
        """Set the thermodynamic ensemble."""

        if ensemble.upper() not in _ensembles:
            warn("Unsupported thermodynamic ensemble. Using default ('NPT').")
            self._ensemble = 'NPT'

        else:
            self._ensemble = ensemble
    @property
    def first_step(self):
        """Return the first time step."""
        return self._first_step

    @first_step.setter
    def first_step(self, first_step):
        """Set the initial time step."""

        if first_step <= 0:
            warn("The initial time step must be positive. Using default (0).")
            self._first_step = 0

        else:
            self._first_step = ceil(first_step)
