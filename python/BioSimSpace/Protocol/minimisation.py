"""
@package biosimspace
@author  Lester Hedges
@brief   A class for storing minimisation protocols.
"""

from .protocol import Protocol, ProtocolType

from warnings import warn

class Minimisation(Protocol):
    """A class for storing minimisation protocols."""

    def __init__(self, method="conjugate-gradient", steps=10000, gas_phase=False):
        """Constructor.

           Keyword arguments:

           steps     -- The maximum number of steps to perform.
           gas_phase -- Whether this a gas phase simulation.
        """

        # Call the base class constructor.
        super().__init__(ProtocolType.MINIMISATION, gas_phase)

        # Set the number of steps.
        self._steps = steps

    @property
    def steps(self):
        """Return the maximum number of steps."""
        return self._steps

    @steps.setter
    def steps(self, steps):
        """Set the maximum number of steps."""

        if steps <= 0:
            warn("Number of steps must be greater than zero. Using default (10000).")
            self._steps = 10000

        else:
            self._steps = steps
