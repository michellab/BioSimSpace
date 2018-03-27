"""
@package biosimspace
@author  Lester Hedges
@brief   A class for storing minimisation protocols.
"""

from .protocol import Protocol

from warnings import warn

class Minimisation(Protocol):
    """A class for storing minimisation protocols."""

    def __init__(self, steps=10000):
        """Constructor.

           Keyword arguments:

           steps -- The maximum number of steps to perform.
        """

        # Set the number of steps.
        self.setSteps(steps)

    def getSteps(self):
        """Return the maximum number of steps."""
        return self._steps

    def setSteps(self, steps):
        """Set the maximum number of steps."""

        if type(steps) is not int:
            raise TypeError("'steps' must be of type 'int'")

        if steps <= 0:
            warn("Number of steps must be greater than zero. Using default (10000).")
            self._steps = 10000

        else:
            self._steps = steps
