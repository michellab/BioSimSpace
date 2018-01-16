"""
@package biosimspace
@author  Lester Hedges
@brief   A class for storing minimisation protocols.
"""

from .protocol import Protocol, ProtocolType

from warnings import warn

class Minimisation(Protocol):
    """A class for storing minimisation protocols."""

    def __init__(self, method="conjugate-gradient", steps=1000, temperature=300, gas_phase=False):
        """Constructor.

           Keyword arguments:

           method      -- The minimisation method.
           step        -- The maximum number of steps to perform.
           temperature -- The system temperature (in Kelvin).
           gas_phase   -- Whether this a gas phase simulation.
        """

        # Call the base class constructor.
        super().__init__(ProtocolType.MINIMISATION, gas_phase)

        # Set the minimisation method.
        self._method = method

        # Set the number of steps.
        self._steps = steps

        # Set the system temperature.
        self._temperature = temperature

    def type(self):
        """Return the protocol type."""
        return self._type

    @property
    def method(self):
        """Return the minimisation method."""
        return self._method

    @method.setter
    def method(self, method):
        """Set the 'method' member variable."""
        # Check that the minimisation method is valid.
        if method not in ["conjugate-gradient", "steepest-descent"]:
            msg = ('Invalid minimisation method: "{x}". '
                   'Allowed values are "conjugate-gradient" or "steepest-descent". '
                   'Using default: "conjugate-gradient".'
                   ).format(x=method)

            warn(msg)
            self._method = "conjugate-gradient"

        else:
            self._method = method

    @property
    def steps(self):
        """Return the maximum number of steps."""
        return self._steps

    @steps.setter
    def steps(self, steps):
        """Set the maximum number of steps."""

        if steps <= 0:
            warn("Number of steps must be greater than zero. Using default (5000).")
            self._steps = 5000

        else:
            self._steps = steps

    @property
    def temperature(self):
        """Return the system temperature."""
        return self._temperature

    @temperature.setter
    def temperature(self, temperature):
        """Set the system temperature."""

        if temperature <= 0:
            warn("Temperature must be positive. Using default (300 K).")
            self._temperature = 300

        else:
            self._temperature = temperature

    def conjugateGradient(self):
        """Whether the minimisation method is conjugate gradient."""
        return self._method == "conjugate-gradient"
