"""
@package biosimspace
@author  Lester Hedges
@brief   A class for storing minimisation protocols.
"""

from .protocol_type import ProtocolType

from warnings import warn

class Minimisation():
    """A class for storing minimisation protocols."""

    def __init__(self, method="conjugate-gradient", steps=1000, temperature=300):
        """Constructor.

           Keyword arguments:

           method      -- The minimisation method.
           step        -- The maximum number of steps to perform.
           temperature -- The system temperature (in Kelvin).
        """

        # Set the protocol type.
        self._type = ProtocolType.MINIMISATION

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

    def type(self):
        """Return the protocol type."""
        return self._type
