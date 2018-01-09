"""
@package biosimspace
@author  Lester Hedges
@brief   A class for storing production protocols.
"""

from .protocol import Protocol, ProtocolType

class Production(Protocol):
    """A class for storing production protocols."""

    def __init__(self):
        """ Constructor. """

        # Call the base class constructor.
        super().__init__(ProtocolType.PRODUCTION)

    def type(self):
        """Return the protocol type."""
        return self._type
