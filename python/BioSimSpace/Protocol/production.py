"""
@package biosimspace
@author  Lester Hedges
@brief   A class for storing production protocols.
"""

from .protocol_type import ProtocolType

class Production():
    """A class for storing production protocols."""

    def __init__(self):
        """ Constructor. """

        # Set the protocol type.
        self._type = ProtocolType.PRODUCTION

    def type(self):
        """Return the protocol type."""
        return self._type
