"""
@package biosimspace
@author  Lester Hedges
@brief   A base class for holding simulation protocols.
"""

from enum import Enum

class ProtocolType(Enum):
    """An enum class containing the list of supported simulation protocols."""
    MINIMISATION  = 1
    EQUILIBRATION = 2
    PRODUCTION    = 3
    CUSTOM        = 4

class Protocol():
    """A base class for holding simulation protocols."""

    def __init__(self, protocol_type, gas_phase=False):
        """Constructor.

           Positional arguments:

           protocol_type -- The type of protocol.
        """

	# Don't allow user to create an instance of this base class.
        if type(self) is Protocol:
            raise Exception("<Protocol> must be subclassed.")

        # Set the protocol type.
        self._type = protocol_type

    def getType(self):
        """Return the protocol type."""
        return self._type
