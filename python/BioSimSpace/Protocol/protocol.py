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

class Protocol():
    """A base class for holding simulation protocols."""

    def __init__(self, protocol_type, gas_phase=False):
        """Constructor.

           Positional arguments:

           protocol_type -- The type of protocol.

           Keyword arguments:

           gas_phase     -- Whether this is a gas phase simulation.
        """

	# Don't allow user to create an instance of this base class.
        if type(self) == Protocol:
            raise Exception("<Protocol> must be subclassed.")

        # Set the protocol type.
        self._type = protocol_type

        # Set the gas phase flag.
        self.gas_phase = gas_phase

    def type(self):
        """Return the protocol type."""
        return self._type

    @property
    def gas_phase(self):
        """Return whether this is a gas phase simulation."""
        return self._gas_phase

    @gas_phase.setter
    def gas_phase(self, gas_phase):
        """Set the gas phase flag."""

        if type(gas_phase) is bool:
            self._gas_phase = gas_phase

        else:
            warn("Non-boolean gas phase flag. Defaulting to False!")
            self._gas_phase = False
