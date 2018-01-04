"""
@package biosimspace
@author  Lester Hedges
@brief   An enum class for holding simulation protocols.
"""

from enum import Enum

class ProtocolType(Enum):
    """ An enum class containing the list of supported simulation protocols. """
    MINIMISATION  = 1
    EQUILIBRATION = 2
    PRODUCTION    = 3
