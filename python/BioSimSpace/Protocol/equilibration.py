"""
@package biosimspace
@author  Lester Hedges
@brief   A class for storing equilibration protocols.
"""

from .protocol_type import ProtocolType

class Equilibration():
    """ A class for storing equilibration protocols. """

    def __init__(self):
        """ Constructor. """

        # Set the protocol type.
        self._type = ProtocolType.EQUILIBRATION

    def type(self):
        """ Return the protocol type. """
        return self._type
