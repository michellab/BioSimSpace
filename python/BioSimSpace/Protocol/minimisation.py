"""
@package biosimspace
@author  Lester Hedges
@brief   A class for storing minimisation protocols.
"""

from .protocol_type import ProtocolType

class Minimisation():
    """ A class for storing minimisation protocols. """

    def __init__(self):
        """ Constructor. """

        # Set the protocol type.
        self._type = ProtocolType.MINIMISATION

    def type(self):
        """ Return the protocol type. """
        return self._type
