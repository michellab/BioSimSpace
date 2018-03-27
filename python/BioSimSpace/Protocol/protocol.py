"""
@package biosimspace
@author  Lester Hedges
@brief   A base class for holding simulation protocols.
"""

class Protocol():
    """A base class for holding simulation protocols."""

    def __init__(self):
        """Constructor."""

	# Don't allow user to create an instance of this base class.
        if type(self) is Protocol:
            raise Exception("<Protocol> must be subclassed.")
