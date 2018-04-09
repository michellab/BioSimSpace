"""
@package biosimspace
@author  Lester Hedges
@brief   A class for storing custom protocols.
"""

from .protocol import Protocol

from os import path

class Custom(Protocol):
    """A class for storing custom protocols."""

    def __init__(self, config):
        """Constructor.

           Keyword arguments:

           config -- The custom protocol configuration.
        """

        # Set the protocol configuration.
        self.setConfig(config)

    def getConfig(self):
        """Return the custom configuration."""
        return self._config.copy()

    def setConfig(self, config):
        """Set the custom configuration."""

        # Check that the passed configuration is a list of strings.
        if _is_list_of_strings(config):
            self._config = config

        # The user has passed a path to a file.
        elif path.isfile(config):

            # Clear the existing config.
            self._config = []

            # Read the contents of the file.
            with open(config, "r") as f:
                for line in f:
                    self._config.append(line.rstrip())
        else:
            raise ValueError("'config' must be a list of strings, or a file path.")

def _is_list_of_strings(lst):
    """Check whether the passed argument is a list of strings."""
    if lst and isinstance(lst, list):
        return all(isinstance(elem, str) for elem in lst)
    else:
        return False
