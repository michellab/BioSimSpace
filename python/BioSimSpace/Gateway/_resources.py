######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2019
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# BioSimSpace is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BioSimSpace. If not, see <http://www.gnu.org/licenses/>.
#####################################################################

"""
Functionality for finding and managing hardware resources.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

import argparse as _argparse

__all__ = ["ResourceManager"]

class ResourceManager():
    """A class for finding and managing hardware resources."""

    def __init__(self):
        """Constructor."""

        # Set default values.
        self._nodes = None
        self._cores = None
        self._gpus = None

        # Create the argument parser.
        self._parser = _argparse.ArgumentParser(description="Command-line parser for hardware resources",
                                                add_help=False, allow_abbrev=False)

        # Add the arguments.
        self._parser.add_argument("--nodes", type=int, help="The number of harwdare nodes.")
        self._parser.add_argument("--cores", type=int, help="The number of harwdare cores.")
        self._parser.add_argument("--gpus",  type=int, help="The number of harwdare graphics processors.")

    def _initialise(self):
        """Initialise the resource manager."""

        # Parse the arguments into a dictionary.
        args = vars(self._parser.parse_known_args()[0])

        # Loop over all of the arguments and set the values.
        for key, value in args.items():
            if key == "nodes":
                if value is not None:
                    self._nodes = int(value)
            elif key == "cores":
                if value is not None:
                    self._cores = int(value)
            elif key == "gpus":
                if value is not None:
                    self._gpus = int(value)

    def getNodes(self):
        """Return the number of nodes.

           Returns
           -------

           nodes : int
               The number of nodes.
        """
        return self._nodes

    def setNodes(self, nodes):
        """Set the number of nodes.

           Parameters
           ----------

           nodes : int
               The number of nodes.
        """

        if type(nodes) is not int:
            raise TypeError("'nodes' must be of type 'int'.")

        if nodes < 0:
            raise ValueError("'nodes' cannot be negative!")

        self._nodes = nodes

    def getCores(self):
        """Return the number of cores.

           Returns
           -------

           cores : int
               The number of cores.
        """
        return self._cores

    def setCores(self, cores):
        """Set the number of cores.

           Parameters
           ----------

           cores : int
               The number of cores.
        """

        if type(cores) is not int:
            raise TypeError("'cores' must be of type 'int'.")

        if cores < 0:
            raise ValueError("'cores' cannot be negative!")

        self._cores = cores

    def getGPUs(self):
        """Return the number of GPUs.

           Returns
           -------

           gpus : int
               The number of GPUs.
        """
        return self._gpus

    def setGPUs(self, gpus):
        """Set the number of GPUs.

           Parameters
           ----------

           gpus : int
               The number of GPUs.
        """

        if type(gpus) is not int:
            raise TypeError("'gpus' must be of type 'int'.")

        if gpus < 0:
            raise ValueError("'gpus' cannot be negative!")

        self._gpus = gpus
