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
Functionality for minimisation protocols.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["Minimisation"]

import warnings as _warnings

from ._protocol import Protocol as _Protocol

class Minimisation(_Protocol):
    """A class for storing minimisation protocols."""

    def __init__(self, steps=10000):
        """Constructor.

           Parameters
           ----------

           steps : int
               The maximum number of steps to perform.
        """

        # Call the base class constructor.
        super().__init__()

        # Set the number of steps.
        self.setSteps(steps)

    def __str__(self):
        """Return a human readable string representation of the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            return "<BioSimSpace.Protocol.Minimisation: steps=%d>" % self._steps

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            return "BioSimSpace.Protocol.Minimisation(steps=%d)" % self._steps

    def getSteps(self):
        """Return the maximum number of steps.

           Returns
           -------

           steps : int
               The maximum number of minimisation steps.
        """
        return self._steps

    def setSteps(self, steps):
        """Set the maximum number of steps.

           Parameters
           ----------

           steps : int
               The maximum number of minimisation steps.
        """
        if type(steps) is not int:
            raise TypeError("'steps' must be of type 'int'")

        if steps <= 0:
            _warnings.warn("Number of steps must be greater than zero. Using default (10000).")
            self._steps = 10000

        else:
            self._steps = steps
