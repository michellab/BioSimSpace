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
Functionality for handling collective variables for metadynamics simulations.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["CollectiveVariable"]

class CollectiveVariable():
    """A base class for holding collective variables."""

    def __init__(self, pbc=True):
        """Constructor.

           Parameters
           ----------

           pbc : bool
              Whether to use periodic boundaries conditions.
        """

	# Don't allow user to create an instance of this base class.
        if type(self) is CollectiveVariable:
            raise Exception("<CollectiveVariable> must be subclassed.")

        # Whether this is a newly created object.
        self._is_new_object = True

        self.setPeriodicBoundaries(pbc)

    def setPeriodicBoundaries(self, pbc):
        """Set whether to use periodic_boundaries when calculating the
           collective variable.

           Parameters
           ----------

           pbc : bool
               Whether to use periodic boundaries conditions.
        """
        if type(pbc) is not bool:
            raise TypeError("'pbc' must be of type 'bool'")
        self._pbc = pbc

    def getPeriodicBoundaries(self):
        """Return whether to take account of periodic boundary conditions
           when computing the collective variable.

           Returns
           -------

           pbc : bool
               Whether to use periodic boundaries conditions.
        """
        return self._pbc
