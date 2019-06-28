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

from .._bound import Bound as _Bound
from .._grid import Grid as _Grid

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

    def setLowerBound(self, lower_bound=None):
        """Set a lower bound on the value of the collective variable.

           Parameters
           ----------

           lower_bound : :class:`Bound <BioSimSpace.Metadynamics.Bound>`
               A lower bound on the value of the collective variable.
        """

        if type(lower_bound) is None:
            self._lower_bound = None
            return

        if type(lower_bound) is not _Bound:
            raise TypeError("'lower_bound' must be of type 'BioSimSpace.Metadynamics.Bound'")

        # Store the existing value.
        old_value = self._lower_bound

        # Set the new value.
        self._lower_bound = lower_bound

        # If we are modifying an existing object, then check for consistency.
        if not self._is_new_object:
            try:
                self._validate()
            except:
                self._lower_bound = old_value
                raise

    def getLowerBound(self):
        """Get the lower bound on the collective variable.

           Returns
           -------

           lower_bound : :class:`Bound <BioSimSpace.Metadynamics.Bound>`
               The lower bound on the value of the collective variable.
        """
        return self._lower_bound

    def setUpperBound(self, upper_bound=None):
        """Set an upper bound on the value of the collective variable.

           Parameters
           ----------

           upper_bound : :class:`Bound <BioSimSpace.Metadynamics.Bound>`
               An upper bound on the value of the collective variable.
        """

        if type(upper_bound) is None:
            self._upper_bound = None
            return

        if type(upper_bound) is not _Bound:
            raise TypeError("'upper_bound' must be of type 'BioSimSpace.Metadynamics.Bound'")

        # Store the existing value.
        old_value = self._upper_bound

        # Set the new value.
        self._upper_bound = upper_bound

        # If we are modifying an existing object, then check for consistency.
        if not self._is_new_object:
            try:
                self._validate()
            except:
                self._upper_bound = old_value
                raise

    def getUpperBound(self):
        """Get the upper bound on the collective variable.

           Returns
           -------

           upper_bound : :class:`Bound <BioSimSpace.Metadynamics.Bound>`
               The upper bound on the value of the collective variable.
        """
        return self._upper_bound

    def setGrid(self, grid=None):
        """Set a grid on which the collective variable will be sampled.
           Call with no arguments to clear the grid.

           Parameters
           ----------

           grid : :class:`Grid <BioSimSpace.Metadynamics.Grid>`
               A grid for the collective variable.
        """

        if grid is None:
            self._grid = None
            return

        if type(grid) is not _Grid:
            raise TypeError("'grid' must be of type 'BioSimSpace.Metadynamics.Grid'")

        # Store the existing value.
        old_value = self._grid

        # Set the new value.
        self._grid = grid

        # If we are modifying an existing object, then check for consistency.
        if not self._is_new_object:
            try:
                self._validate()
            except:
                self._grid = old_value
                raise

    def getGrid(self):
        """Get the grid on which the collective variable is sampled.

           Returns
           -------

           grid : :class:`Grid <BioSimSpace.Metadynamics.Grid>`
               The grid on which the collective variable is sampled.
        """
        return self._grid

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

    def _validate(self):
        """Internal function to validate that the object is in a consistent
           state.
        """
        raise NotImplementedError("'_validate' must be overloaded in derived classes!")
