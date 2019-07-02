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
Functionality for torsion based collective variables.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["Torsion"]

from ._collective_variable import CollectiveVariable as _CollectiveVariable
from ...Types import Angle as _Angle

class Torsion(_CollectiveVariable):
    """A class for torsion based collective variables."""

    def __init__(self, atoms, lower_bound=None, upper_bound=None, grid=None, pbc=True):
        """Constructor.

           Parameters
           ----------

           atoms :  [int, int, int, int]
               The indices of the four atoms involved in the torsion.

           lower_bound : :class:`Bound <BioSimSpace.Metadynamics.Bound>`
               A lower bound on the value of the collective variable.

           upper_bound : :class:`Bound <BioSimSpace.Metadynamics.Bound>`
               An upper bound on the value of the collective variable.

           grid : :class:`Grid <BioSimSpace.Metadynamics.Grid>`
               The grid on which the collective variable will be sampled.
               This can help speed up long metadynamics simulations where
               the number of Guassian kernels can become prohibitive.

           pbc : bool
               Whether to use periodic boundary conditions when computing the
               collective variable.
        """

        # Call the base class constructor.
        super().__init__(pbc)

        # Initialise member data.
        self._atoms = None
        self._lower_bound = None
        self._upper_bound = None
        self._grid = None

        # Set the required parameters.

        self.setAtoms(atoms)

        # Set the optional parameters.

        if lower_bound is not None:
            self.setLowerBound(lower_bound)
        if upper_bound is not None:
            self.setUpperBound(upper_bound)
        if grid is not None:
            self.setGrid(grid)

        # Validate that the state is self-consistent.
        self._validate()

        # Flag that the object has been instantiated, i.e. it is no longer "new".
        self._is_new_object = False

    def __str__(self):
        """Return a human readable string representation of the object."""
        string = "<BioSimSpace.Metadynamics.CollectiveVariable.Distance: "
        string += "atoms=%s" % self._atoms
        if self._lower_bound is not None:
            string += ", lower_bound=%s" % self._lower_bound
        if self._upper_bound is not None:
            string += ", upper_bound=%s" % self._upper_bound
        if self._grid is not None:
            string += ", grid=%s" % self._grid
        string += ", pbc=%s"% self._pbc
        string += ">"
        return string

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return self.__str__()

    def setAtoms(self, atoms):
        """Set the atoms for which the torsion will be calculated.
           measured.

           Parameters
           ----------

           atoms : [int, int, int, int]
               The atoms involved in the torsion.
               will be measured from.
        """

        # Convert tuples to a list.
        if type(atoms) is tuple:
            atoms = list(atoms)

        # List of atom indices.
        if type(atoms) is list and all(isinstance(x, int) for x in atoms):
            pass
        else:
            raise TypeError("'atoms' must be of list of 'int' types.")

        if len(atoms) != 4:
            raise ValueError("'atoms' must contain four indices.")

        # All okay, set the value.
        self._atoms = atoms

    def getAtoms(self):
        """Return list of atom indices involved in the torsion.

           Returns
           -------

           atoms : [int, int, int, int]
               The atom indices involved in the torsion.
        """
        return self._atoms

    def _validate(self):
        """Internal function to check that the object is in a consistent state."""

        if self._lower_bound is not None:
            if type(self._lower_bound.getValue()) is not _Angle:
                raise TypeError("'lower_bound' must be of type 'BioSimSpace.Types.Angle'")
        if self._upper_bound is not None:
            if type(self._upper_bound.getValue()) is not _Angle:
                raise TypeError("'upper_bound' must be of type 'BioSimSpace.Types.Angle'")
        if self._lower_bound is not None and self._upper_bound is not None:
            if self._lower_bound.getValue() >= self._upper_bound.getValue():
                raise TypeError("'lower_bound' must less than 'upper_bound'")

        if self._grid is not None:
            if type(self._grid.getMinimum()) is not _Angle:
                raise TypeError("'grid' minimum must be of type 'BioSimSpace.Types.Angle'")
            if type(self._grid.getMaximum()) is not _Angle:
                raise TypeError("Grid 'maximum' must be of type 'BioSimSpace.Types.Angle'")
            if self._lower_bound is not None and self._grid.getMinimum() > self._lower_bound.getValue():
                raise ValueError("'lower_bound' is less than 'grid' minimum.")
            if self._upper_bound is not None and self._grid.getMaximum() < self._upper_bound.getValue():
                raise ValueError("'upper_bound' is greater than 'grid' maximum.")
