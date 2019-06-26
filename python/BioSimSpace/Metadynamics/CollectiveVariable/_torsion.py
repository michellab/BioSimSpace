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

    def __init__(self, atoms, lower_bound=None, upper_bound=None, pbc=True):
        """Constructor.

           Parameters
           ----------

           atoms :  [int, int, int, int]
               The indices of the four atoms involved in the torsion.

           lower_bound : :class:`Length <BioSimSpace.Types.Angle>`
               The minimum value of the collective variable.

           upper_bound : :class:`Length <BioSimSpace.Types.Angle>`
               The maximum value of the collective variable.

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

        # Set the required parameters.

        self.setAtoms(atoms)

        # Set the optional parameters.

        if lower_bound is not None:
            self.setLowerBound(lower_bound)
        if upper_bound is not None:
            self.setUpperBound(upper_bound)

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

    def setLowerBound(self, lower_bound=None):
        """Set the minimum value of the collective variable. Can be called with
           no arguments to clear the data.

           Parameters
           ----------

           lower_bound : :class: `Angle <BioSimSpace.Types.Angle>`
               The minimum value of the collective variable.
        """
        if type(lower_bound) is None:
            self._lower_bound = None
            return

        if type(lower_bound) is not _Angle:
            raise TypeError("'lower_bound' must be of type 'BioSimSpace.Types.Angle'")

        # Store the existing value.
        old_value = self._lower_bound

        self._lower_bound = lower_bound

        # If we are modifying an existing object, then check for consistency.
        if not self._is_new_object:
            try:
                self._validate()
            except:
                self._lower_bound = old_value
                raise

    def getLowerBound(self):
        """Get the minimum value of the collective variable.

           Returns
           -------

           lower_bound : :class: `Length <BioSimSpace.Types.Angle>`
               The minimum value of the collective variable.
        """
        return self._lower_bound

    def setUpperBound(self, upper_bound=None):
        """Set the maximum value of the collective variable. Can be called with
           no arguments to clear the data.

           Parameters
           ----------

           upper_bound : :class: `Length <BioSimSpace.Types.Angle>`
               The maximum value of the collective variable.
        """
        if type(upper_bound) is None:
            self._upper_bound = None
            return

        if type(upper_bound) is not _Angle:
            raise TypeError("'upper_bound' must be of type 'BioSimSpace.Types.Angle'")

        # Store the existing value.
        old_value = self._lower_bound

        self._upper_bound = upper_bound

        # If we are modifying an existing object, then check for consistency.
        if not self._is_new_object:
            try:
                self._validate()
            except:
                self._upper_bound = old_value
                raise

    def getUpperBound(self):
        """Get the maximum value of the collective variable.

           Returns
           -------

           upper_bound : :class: `Length <BioSimSpace.Types.Angle>`
               The maximum value of the collective variable.
        """
        return self._upper_bound

    def _validate(self):
        """Internal function to check that the object is in a consistent state."""

        if self._lower_bound is not None and self._upper_bound is not None:
            if self._lower_bound >= self._upper_bound:
                raise ValueError("'lower_bound' must be less than 'upper_bound'")
