######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2020
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
Functionality for distance based collective variables.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["Funnel", "funnelMaker"]

from math import ceil as _ceil

from ._collective_variable import CollectiveVariable as _CollectiveVariable
from .._bound import Bound as _Bound
from .._grid import Grid as _Grid
from ..._SireWrappers import Molecule as _Molecule
from ..._SireWrappers import System as _System
from ...Types import Coordinate as _Coordinate
from ...Types import Length as _Length

class Funnel(_CollectiveVariable):
    """A class for a funnel collective variable."""

    def __init__(self, atoms0, atoms1, ligand=None, hill_width=_Length(0.025, "nanometer"),
            lower_bound=None, upper_bound=None, grid=None, pbc=False):
        """Constructor.

           Parameters
           ----------

           atoms0 : [int, int, ...]
               A list of atom indices that define the origin of the funnel.

           atoms1 : [int, int, ...]
               A list of atom indices that define the inflection point of the funnel.

           hill_width : :class:`Length <BioSimSpace.Types.Length>`
               The width of the Gaussian hill used to sample this variable.

           lower_bound : :class:`Bound <BioSimSpace.Metadynamics.Bound>`
               A lower bound on the value of the collective variable. This is
               used to constrain the distance along projection axis of the
               funnel.

           upper_bound : :class:`Bound <BioSimSpace.Metadynamics.Bound>`
               An upper bound on the value of the collective variable. This is
               used to constrain the distance along projection axis of the
               funnel.

           grid : :class:`Grid <BioSimSpace.Metadynamics.Grid>`
               The grid on which the collective variable will be sampled.
               This can help speed up long metadynamics simulations where
               the number of Gaussian kernels can become prohibitive.

           pbc : bool
               Whether to use periodic boundary conditions when computing the
               collective variable. For funnel metadynamics it is advisable
               that this is set to False.
        """

        # Call the base class constructor.
        super().__init__(pbc)

        # This collective variable has two components: the projection along
        # the funnel axis and the extension orthogonal to it.
        self._num_components = 2

        # Initialise member data.
        self._atoms0 = None
        self._atoms1 = None
        self._lower_bound = None
        self._upper_bound = None
        self._grid = None

        # Set the required parameters.
        self.setAtoms0(atoms0)
        self.setAtoms1(atoms1)
        self.setHillWidth(hill_width)

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
        string = "<BioSimSpace.Metadynamics.CollectiveVariable.Funnel: "
        string += "atoms0=%s" % self._atoms0
        string += ", atoms1=%s" % self._atoms1
        string += ", hill_width=%s" % self._hill_width
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

    def setAtoms0(self, atoms0):
        """Set the list of atom indices whose center-of-mass define the origin
           of the funnel.

           Parameters
           ----------

           atoms0 : [int, int, ...]
               The list of atom indices whose center-of-mass define the origin
               of the funnel.
        """

        # Convert tuples to a list.
        if type(atoms0) is tuple:
            atoms0 = list(atoms0)

        # List of atom indices.
        if type(atoms0) is list and all(isinstance(x, int) for x in atoms0):
            pass

        # Invalid type.
        else:
            raise TypeError("'atoms0' must be of list of 'int' types.")

        # All okay, set the value.
        self._atoms0 = atoms0

    def getAtoms0(self):
        """Return indices of the atoms whose center-of-mass defines the origin
           of the funnel.

           Returns
           -------

           atoms0 : [int]
               The list of atom indices whose center-of-mass define the origin
               of the funnel.
        """
        return self._atoms0

    def setAtoms1(self, atoms1):
        """Set the list of atom indices whose center-of-mass define the
           inflection point of the funnel.

           Parameters
           ----------

           atoms1 : [int, int, ...]
               The list of atom indices whose center-of-mass define the
               inflection point of the funnel.
        """

        # Convert tuples to a list.
        if type(atoms1) is tuple:
            atoms1 = list(atoms1)

        # List of atom indices.
        if type(atoms1) is list and all(isinstance(x, int) for x in atoms1):
            pass

        # Invalid type.
        else:
            raise TypeError("'atoms1' must be of list of 'int' types.")

        # All okay, set the value.
        self._atoms1 = atoms1

    def getAtoms1(self):
        """Return indices of the atoms whose center-of-mass defines the
           inflection point of the funnel.

           Returns
           -------

           atoms1 : [int]
               The list of atom indices whose center-of-mass define the
               inflection point of the funnel.
        """
        return self._atoms1

    def setHillWidth(self, hill_width):
        """Set the width of the Gaussian hills used to bias this collective
           variable.

           Parameters
           ----------

           hill_width : :class:`Length <BioSimSpace.Types.Length>`
               The width of the Gaussian hill.
        """
        if type(hill_width) is not _Length:
            raise TypeError("'hill_width' must be of type 'BioSimSpace.Types.Length'")

        if hill_width.magnitude() < 0:
            raise ValueError("'hill_width' must have a magnitude of > 0")

        # Convert to the internal unit.
        self._hill_width = hill_width.nanometers()

    def getHillWidth(self):
        """Return the width of the Gaussian hill used to bias this collective
           variable.

           Returns
           -------

           hill_width : :class:`Length <BioSimSpace.Types.Length>`
               The width of the Gaussian hill.
        """
        return self._hill_width

    def _validate(self):
        """Internal function to check that the object is in a consistent state."""

        if self._lower_bound is not None:
            if type(self._lower_bound.getValue()) is not _Length:
                raise TypeError("'lower_bound' must be of type 'BioSimSpace.Types.Length'")
            # Convert to default unit.
            self._lower_bound.setValue(self._lower_bound.getValue().nanometers())
        if self._upper_bound is not None:
            if type(self._upper_bound.getValue()) is not _Length:
                raise TypeError("'upper_bound' must be of type 'BioSimSpace.Types.Length'")
            # Convert to default unit.
            self._upper_bound.setValue(self._upper_bound.getValue().nanometers())
        if self._lower_bound is not None and self._upper_bound is not None:
            if self._lower_bound.getValue() >= self._upper_bound.getValue():
                raise TypeError("'lower_bound' must less than 'upper_bound'")

        if self._grid is not None:
            if type(self._grid.getMinimum()) is not _Length:
                raise TypeError("'grid' minimum must be of type 'BioSimSpace.Types.Length'")
            # Convert to default unit.
            self._grid.setMinimum(self._grid.getMinimum().nanometers())
            if type(self._grid.getMaximum()) is not _Length:
                raise TypeError("Grid 'maximum' must be of type 'BioSimSpace.Types.Length'")
            # Convert to default unit.
            self._grid.setMaximum(self._grid.getMaximum().nanometers())
            if self._lower_bound is not None and self._grid.getMinimum() > self._lower_bound.getValue():
                raise ValueError("'lower_bound' is less than 'grid' minimum.")
            if self._upper_bound is not None and self._grid.getMaximum() < self._upper_bound.getValue():
                raise ValueError("'upper_bound' is greater than 'grid' maximum.")

            # If the number of bins isn't specified, estimate it out from the hill width.
            if self._grid.getBins() is None:
                grid_range = (self._grid.getMaximum() - self._grid.getMinimum()).magnitude()
                num_bins = _ceil(5.0 * (grid_range / self._hill_width.magnitude()))
                self._grid.setBins(num_bins)

def funnelMaker(system, protein=None, ligand=None,
        search_radius=_Length(10, "A"), alpha_carbon_name="CA"):
    """Constructor.

       Parameters
       ----------

       system : :class:`System <BioSimSpace._SireWrappers.System>`
           The system of interest. This is assumed to be a solvated
           protein-ligand system.

       protein : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, int
           The protein molecule. This can either be a Molecule object, or
           the index of the protein within the passed system. If None is
           passed then we assume that the protein is the largest molecule
           in the system.

       ligand : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, int, \
                :class:`Coordinate <BioSimSpace.Types.Coordinate>`
           The ligand molecule. This can either be a Molecule object, the
           index of the ligand within the passed system, or the position of
           the binding site within the protein. If None is passed then we
           assume that the ligand is the second largest molecule in the system.

       search_radius : :class:`Length <BioSimSpace.Types.Length>`
           The search radius around the binding site used to determine the
           funnel origin and inflection point.

       alpha_carbon_name : str
           The name of alhpa carbon atoms in the system topology.

       Returns
       -------

       atoms0 : [int]
           A list of atom indices that define the origin of the funnel.

       atoms1 : [int]
           A list of atom indices that define the inflection point of the funnel.
    """

    # Validate the input.

    # System.
    if type(system) is not _System:
        raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'.")

    # Protein.
    if type(protein) is _Molecule:
        # Make sure the molecule exists in the system.
        try:
            system.getIndex(protein)
        except KeyError:
            raise(ValueError("The passed 'protein' was not found in the 'system'."))

        protein = protein.copy()

    elif type(protein) is int:
        try:
            protein = system[protein]
        except:
            ValueError(f"Couldn't find 'protein' index {protein} in the 'system'.")

    elif type(protein) is _Coordinate:
        pass

    elif protein is None:
        # Store the number of atoms in each molecule.
        atom_nums = []
        for mol in system:
            atom_nums.append(mol.nAtoms())

        # Sort the molecule indices by the number of atoms they contain.
        sorted_nums = sorted((value, idx) for idx, value in enumerate(atom_nums))

        # Set the protein to the largest molecule.
        protein = system[sorted_nums[-1][1]]

    else:
        raise TypeError("'protein' must be of type 'BioSimSpace._SireWrappers.Molecule' or 'None'.")

    # Ligand.
    if type(ligand) is _Molecule:
        # Make sure the molecule exists in the system.
        try:
            system.getIndex(ligand)
        except KeyError:
            raise(ValueError("The passed 'ligand' was not found in the 'system'."))

    elif type(ligand) is int:
        try:
            ligand = system[ligand]
        except:
            ValueError(f"Couldn't find 'ligand' index {ligand} in the 'system'.")

    elif type(ligand) is _Coordinate:
        pass

    elif ligand is None:
        # Store the number of atoms in each molecule.
        atom_nums = []
        for mol in system:
            atom_nums.append(mol.nAtoms())

        # Sort the molecule indices by the number of atoms they contain.
        sorted_nums = sorted((value, idx) for idx, value in enumerate(atom_nums))

        # Set the ligand to the second largest molecule.
        ligand = system[sorted_nums[-2][1]]

    else:
        raise TypeError("'ligand' must be of type 'BioSimSpace._SireWrappers.Molecule' or 'None'.")

    # Search radius.
    if type(search_radius) is not _Length:
        raise TypeError("'search_radius' must be of type 'BioSimSpace.Types.Length'.")

    # Alpha carbon name.
    if type(alpha_carbon_name) is not str:
        raise TypeError("'alpha_carbon_name' must be of type 'str'.")
