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
A thin wrapper around Sire.Mol.Molecules. This is an internal package and should
not be directly exposed to the user.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["Molecules"]

from Sire import Maths as _SireMaths
from Sire import Mol as _SireMol
from Sire import System as _SireSystem

from BioSimSpace import _isVerbose
from BioSimSpace.Types import Length as _Length
from BioSimSpace import Units as _Units

from ._sire_wrapper import SireWrapper as _SireWrapper

class Molecules(_SireWrapper):
    """An immutable container class for storing molecules."""

    def __init__(self, molecules):
        """Constructor.

           Parameters
           ----------

           molecules : Sire.Mol.Molecules, Sire.System.System, \
                       :class:`System <BioSimSpace._SireWrappers.System>` \
                       [:class:`Molecule <BioSimSpace._SireWrappers.Molecule>`]
               A Sire Molecues object, a Sire or BioSimSpace System object,
               or a list of BioSimSpace Molecule objects.
        """

        # Check that the molecules argument is valid.

        # Convert tuple to list.
        if type(molecules) is tuple:
            molecules = list(molecules)

        # A Sire Molecules object.
        if type(molecules) is _SireMol.Molecules:
            super().__init__(molecules)

        # A Sire System object.
        elif type(molecules) is _SireSystem.System:
            super().__init__(molecules.molecules())

        # A BioSimSpace System object.
        elif type(molecules) is _System:
            super().__init__(molecules._sire_object.molecules())

        # Another BioSimSpace Molecule object.
        elif type(molecules) is list and all(isinstance(x, _Molecule) for x in molecules):
            mols = _SireMol.Molecules()
            for molecule in molecules:
                mols.add(molecule._sire_object)
                super().__init__(mols)

        # Invalid type.
        else:
            raise TypeError("'molecules' must be of type 'Sire.Mol.Molecules' "
                            "'Sire.System.System', 'BioSimSpace._SireWrappers.System', "
                            "or a list of 'BioSimSpace._SireWrappers.Molecule' types.")

        # Store the list of MolNums.
        self._mol_nums = self._sire_object.molNums()

        # Sort the molecule numbers.
        self._mol_nums.sort()

        # Store the number of molecules.
        self._num_mols = self._sire_object.nMolecules()

        # Initialise the iterator counter.
        self._iter_count = 0

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.Molecules: nMolecules=%d>" % self.nMolecules()

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return "<BioSimSpace.Molecules: nMolecules=%d>" % self.nMolecules()

    def __add__(self, other):
        """Addition operator."""

        # Convert tuple to a list.
        if type(other) is tuple:
            other = list(other)

        # Extract the molecules.
        molecules = self._sire_object

        # Validate the input.

        # A System object.
        if type(other) is _System:
            molecules.add(other._sire_object.molecules())

        # A Molecule object.
        elif type(other) is _Molecule:
            molecules.add(other._sire_object)

        # A Molecules object.
        elif type(other) is Molecules:
            molecules.add(other._sire_object)

        # A list of Molecule objects.
        elif type(other) is list and all(isinstance(x, Molecule) for x in other):
            for molecule in other:
                molecules.add(molecule._sire_object)

        # Unsupported.
        else:
            raise TypeError("'other' must be of type 'BioSimSpace._SireWrappers.System', "
                            "'BioSimSpace._SireWrappers.Molecule', 'BioSimSpace._SireWrappers.Molecules' "
                            "or a list of 'BioSimSpace._SireWrappers.Molecule' types")

        # Create and return a new Molecules object.
        return Molecules(molecules)

    def __getitem__(self, key):
        """Get a molecule from the container."""

        # Slice.
        if type(key) is slice:

            # Create a list to hold the molecules.
            molecules = []

            # Iterate over the slice.
            for x in range(*key.indices(self._num_mols)):
                molecules.append(self[x])

            # Convert to a Molecules container and return.
            return Molecules(molecules)

        # Index.
        else:
            try:
                key = int(key)
            except:
                raise TypeError("'key' must be of type 'int'")

            if key < -self._num_mols or key > self._num_mols -1:
                raise IndexError("Molecules index is out of range.")

            if key < 0:
                key = key + self._num_mols

            # Extract and return the corresponding molecule.
            return _Molecule(self._sire_object.molecule(self._mol_nums[key]))

    def __setitem__(self, key, value):
        """Set a molecule in the container."""
        raise TypeError("'Molecules' object does not support assignment.")

    def __iter__(self):
        """An iterator for the object."""
        # Reset the iterator counter and return the object.
        self._iter_count = 0
        return self

    def __next__(self):
        """An iterator for the object."""

        # Stop if we've reached the end of the container.
        if self._iter_count == self._num_mols:
            raise StopIteration

        # Extract the next molecule in the container.
        molecule = self[self._iter_count]

        # Update the iterator counter.
        self._iter_count += 1

        # Return the molecule.
        return molecule

    def __len__(self):
        """Return the size of the container."""
        return self._num_mols

    def nMolecules(self):
        """Return the number of molecules in the system.

           Returns
           -------

           num_molecules : int
               The number of molecules in the system.
        """
        return self._num_mols

    def getMolecule(self, index):
        """Return the molecule at the given index.

           Parameters
           ----------

           index : int
               The index of the molecule.

           Returns
           -------

           molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
               The requested molecule.
        """
        return self[index]

    def toSystem(self):
        """Convert to a System object.

           Returns
           -------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
        """
        return _System(self)

    def charge(self, property_map={}, is_lambda1=False):
        """Return the total molecular charge.

           Parameters
           ----------

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }

           is_lambda1 : bool
              Whether to use the charge at lambda = 1 a molecule is merged.

           Returns
           -------

           charge : :class:`Charge <BioSimSpace.Types.Charge>`
               The molecular charge.
        """

        # Zero the charge.
        charge = 0 * _Units.Charge.electron_charge

        # Loop over all molecules and add the charge.
        for mol in self:
            charge += mol.charge(property_map, is_lambda1)

        # Return the total charge.
        return charge

    def translate(self, vector, property_map={}):
        """Translate each molecule in the container.

           Parameters
           ----------

           vector : [:class:`Length <BioSimSpace.Types.Length>`]
               The translation vector in Angstroms.

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Convert tuple to a list.
        if type(vector) is tuple:
            vector = list(vector)

        # Validate input.
        if type(vector) is list:
            vec = []
            for x in vector:
                if type(x) is int:
                    vec.append(float(x))
                elif type(x) is float:
                    vec.append(x)
                elif type(x) is _Length:
                    vec.append(x.angstroms().magnitude())
                else:
                    raise TypeError("'vector' must contain 'int', 'float', or "
                                    "'BioSimSpace.Types.Length' types only!")
        else:
            raise TypeError("'vector' must be of type 'list' or 'tuple'")

        if type(property_map) is not dict:
            raise TypeError("'property_map' must be of type 'dict'")

        # Translate each of the molecules in the container.
        for n in self._sire_object.molNums():
            # Copy the property map.
            _property_map = property_map.copy()

            # If this is a perturbable molecule, use the coordinates at lambda = 0.
            if self._sire_object.molecule(n).hasProperty("is_perturbable"):
                _property_map["coordinates"] = "coordinates0"

            mol = self._sire_object[n].move().translate(_SireMaths.Vector(vec), _property_map).commit()
            self._sire_object.update(mol)

    def search(self, query):
        """Search the molecules for atoms and residues. Search results will be
           reduced to their minimal representation, i.e. a residue containing
           a single atom will be returned as a atom.

           Parameters
           ----------

           query : str
               The search query.

           Returns
           -------

           results : [:class:`Atom <BioSimSpace._SireWrappers.Atom>`, \
                      :class:`Residue <BioSimSpace._SireWrappers.Residue>`, ...]
               A list of objects matching the search query.

           Examples
           --------

           Search for all residues named ALA.

           >>> result = molecules.search("resname ALA")

           Search for all oxygen or hydrogen atoms.

           >>> result = molecules.search("element oxygen or element hydrogen")

           Search for atom index 23.

           >>> result = molecule.search("atomidx 23")
        """

        if type(query) is not str:
            raise TypeError("'query' must be of type 'str'")

        # Intialise a list to hold the search results.
        results = []

        try:
            # Query the Sire system.
            search_result = self._sire_object.search(query)

        except Exception as e:
            msg = "'Invalid search query: %r" % query
            if _isVerbose():
                raise ValueError(msg) from e
            else:
                raise ValueError(msg) from None

        return _SearchResult(search_result)

    def _getAABox(self, property_map={}):
        """Get the axis-aligned bounding box for the molecular system.

           Parameters
           ----------

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }

           Returns
           -------

           aabox : Sire.Vol.AABox
               The axis-aligned bounding box for the molecule.
        """
        return self.toSystem()._getAABox()

# Import at bottom of module to avoid circular dependency.
from ._molecule import Molecule as _Molecule
from ._search_result import SearchResult as _SearchResult
from ._system import System as _System
