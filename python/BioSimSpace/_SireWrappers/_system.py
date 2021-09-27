######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2021
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
A thin wrapper around Sire.System.System. This is an internal package and should
not be directly exposed to the user.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["System"]

from Sire import IO as _SireIO
from Sire import Maths as _SireMaths
from Sire import Mol as _SireMol
from Sire import System as _SireSystem
from Sire import Units as _SireUnits
from Sire import Vol as _SireVol

from BioSimSpace import _isVerbose
from BioSimSpace._Exceptions import IncompatibleError as _IncompatibleError
from BioSimSpace.Types import Angle as _Angle
from BioSimSpace.Types import Length as _Length
from BioSimSpace import Units as _Units

from ._sire_wrapper import SireWrapper as _SireWrapper

class System(_SireWrapper):
    """A container class for storing molecular systems."""

    def __init__(self, system):
        """Constructor.

           Parameters
           ----------

           system : Sire.System.System, :class:`System <BioSimSpace._SireWrappers.System>`, \
                    Sire.Mol._Mol.Molecule, :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                    :class:`Molecules <BioSimSpace._SireWrappers.Molecules>` \
                    [:class:`Molecule <BioSimSpace._SireWrappers.Molecule>`]
               A Sire or BioSimSpace System object, a Sire or BioSimSpace Molecule object,
               a BioSimSpace Molecules object, or a list of BioSimSpace molecule objects.
        """

        # Check that the system is valid.

        # Convert tuple to a list.
        if type(system) is tuple:
            system = list(system)

        # A Sire System object.
        if type(system) is _SireSystem.System:
            super().__init__(system)

        # Another BioSimSpace System object.
        elif type(system) is System:
            super().__init__(system._sire_object)

        # A Sire Molecule object.
        elif type(system) is _SireMol.Molecule:
            sire_object = _SireSystem.System("BioSimSpace System.")
            super().__init__(sire_object)
            self.addMolecules(_Molecule(system))

        # A BioSimSpace Molecule object.
        elif type(system) is _Molecule:
            sire_object = _SireSystem.System("BioSimSpace System.")
            super().__init__(sire_object)
            self.addMolecules(system)

        # A BioSimSpace Molecules object.
        elif type(system) is _Molecules:
            sire_object = _SireSystem.System("BioSimSpace System.")
            super().__init__(sire_object)
            self.addMolecules(system)

        # A list of BioSimSpace Molecule objects.
        elif type(system) is list:
            if not all(isinstance(x, _Molecule) for x in system):
                raise TypeError("'system' must contain a list of 'BioSimSpace._SireWrappers.Molecule' types.")
            else:
                sire_object = _SireSystem.System("BioSimSpace System.")
                super().__init__(sire_object)
                self.addMolecules(system)

        # Invalid type.
        else:
            raise TypeError("'system' must be of type 'Sire.System.System', 'BioSimSpace._SireWrappers.System', "
                            " Sire.Mol.Molecule', 'BioSimSpace._SireWrappers.Molecule', "
                            "or a list of 'BioSimSpace._SireWrappers.Molecule' types.")

        # Flag that this object holds multiple atoms.
        self._is_multi_atom = True

        # Initialise dictionaries mapping MolNums to the number cumulative
        # number of atoms/residues in the system, up to that MolNum.
        self._atom_index_tally = {}
        self._residue_index_tally = {}

        # Initialise dictionary to map MolNum to MolIdx.
        self._molecule_index = {}

        # Store the molecule numbers.
        self._mol_nums = self._sire_object.molNums()

        # Initialise the iterator counter.
        self._iter_count = 0

        # Copy any fileformat property to each molecule.
        if "fileformat" in self._sire_object.propertyKeys():
            fileformat = self._sire_object.property("fileformat")
            for num in self._mol_nums:
                edit_mol = self._sire_object[num].edit()
                edit_mol = edit_mol.setProperty("fileformat", fileformat)
                self._sire_object.update(edit_mol.commit())
        else:
            # If a molecule has a fileformat property, use the first
            # that we find.
            for mol in self:
                if mol._sire_object.hasProperty("fileformat"):
                    self._sire_object.setProperty("fileformat", mol._sire_object.property("fileformat"))
                    break

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.System: nMolecules=%d>" % self.nMolecules()

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return "<BioSimSpace.System: nMolecules=%d>" % self.nMolecules()

    def __add__(self, other):
        """Addition operator."""

        # Create a copy of the current system.
        system = System(self._sire_object.__deepcopy__())

        # Add the new molecules.
        system.addMolecules(other)

        # Return the combined system.
        return system

    def __sub__(self, other):
        """Subtraction operator."""

        # Create a copy of the current system.
        system = System(self._sire_object.__deepcopy__())

        # Remove the molecules from the other system.
        if type(other) is System:
            system.removeMolecules(other.getMolecules())
        else:
            system.removeMolecules(other)

        # Return the new system.
        return system

    def __getitem__(self, key):
        """Get a molecule from the system."""

        # Slice.
        if type(key) is slice:

            # Create a list to hold the molecules.
            molecules = []

            # Iterate over the slice.
            for x in range(*key.indices(self.nMolecules())):
                molecules.append(self[x])

            # Convert to a Molecules container and return.
            return _Molecules(molecules)

        # Index.
        else:
            try:
                key = int(key)
            except:
                raise TypeError("'key' must be of type 'int'")

            if key < -self.nMolecules() or key > self.nMolecules() -1:
                raise IndexError("Molecules index is out of range.")

            if key < 0:
                key = key + self.nMolecules()

            # Extract and return the corresponding molecule.
            return _Molecule(self._sire_object[_SireMol.MolIdx(key)])

    def __setitem__(self, key, value):
        """Set a molecule in the container."""
        raise TypeError("'System' object does not support assignment.")

    def __iter__(self):
        """An iterator for the object."""
        # Reset the iterator counter and return the object.
        self._iter_count = 0
        return self

    def __next__(self):
        """An iterator for the object."""

        # Stop if we've reached the end of the container.
        if self._iter_count == self.nMolecules():
            raise StopIteration

        # Extract the next molecule in the container.
        molecule = self[self._iter_count]

        # Update the iterator counter.
        self._iter_count += 1

        # Return the molecule.
        return molecule

    def __len__(self):
        """Return the number of molecules in the system."""
        return self.nMolecules()

    def nMolecules(self):
        """Return the number of molecules in the system.

           Returns
           -------

           num_molecules : int
               The number of molecules in the system.
        """
        return self._sire_object.nMolecules()

    def nResidues(self):
        """Return the number of residues in the system.

           Returns
           -------

           num_residues : int
               The number of residues in the system.
        """

        tally = 0

        for n in self._sire_object.molNums():
            tally += self._sire_object[n].nResidues()

        return tally

    def nChains(self):
        """Return the number of chains in the system.

           Returns
           -------

           num_chains : int
               The number of chains in the system.
        """

        tally = 0

        for n in self._sire_object.molNums():
            tally += self._sire_object[n].nChains()

        return tally

    def nAtoms(self):
        """Return the number of atoms in the system.

           Returns
           -------

           num_atoms : int
               The number of atoms in the system.
        """

        tally = 0

        for n in self._sire_object.molNums():
            tally += self._sire_object[n].nAtoms()

        return tally

    def charge(self, property_map={}, is_lambda1=False):
        """Return the total molecular charge.

           Parameters
           ----------

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }

           is_lambda1 : bool
              Whether to use the charge at lambda = 1 if the molecule is merged.

           Returns
           -------

           charge : :class:`Charge <BioSimSpace.Types.Charge>`
               The molecular charge.
        """

        # Zero the charge.
        charge = 0 * _Units.Charge.electron_charge

        # Loop over all molecules and add the charge.
        for mol in self.getMolecules():
            charge += mol.charge(property_map, is_lambda1)

        # Return the total charge.
        return charge

    def fileFormat(self, property_map={}):
        """Return the file formats associated with the system.

           Parameters
           ----------

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }

           Returns
           -------

           format : str
              The file formats associated with the system.
        """

        if type(property_map) is not dict:
            raise TypeError("'property_map' must be of type 'dict'")

        prop = property_map.get("fileformat", "fileformat")

        try:
            return self._sire_object.property(prop).value()
        except:
            return None

    def addMolecules(self, molecules):
        """Add a molecule, or list of molecules to the system.

           Parameters
           ----------

           molecules : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                       :class:`Molecules <BioSimSpace._SireWrappers.Molecules>`, \
                       [:class:`Molecule <BioSimSpace._SireWrappers.Molecule>`], \
                       :class:`System <BioSimSpace._SireWrappers.System>`
              A Molecule, Molecules object, a list of Molecule objects, or a System containing molecules.
        """

        # Whether the molecules are in a Sire container.
        is_sire_container = False

        # Convert tuple to a list.
        if type(molecules) is tuple:
            molecules = list(molecules)

        # A Molecule object.
        if type(molecules) is _Molecule:
            molecules = [molecules]

        # A Molecules object.
        if type(molecules) is _Molecules:
            is_sire_container = True

        # A System object.
        elif type(molecules) is System:
            molecules = molecules.getMolecules()
            is_sire_container = True

        # A list of Molecule objects.
        elif type(molecules) is list and all(isinstance(x, _Molecule) for x in molecules):
            pass

        # Invalid argument.
        else:
            raise TypeError("'molecules' must be of type 'BioSimSpace._SireWrappers.Molecule', "
                            ", 'BioSimSpace._SireWrappers.System', or a list of "
                            "'BioSimSpace._SireWrappers.Molecule' types.")

        # The system is empty: create a new Sire system from the molecules.
        if self._sire_object.nMolecules() == 0:
            self._sire_object = self._createSireSystem(molecules)

        # Otherwise, add the molecules to the existing "all" group.
        else:
            if is_sire_container:
                if type(molecules) is _Molecules:
                    molecules = molecules._sire_object
                else:
                    try:
                        molecules = molecules._sire_object.at(_SireMol.MGNum(1))
                    except:
                        molecules = molecules._sire_object.molecules()
                self._sire_object.add(molecules, _SireMol.MGName("all"))
            else:
                for mol in molecules:
                    if mol._sire_object.number() in self._mol_nums:
                        raise ValueError("'BioSimSpace._SireWrappers.System' can only "
                                         "contain unique molecules. Use the 'copy' method "
                                         "of 'BioSimSpace._SireWrappers.Molecule' to "
                                         "create a new version of a molecule.")
                    self._sire_object.add(mol._sire_object, _SireMol.MGName("all"))

        # Reset the index mappings.
        self._reset_mappings()

        # Update the molecule numbers.
        self._mol_nums = self._sire_object.molNums()

    def removeMolecules(self, molecules):
        """Remove a molecule, or list of molecules from the system.

           Parameters
           ----------

           molecules : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                       :class:`Molecules <BioSimSpace._SireWrappers.Molecules>`, \
                       [:class:`Molecule <BioSimSpace._SireWrappers.Molecule>`]
              A Molecule, Molecules object, or list of Molecule objects.
        """

        # Whether the molecules are in a Sire container.
        is_sire_container = False

        # Convert tuple to a list.
        if type(molecules) is tuple:
            molecules = list(molecules)

        # A Molecule object.
        if type(molecules) is _Molecule:
            molecules = [molecules]

        # A Molecules object.
        if type(molecules) is _Molecules:
            is_sire_container = True

        # A list of Molecule objects.
        elif type(molecules) is list and all(isinstance(x, _Molecule) for x in molecules):
            pass

        # Invalid argument.
        else:
            raise TypeError("'molecules' must be of type 'BioSimSpace._SireWrappers.Molecule' "
                            "or a list of 'BioSimSpace._SireWrappers.Molecule' types.")

        # Remove the molecules in the system.
        if is_sire_container:
            self._sire_object.remove(molecules._sire_object, _SireMol.MGName("all"))
        else:
            for mol in molecules:
                self._sire_object.remove(mol._sire_object.number())

        # Reset the index mappings.
        self._reset_mappings()

        # Update the molecule numbers.
        self._mol_nums = self._sire_object.molNums()

    def removeWaterMolecules(self):
        """Remove all of the water molecules from the system."""

        # Get the list of water molecules.
        waters = self.getWaterMolecules()

        # Remove the molecules in the system.
        self._sire_object.remove(waters._sire_object, _SireMol.MGName("all"))

        # Reset the index mappings.
        self._reset_mappings()

        # Update the molecule numbers.
        self._mol_nums = self._sire_object.molNums()

    def updateMolecule(self, index, molecule):
        """Updated the molecule at the given index.

           Parameters
           ----------

           index : int
               The index of the molecule.

           molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
               The updated (or replacement) molecule.
        """

        if type(index) is not int:
            raise TypeError("'index' must be of type 'int'")

        if index < -self.nMolecules() or index >= self.nMolecules():
            raise IndexError("The molecule 'index' is out of range.")

        if type(molecule) is not _Molecule:
            raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'.")


        # The molecule numbers match.
        if self[index].number() == molecule.number():
            self.updateMolecules(molecule)

        # The molecule numbers don't match.
        else:
            # Store the original molecule number.
            mol_num = self[index]._sire_object.number()

            # Store a copy of the passed molecule.
            mol_copy = molecule._sire_object.__deepcopy__()

            # Delete the existing molecule.
            self._sire_object.remove(mol_num)

            # Renumber the new molecule.
            edit_mol = mol_copy.edit()
            edit_mol.renumber(mol_num)
            mol_copy = edit_mol.commit()

            # Add the renumbered molecule to the system.
            self._sire_object.add(mol_copy, _SireMol.MGName("all"))

            # Reset the index mappings.
            self._reset_mappings()

            # Update the molecule numbers.
            self._mol_nums = self._sire_object.molNums()

    def updateMolecules(self, molecules):
        """Update a molecule, or list of molecules in the system.

           Parameters
           ----------

           molecules : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                       [:class:`Molecule <BioSimSpace._SireWrappers.Molecule>`]
              A Molecule, or list of Molecule objects.
        """

        # Convert tuple to a list.
        if type(molecules) is tuple:
            molecules = list(molecules)

        # A Molecule object.
        if type(molecules) is _Molecule:
            molecules = [molecules]

        # A list of Molecule objects.
        elif type(molecules) is list and all(isinstance(x, _Molecule) for x in molecules):
            pass

        # Invalid argument.
        else:
            raise TypeError("'molecules' must be of type 'BioSimSpace._SireWrappers.Molecule' "
                            "or a list of 'BioSimSpace._SireWrappers.Molecule' types.")

        # Update each of the molecules.
        # TODO: Currently the Sire.System.update method doesn't work correctly
        # for certain changes to the Molecule molInfo object. As such, we remove
        # the old molecule from the system, then add the new one in.
        for mol in molecules:
            try:
                self._sire_object.update(mol._sire_object)
            except:
                self._sire_object.remove(mol._sire_object.number())
                self._sire_object.add(mol._sire_object, _SireMol.MGName("all"))

        # Reset the index mappings.
        self._reset_mappings()

        # Update the molecule numbers.
        self._mol_nums = self._sire_object.molNums()

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

    def getMolecules(self, group="all"):
        """Return a list containing all of the molecules in the specified group.

           Parameters
           ----------

           group : str
               The name of the molecule group.

           Returns
           -------

           molecules : [:class:`Molecule <BioSimSpace._SireWrappers.Molecule>`]
               The list of molecules in the group.
        """

        if type(group) is not str:
            raise TypeError("'group' must be of type 'str'")

        # Try to extract the molecule group.
        try:
            molgrp = self._sire_object.group(_SireMol.MGName(group))
        except:
            raise ValueError("No molecules in group '%s'" % group)

        # Return a molecules container.
        return _Molecules(molgrp)

    def getWaterMolecules(self):
        """Return a list containing all of the water molecules in the system.

           Returns
           -------

           molecules : :class:`Molecules <BioSimSpace._SireWrappers.Molecules>`
               A container of water molecule objects.
        """

        return _Molecules(self._sire_object.search("water").toGroup())

    def nWaterMolecules(self):
        """Return the number of water molecules in the system.

           Returns
           -------

           num_waters : int
               The number of water molecules in the system.
        """
        return len(self.getWaterMolecules())

    def getPerturbableMolecules(self):
        """Return a list containing all of the perturbable molecules in the system.

           Returns
           -------

           molecules : [:class:`Molecule <BioSimSpace._SireWrappers.Molecule>`]
               A list of perturbable molecules.
        """

        molecules = []

        for mol in self._sire_object.search("perturbable"):
            molecules.append(_Molecule(mol))

        return molecules

    def nPerturbableMolecules(self):
        """Return the number of perturbable molecules in the system.

           Returns
           -------

           num_perturbable : int
               The number of perturbable molecules in the system.
        """
        return len(self.getPerturbableMolecules())

    def search(self, query, property_map={}):
        """Search the system for atoms, residues, and molecules. Search results
           will be reduced to their minimal representation, i.e. a molecule
           containing a single residue will be returned as a residue.

           Parameters
           ----------

           query : str
               The search query.

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }

           Returns
           -------

           results : [:class:`Atom <BioSimSpace._SireWrappers.Atom>`, \
                      :class:`Residue <BioSimSpace._SireWrappers.Residue>`, \
                      :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, ...]
               A list of objects matching the search query.

           Examples
           --------

           Search for residues names ALA.

           >>> result = system.search("resname ALA")

           Search for residues names ALA using a case insensitive search.

           >>> result = system.search("resname /ala/i")

           Search for all carbon atoms in residues named ALA.

           >>> result = system.search("resname ALA and atomname C")

           Search for all oxygen or hydrogen atoms.

           >>> result = system.search("element oxygen or element hydrogen")

           Search for atom index 23 in molecule index 10.

           >>> result = system.search("molidx 10 and atomidx 23")
        """

        if type(query) is not str:
            raise TypeError("'query' must be of type 'str'")

        if type(property_map) is not dict:
            raise TypeError("'property_map' must be of type 'dict'")

        # Initialise a list to hold the search results.
        results = []

        try:
            # Query the Sire system.
            search_result = _SireMol.Select(query)(self._sire_object, property_map)

        except Exception as e:
            msg = "'Invalid search query: %r" % query
            if _isVerbose():
                raise ValueError(msg) from e
            else:
                raise ValueError(msg) from None

        return _SearchResult(search_result)

    def getIndex(self, item):
        """Convert indices of atoms and residues to their absolute values in
           this system.

           Parameters
           ----------

           item : :class:`Atom <BioSimSpace._SireWrappers.Atom>`, \
                  :class:`Residue <BioSimSpace._SireWrappers.Residue>` \
                  :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
               An Atom, Residue, or Molecule object from the System, or a
               list containing objects of these types.

           Returns
           -------

           index : int, [int]
               The absolute index of the atom/residue/molecule in the system.
        """

        # Convert single object to list.
        if type(item) is not tuple and type(item) is not list:
            item = [item]

        # Create a list to hold the indices.
        indices = []

        # Loop over all of the items.
        for x in item:
            # Atom.
            if type(x) is _Atom:
                # Only compute the atom index mapping if it hasn't already
                # been created.
                if len(self._atom_index_tally) == 0:
                    # Loop over all molecules in the system and keep track
                    # of the cumulative number of atoms.
                    num_atoms = 0
                    for num in self._sire_object.molNums():
                        self._atom_index_tally[num] = num_atoms
                        num_atoms += self._sire_object.molecule(num).nAtoms()

                # Get the AtomIdx and MolNum from the Atom.
                index = x.index()
                mol_num = x._sire_object.molecule().number()

                try:
                    index += self._atom_index_tally[mol_num]
                except KeyError:
                    raise KeyError("The atom belongs to molecule '%s' that is not part of "
                                "this system!" % mol_num)

                indices.append(index)

            # Residue.
            elif type(x) is _Residue:
                # Only compute the residue index mapping if it hasn't already
                # been created.
                if len(self._residue_index_tally) == 0:
                    # Loop over all molecules in the system and keep track
                    # of the cumulative number of residues.
                    num_residues = 0
                    for num in self._sire_object.molNums():
                        self._residue_index_tally[num] = num_residues
                        num_residues += self._sire_object.molecule(num).nResidues()

                # Get the ResIdx and MolNum from the Residue.
                index = x.index()
                mol_num = x._sire_object.molecule().number()

                try:
                    index += self._residue_index_tally[mol_num]
                except KeyError:
                    raise KeyError("The residue belongs to molecule '%s' that is not part of "
                                "this system!" % mol_num)

                indices.append(index)

            # Residue.
            elif type(x) is _Molecule:
                # Only compute the molecule index mapping if it hasn't already
                # been created.
                if len(self._molecule_index) == 0:
                    # Loop over all molecules in the system by index and map
                    # the MolNum to index.
                    for idx in range(0, self.nMolecules()):
                        self._molecule_index[self._sire_object[_SireMol.MolIdx(idx)].number()] = idx

                # Get the MolNum from the molecule.
                mol_num = x._sire_object.molecule().number()

                try:
                    index = self._molecule_index[mol_num]
                except KeyError:
                    raise KeyError("The molecule '%s' is not part of this system!" % mol_num)

                indices.append(index)

            # Unsupported.
            else:
                raise TypeError("'item' must be of type 'BioSimSpace._SireWrappers.Atom' "
                                "or 'BioSimSpace._SireWrappers.Residue'")

        # If this was a single object, then return a single index.
        if len(indices) == 1:
            return indices[0]
        # Return the list of indices.
        else:
            return indices

    def setBox(self, box, angles=3*[_Angle(90, "degree")], property_map={}):
        """Set the size of the periodic simulation box.

           Parameters
           ----------

           box : [:class:`Length <BioSimSpace.Types.Length>`]
               The box vector magnitudes in each dimension.

           angles : [:class:`Angle <BioSimSpace.Types.Angle>`]
               The box vector angles: yz, xz, and xy.

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }

        """

        # Convert tuples to lists.
        if type(box) is tuple:
            box = list(box)
        if type(angles) is tuple:
            angles = list(angles)

        # Validate input.
        if type(box) is not list or not all(isinstance(x, _Length) for x in box):
            raise TypeError("'box' must be a list of 'BioSimSpace.Types.Length' objects.")

        if type(angles) is not list or not all(isinstance(x, _Angle) for x in angles):
            raise TypeError("'angles' must be a list of 'BioSimSpace.Types.Angle' objects.")

        if len(box) != 3:
            raise ValueError("'box' must contain three items.")

        if len(box) != 3:
            raise ValueError("'angles' must contain three items.")

        # Convert sizes to Anstrom.
        vec = [x.angstroms().magnitude() for x in box]

        # Whether the box is triclinic.
        is_triclinic = False

        # If any angle isn't 90 degrees, then the box is triclinic.
        for x in angles:
            if abs(x.degrees().magnitude() - 90) > 1e-6:
                is_triclinic = True
                break

        if is_triclinic:
            # Convert angles to Sire units.
            ang = [x.degrees().magnitude()*_SireUnits.degree for x in angles]

            # Create a triclinic box object.
            space = _SireVol.TriclinicBox(vec[0], vec[1], vec[2],
                                        ang[0], ang[1], ang[2])

        else:
            # Create a periodic box object.
            space = _SireVol.PeriodicBox(_SireMaths.Vector(vec))

        # Set the "space" property.
        self._sire_object.setProperty(property_map.get("space", "space"), space)

    def getBox(self, property_map={}):
        """Get the size of the periodic simulation box.

           Parameters
           ----------

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }

           Returns
           -------

           box_size : [:class:`Length <BioSimSpace.Types.Length>`]
               The box vector magnitudes in each dimension.

           angles : [:class:`Angle <BioSimSpace.Types.Angle>`]
               The box vector angles: yz, xz, and xy.
       """

        # Get the "space" property and convert to a list of BioSimSpace.Type.Length
        # objects.
        try:
            space = self._sire_object.property(property_map.get("space", "space"))

            # Periodic box.
            if type(space) is _SireVol.PeriodicBox:
                box = [ _Length(x, "Angstrom") for x in space.dimensions() ]
                angles = 3*[_Angle(90, "degrees")]

            # TriclinicBox box.
            elif type(space) is _SireVol.TriclinicBox:
                box = [_Length(space.vector0().magnitude(), "Angstrom"),
                       _Length(space.vector1().magnitude(), "Angstrom"),
                       _Length(space.vector2().magnitude(), "Angstrom")]
                angles = [_Angle(space.alpha(), "degree"),
                          _Angle(space.beta(), "degree"),
                          _Angle(space.gamma(), "degree")]
        except:
            box = None
            angles = None

        return box, angles

    def translate(self, vector, property_map={}):
        """Translate the system.

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

        # Translate each of the molecules in the system.
        for n in self._sire_object.molNums():
            # Copy the property map.
            _property_map = property_map.copy()

            # If this is a perturbable molecule, use the coordinates at lambda = 0.
            if self._sire_object.molecule(n).hasProperty("is_perturbable"):
                _property_map["coordinates"] = "coordinates0"

            mol = self._sire_object[n].move().translate(_SireMaths.Vector(vec), _property_map).commit()
            self._sire_object.update(mol)

    def getRestraintAtoms(self, restraint, mol_index=None, is_absolute=True,
            allow_zero_matches=False, property_map={}):
        """Get the indices of atoms involved in a restraint.

           Parameters
           ----------

           restraint : str
               The type of restraint.

           mol_index : int
               The index of the molecule of interest. If None, then the
               entire system is searched.

           is_absolute : bool
               Whether the indices are absolute, i.e. indexed within the
               entire system. If False, then indices are relative to the
               molecule in which the atom is found.

           allow_zero_matches : bool
               Whether to raise an exception when no atoms match the restraint.
               Setting to false is useful if you need to determine restraints
               on a per-molecule basis, i.e. so molecules won't contain atoms
               matching the restraint, but the entire system will.

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their

           Returns
           -------

           indices : [int]
               A list of the backbone atom indices.
        """

        if type(restraint) is not str:
            raise TypeError("'restraint' must be of type 'str'.")

        # Allowed keyword options.
        allowed = ["backbone", "heavy", "all"]

        if type(restraint) is str:
            # Convert to lower case and strip whitespace.
            restraint = restraint.lower().replace(" ", "")
            if restraint not in allowed:
                raise ValueError(f"'restraint' must be one of: {allowed}")
        else:
            raise TypeError("'restraint' must be of type 'str'.")

        if mol_index is not None and type(mol_index) is not int:
            raise TypeError("'mol_index' must be of type 'int'.")

        if type(is_absolute) is not bool:
            raise TypeError("'is_absolute' must be of type 'bool'.")

        if type(property_map) is not dict:
            raise TypeError("'property_map' must be of type 'dict'.")

        # Initialise the list of indices.
        indices = []

        # Get the element property.
        element = property_map.get("element", "element")

        # A set of protein residues. Taken from MDAnalysis.
        prot_res = {
            # CHARMM top_all27_prot_lipid.rtf
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HSD",
            "HSE", "HSP", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR",
            "TRP", "TYR", "VAL", "ALAD",
            ## 'CHO','EAM', # -- special formyl and ethanolamine termini of gramicidin
            # PDB
            "HIS", "MSE",
            # from Gromacs 4.5.3 oplsaa.ff/aminoacids.rtp
            "ARGN", "ASPH", "CYS2", "CYSH", "QLN", "PGLU", "GLUH", "HIS1", "HISD",
            "HISE", "HISH", "LYSH",
            # from Gromacs 4.5.3 gromos53a6.ff/aminoacids.rtp
            "ASN1", "CYS1", "HISA", "HISB", "HIS2",
            # from Gromacs 4.5.3 amber03.ff/aminoacids.rtp
            "HID", "HIE", "HIP", "ORN", "DAB", "LYN", "HYP", "CYM", "CYX", "ASH",
            "GLH", "ACE", "NME",
            # from Gromacs 2016.3 amber99sb-star-ildn.ff/aminoacids.rtp
            "NALA", "NGLY", "NSER", "NTHR", "NLEU", "NILE", "NVAL", "NASN", "NGLN",
            "NARG", "NHID", "NHIE", "NHIP", "NTRP", "NPHE", "NTYR", "NGLU", "NASP",
            "NLYS", "NPRO", "NCYS", "NCYX", "NMET", "CALA", "CGLY", "CSER", "CTHR",
            "CLEU", "CILE", "CVAL", "CASF", "CASN", "CGLN", "CARG", "CHID", "CHIE",
            "CHIP", "CTRP", "CPHE", "CTYR", "CGLU", "CASP", "CLYS", "CPRO", "CCYS",
            "CCYX", "CMET", "CME", "ASF",
        }

        # A list of ion elements.
        ions = [
            "F", "Cl", "Br", "I", "Li", "Na", "K", "Rb", "Cs", "Mg", "Tl", "Cu", "Ag",
            "Be", "Cu", "Ni", "Pt", "Zn", "Co", "Pd", "Ag", "Cr", "Fe", "Mg", "V", "Mn",
            "Hg", "Cd", "Yb", "Ca", "Sn", "Pb", "Eu", "Sr", "Sm", "Ba", "Ra", "Al", "Fe",
            "Cr", "In", "Tl", "Y", "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy",
            "Er", "Tm", "Lu", "Hf", "Zr", "Ce", "U", "Pu", "Th",
        ]

        # Search the entire system.
        if mol_index is None:
            # Backbone restraints.
            if restraint == "backbone":
                # Find all N, CA, C, and O atoms in protein residues.
                string = "(atoms in not water) and (atoms in resname " + ",".join(prot_res) + " and atomname N,CA,C,O)"
                search = self.search(string)

            elif restraint == "heavy":
                # Convert to a formatted string for the search.
                ion_string = ",".join(ions + ["H"])
                # Find all non-water, non-hydrogen, non-ion elements.
                string = f"(atoms in not water) and not element {ion_string}"
                search = self.search(string)

            elif restraint == "all":
                # Convert to a formatted string for the search.
                ion_string = ",".join(ions)
                # Find all non-water, non-ion elements.
                string = f"(atoms in not water) and not element {ion_string}"
                search = self.search(string)

        # Search the chosen molecule.
        else:
            if restraint == "backbone":
                mol = self[mol_index]
                if mol.isWater():
                    search = []

                else:
                    # Find all N, CA, C, and O atoms in protein residues.
                    string = "atoms in resname " + ",".join(prot_res) + " and atomname N,CA,C,O"
                    search = mol.search(string)

            elif restraint == "heavy":
                mol = self[mol_index]
                if mol.isWater():
                    search = []
                else:
                    # Convert to a formatted string for the search.
                    ion_string = ",".join(ions + ["H"])
                    # Find all non-water, non-hydrogen, non-ion elements.
                    string = f"not element {ion_string}"
                    search = mol.search(string)

            elif restraint == "all":
                mol = self[mol_index]
                if mol.isWater():
                    search = []
                else:
                    # Convert to a formatted string for the search.
                    ion_string = ",".join(ions)
                    # Find all non-water, non-ion elements.
                    string = f"not element {ion_string}"
                    search = mol.search(string)

        # Raise an exception if no atoms match the restraint.
        if len(search) == 0 and not allow_zero_matches:
            msg = "No atoms matched the restraint!"
            if restraint == "backbone":
                msg += " Backbone restraints only apply to atoms in protein residues."
            raise _IncompatibleError(msg)

        # Now loop over all matching atoms and get their indices.
        for atom in search:
            if is_absolute:
                indices.append(self.getIndex(atom))
            else:
                indices.append(atom.index())

        # The indices should be sorted, but do so just in case.
        indices.sort()

        return indices

    def _isParameterised(self, property_map={}):
        """Whether the system is parameterised, i.e. can we run a simulation
           using this system. Essentially we check whether every molecule in
           the system has "bond" and "LJ" properties and assume that is
           parameterised if so.

           Parameters
           ----------

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their

           Returns
           -------

           is_parameterised : bool
               Whether the system is parameterised.
        """

        # Get the "bond" property.
        bond = property_map.get("bond", "bond")
        # Test for perturbable molecules too.
        bond0 = bond + "0"

        # Get the "LJ" property.
        LJ = property_map.get("LJ", "LJ")
        # Test for perturbable molecules too.
        LJ0 = LJ + "0"

        # Check each molecule for "bond" and "LJ" properties.
        for mol in self.getMolecules():
            props = mol._sire_object.propertyKeys()
            if (bond not in props and bond0 not in props) or \
               (LJ   not in props and LJ0   not in props):
                return False

        # If we get this far, then all molecules are okay.
        return True

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

        # Initialise the coordinates vector.
        coord = []

        # Get all of the molecules in the system.
        mols = self.getMolecules()

        # Loop over all of the molecules.
        for idx, mol in enumerate(mols):

            # Extract the atomic coordinates and append them to the vector.
            try:
                if "coordinates" in property_map:
                    prop = property_map["coordinates"]
                else:
                    if mol.isPerturbable():
                        prop = "coordinates0"
                    else:
                        prop = "coordinates"
                coord.extend(mol._sire_object.property(prop).toVector())

            except UserWarning as e:
                msg = "Unable to compute the axis-aligned bounding " + \
                      "box since a molecule has no 'coordinates' property."
                if _isVerbose():
                    raise _IncompatibleError(msg) from e
                else:
                    raise _IncompatibleError(msg) from None

        # Return the AABox for the coordinates.
        return _SireVol.AABox(coord)

    def _renumberMolecules(self, molecules, is_rebuild=False):
        """Helper function to renumber the molecules to be consistent with the
           system.

           Parameters
           ----------

           molecules : [:class:`Molecule <BioSimSpace._SireWrappers.Molecule>`]
               A list of molecule objects.

           Returns
           -------

           molecules : [:class:`Molecule <BioSimSpace._SireWrappers.Molecule>`]
               The renumber list of molecule objects.
        """

        # Renumber everything.
        if is_rebuild:
            num_molecules = 0
            num_residues = 0
            num_atoms = 0

        # Get the current number of molecules, residues, and atoms.
        else:
            num_molecules = self.nMolecules()
            num_residues = self.nResidues()
            num_atoms = self.nAtoms()

        # Create a list to hold the modified molecules.
        new_molecules = []

        # Loop over all of the molecules.
        for mol in molecules:

            # Create a copy of the molecule.
            new_mol = _Molecule(mol)

            # Get the Sire molecule and make it editable.
            edit_mol = new_mol._sire_object.edit()

            # Renumber the molecule.
            edit_mol = edit_mol.renumber(_SireMol.MolNum(num_molecules+1)).molecule()
            num_molecules += 1

            # A hash mapping between old and new numbers.
            num_hash = {}

            # Loop over all residues and add them to the hash.
            for res in edit_mol.residues():
                num_hash[res.number()] = _SireMol.ResNum(num_residues+1)
                num_residues += 1

            # Renumber the residues.
            edit_mol = edit_mol.renumber(num_hash).molecule()

            # Clear the hash.
            num_hash = {}

            # Loop over all of the atoms and add them to the hash.
            for atom in edit_mol.atoms():
                num_hash[atom.number()] = _SireMol.AtomNum(num_atoms+1)
                num_atoms += 1

            # Renumber the atoms.
            edit_mol = edit_mol.renumber(num_hash).molecule()

            # Commit the changes and replace the molecule.
            new_mol._sire_object = edit_mol.commit()

            # Append to the list of molecules.
            new_molecules.append(new_mol)

        # Return the renumbered molecules.
        return new_molecules

    def _updateCoordinates(self, system, property_map0={}, property_map1={},
            is_lambda1=False):
        """Update the coordinates of atoms in the system.

           Parameters
           ----------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               A system containing the updated coordinates.

           property_map0 : dict
               A dictionary that maps system "properties" to their user defined
               values in this system.

           property_map1 : dict
               A dictionary that maps system "properties" to their user defined
               values in the passed system.

           is_lambda1 : bool
              Whether to update coordinates of perturbed molecules at lambda = 1.
              By default, coordinates at lambda = 0 are used.
        """

        # Validate the system.
        if type(system) is not System:
            raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")

        # Check that the passed system contains the same number of molecules.
        if self.nMolecules() != system.nMolecules():
            raise _IncompatibleError("The passed 'system' contains a different number of "
                                     "molecules. Expected '%d', found '%d'"
                                     % (self.nMolecules(), system.nMolecules()))

        if type(property_map0) is not dict:
            raise TypeError("'property_map0' must be of type 'dict'.")

        if type(property_map1) is not dict:
            raise TypeError("'property_map1' must be of type 'dict'.")

        if type(is_lambda1) is not bool:
            raise TypeError("'is_lambda1' must be of type 'bool'.")

        # Check that each molecule in the system contains the same number of atoms.
        for idx in range(0, self.nMolecules()):
            # Extract the number of atoms in the molecules.
            num_atoms0 = self._sire_object.molecule(_SireMol.MolIdx(idx)).nAtoms()
            num_atoms1 = system._sire_object.molecule(_SireMol.MolIdx(idx)).nAtoms()

            if num_atoms0 != num_atoms1:
                raise _IncompatibleError("Mismatch in atom count for molecule '%d': "
                                         "Expected '%d', found '%d'" % (num_atoms0, num_atoms1))

        # Work out the name of the "coordinates" property.
        prop0 = property_map0.get("coordinates0", "coordinates")
        prop1 = property_map1.get("coordinates1", "coordinates")

        # Loop over all molecules and update the coordinates.
        for idx in range(0, self.nMolecules()):
            # Extract the molecules from each system.
            mol0 = self._sire_object.molecule(_SireMol.MolIdx(idx))
            mol1 = system._sire_object.molecule(_SireMol.MolIdx(idx))

            # Check whether the molecule is perturbable.
            if mol0.hasProperty("is_perturbable"):
                if is_lambda1:
                    prop = "coordinates1"
                else:
                    prop = "coordinates0"
            else:
                prop = prop0

            # Try to update the coordinates property.
            try:
                mol0 = mol0.edit().setProperty(prop, mol1.property(prop1)).molecule().commit()
            except Exception as e:
                msg = "Unable to update 'coordinates' for molecule index '%d'" % idx
                if _isVerbose():
                    raise _IncompatibleError(msg) from e
                else:
                    raise _IncompatibleError(msg) from None

            # Update the molecule in the original system.
            self._sire_object.update(mol0)

    @staticmethod
    def _createSireSystem(molecules):
        """Create a Sire system from a Molecules object or a list of Molecule
           objects.

           Parameters
           ----------

           molecules : :class:`Molecules <BioSimSpace._SireWrappers.Molecules>` \
                       [:class:`Molecule <BioSimSpace._SireWrappers.Molecule>`]

           Returns
           -------

           system : Sire.System.System
               A Sire system object.
        """

        # Create an empty Sire System.
        system = _SireSystem.System("BioSimSpace System")

        # Create a new "all" molecule group.
        molgrp = _SireMol.MoleculeGroup("all")

        # Add the molecules to the group.
        if type(molecules) is _Molecules:
            molgrp.add(molecules._sire_object)
        else:
            for mol in molecules:
                molgrp.add(mol._sire_object)

        # Add the molecule group to the system.
        system.add(molgrp)

        return system

    def _getRelativeIndices(self, abs_index):
        """Given an absolute index, get the number of the molecule to which it
           belongs and the relative index within that molecule.

           Parameters
           ----------

           abs_index : int
               The absolute index of the atom in the system.

           Returns
           -------

           mol_index : int
               The molecule index to which the atom belongs.

           rel_index : int
               The relative index of the atom in the molecule to which it
               belongs.
        """
        # Make sure the atom index tally has been created.
        self.getIndex(self[0].getAtoms()[0])

        # Set tally counters for the total number of atoms to date
        # and the total up to the previous molecule.
        tally = 0
        tally_last = 0

        # Set the current and previous MolNum.
        mol_num = self._mol_nums[0]
        mol_num_last = self._mol_nums[0]

        # Loop over each molecule until the total atom number is equal
        # to or greater than the absolute index. When it is, return the
        # absolute index minus the tally up to the previous molecule.
        for num, num_atoms in self._atom_index_tally.items():
            if abs_index >= tally:
                tally_last = tally
                tally = num_atoms
                mol_num_last = mol_num
                mol_num = num
            else:
                mol_idx = self._molecule_index[mol_num_last]
                return mol_idx, abs_index - tally_last

        raise ValueError("'abs_index' exceeded system atom tally!")

    def _reset_mappings(self):
        """Internal function to reset index mapping dictionaries."""

        # Clear dictionaries.
        self._molecule_index = {}
        self._atom_index_tally = {}
        self._residue_index_tally = {}

        # Rebuild the MolNum to index mapping.
        for idx in range(0, self.nMolecules()):
            self._molecule_index[self._sire_object[_SireMol.MolIdx(idx)].number()] = idx

    def _set_water_topology(self, format, property_map={}):
        """Internal function to swap the water topology to AMBER or GROMACS format.

            Parameters
            ----------

            format : string
                The format to convert to: either "AMBER" or "GROMACS".

            property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Validate input.

        if type(format) is not str:
            raise TypeError("'format' must be of type 'str'")

        if type(property_map) is not dict:
            raise TypeError("'property_map' must be of type 'dict'")

        # Strip whitespace and convert to upper case.
        format = format.replace(" ", "").upper()

        # We allow conversion to AMBER or GROMACS format water toplogies. While
        # AMBER requires a specific topology, GROMACS can handle whatever. We
        # leave this format as an option in case we want to allow the user to
        # swap the topology in future.
        if format not in ["AMBER", "GROMACS"]:
            raise ValueError("'format' must be 'AMBER' or 'GROMACS'.")

        # Get the water molecules.
        waters = self.getWaterMolecules()

        if len(waters) > 0:

            # Don't perform conversion if the topology already matches the
            # the template for the desired format.

            if format == "AMBER":
                if waters[0].isAmberWater():
                    return

            elif format == "GROMACS":
                if waters[0].isGromacsWater():
                    return

            # There will be a "water_model" system property if this object was
            # solvated by BioSimSpace.
            if "water_model" in self._sire_object.propertyKeys():
                water_model = self._sire_object.property("water_model").toString()

            # Otherwise, convert to an appropriate topology.
            else:
                num_point = waters[0].nAtoms()

                if num_point == 3:
                    # TODO: Assume TIP3P. Not sure how to detect SPC/E at present.
                    water_model = "TIP3P"
                elif num_point == 4:
                    water_model = "TIP4P"
                elif num_point == 5:
                    water_model = "TIP5P"

            if format == "AMBER":
                new_waters = _SireIO.setAmberWater(self._sire_object.search("water"), water_model)
            else:
                new_waters = _SireIO.setGromacsWater(self._sire_object.search("water"), water_model)

            # Create a new system and molecule group.
            system = _SireSystem.System("BioSimSpace System")
            molgrp = _SireMol.MoleculeGroup("all")

            # Create a dictionary of water molecule numbers to list index.
            wat_num_to_idx = {}
            for idx, wat in enumerate(waters):
                wat_num_to_idx[wat.number()] = idx

            # Add all moleclules from the original system into the new one,
            # preserving the order.
            for mol in self.getMolecules():
                try:
                    # Water molecule.
                    molgrp.add(new_waters[wat_num_to_idx[mol.number()]])
                except:
                    # Non-water molecule.
                    molgrp.add(mol._sire_object)

            # Add the group to the system.
            system.add(molgrp)

            # Copy across system properties.
            for prop in self._sire_object.propertyKeys():
                system.setProperty(prop, self._sire_object.property(prop))

            # Set the system.
            self._sire_object = system

            # Reset the index mappings.
            self._reset_mappings()

            # Update the molecule numbers.
            self._mol_nums = self._sire_object.molNums()

# Import at bottom of module to avoid circular dependency.
from ._atom import Atom as _Atom
from ._molecule import Molecule as _Molecule
from ._molecules import Molecules as _Molecules
from ._residue import Residue as _Residue
from ._search_result import SearchResult as _SearchResult
