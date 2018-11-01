######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2018
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
A thin wrapper around Sire.System. This is an internal package and should
not be directly exposed to the user.

Author: Lester Hedges <lester.hedges@gmail.com>
"""

import Sire.Maths as _SireMaths
import Sire.Mol as _SireMol
import Sire.System as _SireSystem
import Sire.Vol as _SireVol

from ..Types import Length as _Length

import BioSimSpace.Units as _Units

__all__ = ["System"]

class _MolWithResName(_SireMol.MolWithResID):
    def __init__(self, resname):
        super().__init__(_SireMol.ResName(resname))

class System():
    """A container class for storing molecular systems."""

    def __init__(self, system):
        """Constructor.


           Positional arguments
           --------------------

           system : Sire.System.System, BioSimSpace._SireWrappers.System, Sire.Mol.Molecule,
                    BioSimSpace._SireWrappers.Molecule, [ BioSimSpace._SireWrappers.Molecule ]
               A Sire System object.
        """

        # Check that the system is valid.

        # Convert tuple to a list.
        if type(system) is tuple:
            system = list(system)

        # A Sire System object.
        if type(system) is _SireSystem.System:
            self._sire_system = system.__deepcopy__()

        # Another BioSimSpace System object.
        elif type(system) is System:
            self._sire_system = system._sire_system.__deepcopy__()

        # A Sire Molecule object.
        elif type(system) is _SireMol.Molecule:
            self._sire_system = _SireSystem.System("BioSimSpace System.")
            self.addMolecules(_Molecule(system))

        # A BioSimSpace Molecule object.
        elif type(system) is _Molecule:
            self._sire_system = _SireSystem.System("BioSimSpace System.")
            self.addMolecules(system)

        # A list of BioSimSpace Molecule objects.
        elif type(system) is list:
            if not all(isinstance(x, _Molecule) for x in system):
                raise TypeError("'system' must contain a list of 'BioSimSpace._SireWrappers.Molecule' types.")
            else:
                self._sire_system = _SireSystem.System("BioSimSpace System.")
                self.addMolecules(system)

        # Invalid type.
        else:
            raise TypeError("'system' must be of type 'Sire.System._System.System', 'BioSimSpace._SireWrappers.System', "
                            " Sire.Mol._Mol.Molecule', 'BioSimSpace._SireWrappers.Molecule', "
                            "or a list of 'BioSimSpace._SireWrappers.Molecule' types.")

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.System: nMolecules=%d>" % self.nMolecules()

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return "<BioSimSpace.System: nMolecules=%d>" % self.nMolecules()

    def __add__(self, other):
        """Addition operator."""

        # Create a copy of the current system.
        system = System(self._sire_system.__deepcopy__())

        # Add the new molecules.
        if type(other) is System:
            system.addMolecules(other.getMolecules())
        else:
            system.addMolecules(other)

        # Return the combined system.
        return system

    def __sub__(self, other):
        """Subtraction operator."""

        # Create a copy of the current system.
        system = System(self._sire_system.__deepcopy__())

        # Remove the molecules from the other system.
        if type(other) is System:
            system.removeMolecules(other.getMolecules())
        else:
            system.removeMolecules(other)

        # Return the new system.
        return system

    def copy(self):
        """Create a copy of this system."""
        return System(self)

    def nMolecules(self):
        """Return the number of molecules in the system."""
        return self._sire_system.nMolecules()

    def nResidues(self):
        """Return the number of residues in the system."""

        tally = 0

        for n in self._sire_system.molNums():
            tally += self._sire_system[n].nResidues()

        return tally

    def nChains(self):
        """Return the number of chains in the system."""

        tally = 0

        for n in self._sire_system.molNums():
            tally += self._sire_system[n].nChains()

        return tally

    def nAtoms(self):
        """Return the number of atoms in the system."""

        tally = 0

        for n in self._sire_system.molNums():
            tally += self._sire_system[n].nAtoms()

        return tally

    def charge(self, property_map={}, is_lambda1=False):
        """Return the total molecular charge.


           Keyword arguments
           -----------------

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }

           is_lambda1 : bool
              Whether to use the charge at lambda = 1 if the molecule is merged.
        """

        # Zero the charge.
        charge = 0 * _Units.Charge.electron_charge

        # Loop over all molecules and add the charge.
        for mol in self.getMolecules():
            # Reset the map.
            _property_map = property_map.copy()

            # If the molecule is merged, then re-map the charge property.
            if mol.isMerged():
                if is_lambda1:
                    _property_map = { "charge" : "charge1" }
                else:
                    _property_map = { "charge" : "charge0" }

            # Add the charge.
            try:
                charge += mol.charge(_property_map)
            except:
                pass

        # Return the total charge.
        return charge

    def fileFormat(self, property_map={}):
        """Return the file formats associated with the system.


           Keyword arguments
           -----------------

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

        if "fileformat" in property_map:
            prop = property_map["fileformat"]
        else:
            prop = "fileformat"

        return self._sire_system.property(prop).value()

    def addMolecules(self, molecules):
        """Add a molecule, or list of molecules to the system.


           Positional arguments
           --------------------

           molecules : BioSimSpace._SireWrappers.Molecule, [ BioSimSpace._SireWrappers.Molecule ]
              A Molecule, or list of Molecule objects.
        """

        # Convert tuple to a list.
        if type(molecules) is tuple:
            system = list(molecules)

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

        # The system is empty: create a new Sire system from the molecules.
        if self._sire_system.nMolecules() == 0:
            self._sire_system = self._createSireSystem(molecules)

        # Otherwise, add the molecules to the existing "all" group.
        else:
            for mol in molecules:
                self._sire_system.add(mol._sire_molecule, _SireMol.MGName("all"))

    def removeMolecules(self, molecules):
        """Remove a molecule, or list of molecules from the system.


           Positional arguments
           --------------------

           molecules : BioSimSpace._SireWrappers.Molecule, [ BioSimSpace._SireWrappers.Molecule ]
              A Molecule, or list of Molecule objects.
        """

        # Convert tuple to a list.
        if type(molecules) is tuple:
            system = list(molecules)

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

        # Remove the molecules in the system.
        for mol in molecules:
            self._sire_system.remove(mol._sire_molecule.number())

    def removeWaterMolecules(self):
        """Remove all of the water molecules from the system."""

        # Get the list of water molecules.
        waters = self.getWaterMolecules()

        # Remove the molecules in the system.
        for mol in waters:
            self._sire_system.remove(mol._sire_molecule.number())

    def updateMolecules(self, molecules):
        """Update a molecule, or list of molecules in the system.


           Positional arguments
           --------------------

           molecules : BioSimSpace._SireWrappers.Molecule, [ BioSimSpace._SireWrappers.Molecule ]
              A Molecule, or list of Molecule objects.
        """

        # Convert tuple to a list.
        if type(molecules) is tuple:
            system = list(molecules)

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
            self._sire_system.remove(mol._sire_molecule.number())
            self._sire_system.add(mol._sire_molecule, _SireMol.MGName("all"))

        # Create a new Sire system from the molecules.
        self._sire_system = self._createSireSystem(molecules)

    def getMolecules(self, group="all"):
        """Return a list containing all of the molecules in the specified group.


           Keyword arguments
           -----------------

           group : str
               The name of the molecule group.


           Returns
           -------

           molecule : [ BioSimSpace._SireWrappers.Molecule ]
               The list of molecules in the group.
        """

        if type(group) is not str:
            raise TypeError("'group' must be of type 'str'")

        # Try to extract the molecule group.
        try:
            molgrp = self._sire_system.group(_SireMol.MGName(group))
        except:
            raise ValueError("No molecules in group '%s'" % group)

        # Create a list to store the molecules.
        mols = []

        # Get a list of the MolNums in the group and sort them.
        nums = molgrp.molNums()

        # Loop over all of the molecules in the group and append to the list.
        for num in nums:
            mols.append(_Molecule(molgrp.molecule(num)))

            # This is a merged molecule.
            if mols[-1]._sire_molecule.hasProperty("is_perturbable"):
                mols[-1]._is_merged = True

        return mols

    def getWaterMolecules(self):
        """Return a list containing all of the water molecules in the system.

           Returns
           -------

           waters : [ BioSimSpace._SireWrappers.Molecule ]
               A list of water molecule objects.
        """

        waters = []

        for mol in self._sire_system.search("water"):
            waters.append(_Molecule(mol))

        return waters

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

           molecules : [ BioSimSpace._SireWrappers.Molecule ]
               A list of perturbable molecules.
        """

        molecules = []

        for mol in self._sire_system.search("perturbable"):
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

    def getMolWithResName(self, resname):
        """Return the molecule containing the given residue.


           Positional arguments
           --------------------

           resname : str
               The name of a residue unique to the molecule.


           Returns
           -------

           molecule : BioSimSpace._SireWrappers.Molecule
               The matching molecule.
        """
        try:
            return _Molecule(self._sire_system[_MolWithResName(resname)])
        except:
            raise KeyError("System does not contain residue '%s'" % resname)

    def translate(self, vector, property_map={}):
        """Translate the system.


           Positional arguments
           --------------------

           vector : list, tuple
               The translation vector (in Angstroms).


           Keyword arguments
           -----------------

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
        for n in self._sire_system.molNums():
            # Copy the property map.
            _property_map = property_map.copy()

            # If this is a perturbable molecule, use the coordinates at lambda = 0.
            if self._sire_system.molecule(n).hasProperty("is_perturbable"):
                _property_map["coordinates"] = "coordinates0"

            mol = self._sire_system[n].move().translate(_SireMaths.Vector(vec), _property_map).commit()
            self._sire_system.update(mol)

    def _getSireSystem(self):
        """Return the full Sire System object."""
        return self._sire_system

    def _getBoxSize(self, property_map={}):
        """Get the size of the periodic box.


           Keyword arguments
           -----------------

           system : Sire.System.System
               A Sire molecular system.

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }


           Returns
           -------

           box_size : tuple
               The size of the box in each dimension (in Angstroms).
        """

        try:
            if "space" in property_map:
                prop = property_map["space"]
            else:
                prop = "space"
            box = self._sire_system.property(prop)
            return box.dimensions()

        except UserWarning:
            return None

    def _getAABox(self, property_map={}):
        """Get the axis-aligned bounding box for the molecular system.


           Keyword arguments
           -----------------

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
                    if mol.isMerged():
                        prop = "coordinates0"
                    else:
                        prop = "coordinates"
                coord.extend(mol._sire_molecule.property(prop).toVector())

            except UserWarning:
                raise UserWarning("Molecule %d has no 'coordinates' property." % idx) from None

        # Return the AABox for the coordinates.
        return _SireVol.AABox(coord)

    def _renumberMolecules(self, molecules, is_rebuild=False):
        """Helper function to renumber the molecules to be consistent with the
           system.


           Positional arguments
           --------------------

           molecules : [ BioSimSpace._SireWrappers.Molecule ]
               A list of molecule objects.


           Returns
           -------

           molecules : [ BioSimSpace._SireWrappers.Molecule ]
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
            edit_mol = new_mol._sire_molecule.edit()

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
            new_mol._sire_molecule = edit_mol.commit()

            # Append to the list of molecules.
            new_molecules.append(new_mol)

        # Return the renumbered molecules.
        return new_molecules

    @staticmethod
    def _createSireSystem(molecules):
        """Create a Sire system from a list of molecules.


           Positional arguments
           --------------------

           molecules : [ BioSimSpace._SireWrappers.Molecule ]


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
        for mol in molecules:
            molgrp.add(mol._sire_molecule)

        # Add the molecule group to the system.
        system.add(molgrp)

        return system

# Import at bottom of module to avoid circular dependency.
from ._molecule import Molecule as _Molecule
