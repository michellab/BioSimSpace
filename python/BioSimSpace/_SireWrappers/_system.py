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
import Sire.Move as _SireMove
import Sire.System as _SireSystem
import Sire.Vol as _SireVol

from ..Types import Length as _Length

__all__ = ["System"]

class _MolWithResName(_SireMol.MolWithResID):
    def __init__(self, resname):
        super().__init__(_SireMol.ResName(resname))

class System():
    """A container class for storing molecular systems."""

    def __init__(self, system):
        """Constructor.

           Positional arguments:

           system -- A Sire System object.
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

        # A BioSimSpace Molecule object.
        elif type(system) is _Molecule:
            self._sire_system = _SireSystem.System("BioSimSpace System.")
            self.addMolecules(system)

        # A list of BioSimSpace Molecule objects.
        elif type(system) is list:
            if not all(isinstance(x, _Molecule) for x in system):
                raise TypeError("'system' must contain a of 'BioSimSpace.Molecule' types.")
            else:
                self._sire_system = _SireSystem.System("BioSimSpace System.")
                self.addMolecules(system)

        # Invalid type.
        else:
            raise TypeError("'system' must be of type 'Sire.System._System.System' "
                + "or a list of 'BioSimSpace.Molecule' types.")

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

    def nMolecules(self):
        """Return the number of molecules in the system."""
        return self._sire_system.nMolecules()

    def fileFormat(self, map={}):
        """Return the file formats associated with the system.

           map -- A dictionary that maps system "properties" to their user defined
                  values. This allows the user to refer to properties with their
                  own naming scheme, e.g. { "charge" : "my-charge" }
        """

        if type(map) is not dict:
            raise TypeError("'map' must be of type 'dict'")

        if "fileformat" in map:
            prop = map["fileformat"]
        else:
            prop = "fileformat"

        return self._sire_system.property(prop).value()

    def addMolecules(self, molecules):
        """Add a molecule, or list of molecules to the system.

           Positional arguments:

           molecules -- A Molecule, or list of Molecule objects.
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
            raise TypeError("'molecules' must be of type 'BioSimSpace.Molecule' "
                + "or a list of 'BioSimSpace.Molecule' types.")

        if self._sire_system.nMolecules() == 0:
            # Create a new "all" molecule group.
            molgrp = _SireMol.MoleculeGroup("all")

            # Add the molecules to the group.
            for num_mols, mol in enumerate(molecules):
                # Get the existing molecule.
                m = mol._getSireMolecule()

                # Renumber the molecule, starting from 1.
                m = m.edit().renumber(_SireMol.MolNum(num_mols+1)).molecule().commit()

                # Add the molecule.
                molgrp.add(m)

            # Add the molecule group to the system.
            self._sire_system.add(molgrp)

        # Otherwise, add the molecules to the existing "all" group.
        else:
            for mol in molecules:
                self._sire_system.add(mol._getSireMolecule(), _SireMol.MGName("all"))

    def removeMolecules(self, molecules):
        """Remove a molecule, or list of molecules from the system.

           Positional arguments:

           molecules -- A Molecule, or list of Molecule objects.
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
            raise TypeError("'molecules' must be of type 'BioSimSpace.Molecule' "
                + "or a list of 'BioSimSpace.Molecule' types.")

        # Remove the molecules to the system.
        for mol in molecules:
            self._sire_system.remove(mol._getSireMolecule().number())

    def udpateMolecules(self, molecules):
        """Update a molecule, or list of molecules in the system.

           Positional arguments:

           molecules -- A Molecule, or list of Molecule objects.
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
            raise TypeError("'molecules' must be of type 'BioSimSpace.Molecule' "
                + "or a list of 'BioSimSpace.Molecule' types.")

        # Update each of the molecules.
        for mol in molecule:
            self._sire_system.update(mol._getSireMolecule())

    def getMolecules(self, group="all"):
        """Return a list containing all of the molecules in the specified group.

           Keyword arguments:

           group -- The name of the molecule group.
        """

        if type(group) is not str:
            raise TypeError("'group' must be of type 'str'")

        # Try to extract the molecule group.
        try:
            molgrp = self._sire_system.group(_SireMol.MGName(group)).molecules()
        except:
            raise ValueError("No molecules in group '%s'" % group)

        # Create a list to store the molecules.
        mols = []

        # Get a list of the MolNums in the group and sort them.
        nums = molgrp.molNums()
        nums.sort()

        # Loop over all of the molecules in the group and append to the list.
        for num in nums:
            mols.append(_Molecule(molgrp[num]))

        return mols

    def getMolWithResName(self, resname):
        """Return the molecule containing the given residue.

           Positional arguments:

           resname -- The name of a residue unique to the molecule.
        """
        try:
            return _Molecule(self._sire_system[_MolWithResName(resname)])
        except:
            raise KeyError("System does not contain residue '%s'" % resname)

    def translate(self, vector):
        """Translate the system.

           Positional arguments:

           vector -- The translation vector (in Angstroms).
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
                        + "'BioSimSpace.Types.Length' types only!")
        else:
            raise TypeError("'vector' must be of type 'list' or 'tuple'")

        # Translate each of the molecules in the system.
        for n in self._sire_system.molNums():
            mol = self._sire_system[n].move().translate(_SireMaths.Vector(vec)).commit()
            self._sire_system.update(mol)

    def _getSireSystem(self):
        """Return the full Sire System object."""
        return self._sire_system

    def _getBoxSize(self, map={}):
        """Get the size of the periodic box.

           Keyword arguments:

           system -- A Sire molecular system.
           map    -- A dictionary that maps system "properties" to their user defined
                     values. This allows the user to refer to properties with their
                     own naming scheme, e.g. { "charge" : "my-charge" }
        """

        try:
            if "space" in map:
                prop = map["space"]
            else:
                prop = "space"
            box = self._sire_system.property(prop)
            return box.dimensions()

        except UserWarning:
            return None

    def _getAABox(self, map={}):
        """Get the axis-aligned bounding box for the molecular system.

           Keyword arguments:

           map -- A dictionary that maps system "properties" to their user defined
                  values. This allows the user to refer to properties with their
                  own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Initialise the coordinates vector.
        coord = []

        # Loop over all of the molecules.
        for n in self._sire_system.molNums():

            # Extract the atomic coordinates and append them to the vector.
            try:
                if "coordinates" in map:
                    prop = map["coordinates"]
                else:
                    prop = "coordinates"
                coord.extend(self._sire_system[n].property(prop).toVector())

            except UserWarning:
                raise

        # Return the AABox for the coordinates.
        return _SireVol.AABox(coord)

# Import at bottom of module to avoid circular dependency.
from ._molecule import Molecule as _Molecule
