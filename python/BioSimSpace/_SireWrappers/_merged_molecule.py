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
Functionality for creating and storing merged molecules for use in
dual-topology free energy calculations.

Author: Lester Hedges <lester.hedges@gmail.com>
"""

import Sire.Maths as _SireMaths
import Sire.Mol as _SireMol
import Sire.Vol as _SireVol

from ..Types import Length as _Length

__all__ = ["MergedMolecule"]

class MergedMolecule():
    """A container class for storing merged molecules."""

    def __init__(self, molecule0, molecule1, mapping):
        """Constructor.

           Positional arguments
           --------------------

           molecule0 : Sire.Mol.Molecule
               The initial molecule.

           molecule1 : Sire.Mol.Molecule
               The final molecule.

           mapping : dict
               The mapping between matching atom indices in the two molecules.
        """

        # Check that molecule0 is valid.

        # A Sire Molecule object.
        if type(molecule0) is _SireMol.Molecule:
            self._sire_molecule0 = molecule0.__deepcopy__()

        # Another BioSimSpace Molecule object.
        elif type(molecule0) is _Molecule:
            self._sire_molecule0 = molecule0._sire_molecule.__deepcopy__()

        # Invalid type.
        else:
            raise TypeError("'molecule0' must be of type 'Sire.Mol._Mol.Molecule' "
                + "or 'BioSimSpace._SireWrappers.Molecule'.")

        # Check that molecule1 is valid.

        # A Sire Molecule object.
        if type(molecule1) is _SireMol.Molecule:
            self._sire_molecule1 = molecule1.__deepcopy__()

        # Another BioSimSpace Molecule object.
        elif type(molecule1) is _Molecule:
            self._sire_molecule1 = molecule1._sire_molecule.__deepcopy__()

        # Invalid type.
        else:
            raise TypeError("'molecule1' must be of type 'Sire.Mol._Mol.Molecule' "
                + "or 'BioSimSpace._SireWrappers.Molecule'.")

        # Check that the mapping is valid.
        if type(mapping) is dict:
            for idx0, idx1 in mapping.items():
                if type(idx0) is not _SireMol.AtomIdx or type(idx1) is not _SireMol.AtomIdx:
                        raise TypeError("'mapping' dictionary key:value pairs must be of type 'Sire.Mol.AtomIdx'")
            self._mapping = mapping.copy()
        else:
            raise TypeError("'mapping' must be of type 'dict'.")

        # Create the merged molecule.
        self._sire_molecule = self._merge()

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.MergedMolecule: nAtoms=%d, nResidues=%d>" % (self.nAtoms(), self.nResidues())

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return "<BioSimSpace.MergedMolecule: nAtoms=%d, nResidues=%d>" % (self.nAtoms(), self.nResidues())

    def nAtoms(self):
        """Return the number of atoms in the molecule."""
        return self._sire_molecule.nAtoms()

    def nResidues(self):
        """Return the number of residues in the molecule."""
        return self._sire_molecule.nResidues()

    def nChains(self):
        """Return the number of chains in the molecule."""
        return self._sire_molecule.nChains()

    def translate(self, vector):
        """Translate the molecule.

           Positional arguments
           --------------------

           vector : list, tuple
               The translation vector (in Angstroms).
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

        # Perform the translation.
        self._sire_molecule = self._sire_molecule.move().translate(_SireMaths.Vector(vec)).commit()

    def _getSireMolecule(self):
        """Return the full Sire Molecule object."""
        return self._sire_molecule

    def _merge(self):
        """Create the merged molecule."""

        # Get the atom indices from the mapping.
        idx0 = self._mapping.keys()
        idx1 = self._mapping.values()

        # Create lists to store the atoms that are unique to each molecule.
        atoms0 = []
        atoms1 = []

        # Loop over each molecule to find the unique atom indices.

        # molecule0
        for atom in self._sire_molecule0.atoms():
            if atom.index() not in idx0:
                atoms0.append(atom)

        # molecule1
        for atom in self._sire_molecule1.atoms():
            if atom.index() not in idx1:
                atoms1.append(atom)

        # Create a new molecule to hold the merged molecule.
        self._molecule = _SireMol.Molecule()

        # Add a single residue called RES.
        res = self._molecule.edit().add(_SireMol.ResNum(1))
        res.rename(_SireMol.ResName("RES"))

        # Create a single cut-group.
        cg = res.molecule().add(_SireMol.CGName("1"))

        # Counter for the number of atoms.
        num = 1

        # First add all of the atoms from molecule0.
        for atom in self._sire_molecule0.atoms():
            added = cg.add(atom.name())
            added.renumber(_SireMol.AtomNum(num))
            added.reparent(_SireMol.ResIdx(0))
            num += 1

        # Now add all of the atoms from molecule1 that aren't in molecule0.
        for atom in atoms1:
            added = cg.add(atom.name())
            added.renumber(_SireMol.AtomNum(num))
            added.reparent(_SireMol.ResIdx(0))
            num += 1

        # Return the merged molecule.
        return cg.molecule().commit()

# Import at bottom of module to avoid circular dependency.
from ._molecule import Molecule as _Molecule
