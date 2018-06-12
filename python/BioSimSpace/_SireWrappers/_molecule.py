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
A thin wrapper around Sire.Mol. This is an internal package and should
not be directly exposed to the user.

Author: Lester Hedges <lester.hedges@gmail.com>
"""

import Sire.Mol as _SireMol
import Sire.MM as _SireMM

__all__ = ["Molecule"]

class IncompatibleError(Exception):
    pass

class Molecule():
    """A container class for storing a molecule."""

    def __init__(self, molecule):
        """Constructor.

           Positional arguments:

           molecule -- A Sire Molecule object.
        """

        # Check that the molecule is valid.
        if not isinstance(molecule, _SireMol.Molecule):
            raise TypeError("'molecule' must be of type 'Sire.Mol._Mol.Molecule'")

        # Set the molecule.
        self._molecule = molecule

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.Molecule: nAtoms=%d, nResidues=%d>" % (self.nAtoms(), self.nResidues())

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return "<BioSimSpace.Molecule: nAtoms=%d, nResidues=%d>" % (self.nAtoms(), self.nResidues())

    def nAtoms(self):
        """Return the number of atoms in the molecule."""
        return self._molecule.nAtoms()

    def nResidues(self):
        """Return the number of residues in the molecule."""
        return self._molecule.nResidues()

    def _getSireMolecule(self):
        """Return the full Sire Molecule object."""
        return self._molecule

    def _makeCompatibleWith(self, molecule, overwrite=False):
        """Make this molecule compatible with this one, i.e. match atoms and
           add all additional properties.

           Positional arguments:

           molecule  -- The molecule to match with.

           Keyword arguments:

           overwrite -- Whether to overwrite any duplicate properties.
        """

        # Check that the molecule is valid.
        if isinstance(molecule, _SireMol.Molecule):
            mol1 = molecule
        if type(molecule) is Molecule:
            mol1 = molecule._molecule
        else:
            raise TypeError("'molecule' must be of type 'BioSimSpace.Molecule', or 'Sire.Mol._Mol.Molecule'")

        # Get the two Sire molecules.
        mol0 = self._molecule

        # Store the number of atoms to match.
        num_atoms = mol0.nAtoms()

        # The new molecule must have at least as many atoms.
        if mol1.nAtoms() < num_atoms:
            raise IncompatibleError("The passed molecule does not contain enough atoms!")

        # Match the atoms based on residue index and atom name.
        matches = _SireMol.ResIdxAtomNameMatcher().match(mol0, mol1)

        # Have we matched all of the atoms?
        if len(matches) < num_atoms:
            # Atom names might have changed. Try to match using residue index
            # and maximum common substructure (matching light atoms too).
            matches = _SireMol.ResIdxAtomMCS(True, False).match(mol0, mol1)

            # Have we matched all of the atoms?
            if len(matches) < num_atoms:
                raise IncompatibleError("Failed to match all atoms!")

        # Get a list of the property keys for each molecule.
        props0 = mol0.propertyKeys()
        props1 = mol1.propertyKeys()

        # Create a dictionary to flag whether a property has been seen.
        seen_prop = {}
        for prop in props1:
            seen_prop[prop] = False

        # Make the molecule editable.
        edit_mol = mol0.edit()

        # First set atom based properties.

        # Loop over all of the keys in the new molecule.
        for prop in props1:
            # This is a new property, or we are allowed to overwrite.
            if (not mol0.hasProperty(prop)) or overwrite:
                # Loop over all of the atom mapping pairs and set the property.
                for idx0, idx1 in matches.items():
                    # Does the atom have this property?
                    # If so, add it to the matching atom in this molecule.
                    if mol1.atom(idx1).hasProperty(prop):
                        try:
                            edit_mol = edit_mol.atom(idx0).setProperty(prop, mol1.atom(idx1).property(prop)).molecule()
                            seen_prop[prop] = True
                        except:
                            raise IncompatibleError("Failed to map property '%s' from AtomIdx(%d) to AtomIdx(%d)."
                                % (prop, idx1, idx0))

        # Now deal with all unseen properties. These will be non atom-based
        # properties, such as TwoAtomFunctions, StringProperty, etc.
        for prop in seen_prop:
            # Skip parameters and intrascale, since these contain an implicit
            # reference to a molinfo object.
            if not seen_prop[prop] and prop != "parameters":
                # This is a new property, or we are allowed to overwrite.
                if (not mol0.hasProperty(prop)) or overwrite:
                    # Try to match with atoms in the original molecule.
                    try:
                        edit_mol.setProperty(prop, mol1.property(prop).makeCompatibleWith(mol0, matches))

                    # This is probably just a molecule property, such as a forcefield definition.
                    except AttributeError:
                        try:
                            edit_mol.setProperty(prop, mol1.property(prop))
                        except:
                            raise IncompatibleError("Failed to set property '%s'" % prop)

        # Commit the changes.
        self._molecule = edit_mol.commit()
