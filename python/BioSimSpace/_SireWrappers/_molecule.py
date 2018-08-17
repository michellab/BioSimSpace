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

import Sire.Base as _SireBase
import Sire.Maths as _SireMaths
import Sire.MM as _SireMM
import Sire.Mol as _SireMol
import Sire.System as _SireSystem
import Sire.Units as _SireUnits
import Sire.Vol as _SireVol

from .._Exceptions import IncompatibleError as _IncompatibleError
from ..Types import Length as _Length

import BioSimSpace.Units as _Units

from pytest import approx as _approx

__all__ = ["Molecule"]

class Molecule():
    """A container class for storing a molecule."""

    def __init__(self, molecule):
        """Constructor.

           Positional arguments
           --------------------

           molecule : Sire.Mol.Molecule
               A Sire Molecule object.
        """

        # Set the force field variable. This records the force field with which
        # the molecule has been parameterised, i.e. by BSS.Parameters.
        self._forcefield = None

        # Set the molecule as un-merged.
        self._is_merged = False

        # Set the components of the merged molecule to None.
        self._molecule0 = None
        self._molecule1 = None

        # Check that the molecule is valid.

        # A Sire Molecule object.
        if type(molecule) is _SireMol.Molecule:
            self._sire_molecule = molecule.__deepcopy__()
            if self._sire_molecule.hasProperty("is_perturbable"):
                self._convertFromMergedMolecule()

        # Another BioSimSpace Molecule object.
        elif type(molecule) is Molecule:
            self._sire_molecule = molecule._sire_molecule.__deepcopy__()
            if molecule._molecule0 is not None:
                self._molecule0 = Molecule(molecule._molecule0)
            if molecule._molecule1 is not None:
                self._molecule1 = Molecule(molecule._molecule1)
            self._forcefield = molecule._forcefield
            self._is_merged = molecule._is_merged

        # Invalid type.
        else:
            raise TypeError("'molecule' must be of type 'Sire.Mol._Mol.Molecule' "
                + "or 'BioSimSpace._SireWrappers.Molecule'.")

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.Molecule: nAtoms=%d, nResidues=%d>" % (self.nAtoms(), self.nResidues())

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return "<BioSimSpace.Molecule: nAtoms=%d, nResidues=%d>" % (self.nAtoms(), self.nResidues())

    def __add__(self, other):
        """Addition operator."""

        # Convert tuple to a list.
        if type(other) is tuple:
            other = list(other)

        # Create a list of molecules.
        molecules = [self]

        # Validate the input.

        # A single Molecule object.
        if type(other) is Molecule:
            molecules.append(other)

        # A System object.
        elif type(other) is _System:
            system = _System(other)
            system.addMolecules(self)
            return system

        # A list of Molecule objects.
        elif type(other) is list and all(isinstance(x, Molecule) for x in other):
            molecules.extend(other)

        # Unsupported.
        else:
            raise TypeError("'other' must be of type 'BioSimSpace._SireWrappers.System', "
                + "'BioSimSpace._SireWrappers.Molecule', or a list "
                + "of 'BioSimSpace._SireWrappers.Molecule' types")

        # Create and return a new system.
        return _System(molecules)

    def molecule0(self):
        """Return the component of the merged molecule at lambda = 0.

           Returns
           -------

           molecule : BioSimSpace._SireWrappers.Molecule, None
               The component of the merged molecule at lambda = 0.
               Returns None if this isn't a merged molecule.
        """
        return self._molecule0

    def molecule1(self):
        """Return the component of the merged molecule at lambda = 1.

           Returns
           -------

           molecule : BioSimSpace._SireWrappers.Molecule, None
               The component of the merged molecule at lambda = 1.
               Returns None if this isn't a merged molecule.
        """
        return self._molecule1

    def nAtoms(self):
        """Return the number of atoms in the molecule."""
        return self._sire_molecule.nAtoms()

    def nResidues(self):
        """Return the number of residues in the molecule."""
        return self._sire_molecule.nResidues()

    def nChains(self):
        """Return the number of chains in the molecule."""
        return self._sire_molecule.nChains()

    def isMerged(self):
        """Whether this molecule has been merged with another."""
        return self._is_merged

    def charge(self, map={}, is_lambda1=False):
        """Return the total molecular charge.

           Keyword argument
           ----------------

           map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }

           is_lambda1 : bool
              Whether to use the charge at lambda = 1 if the molecule is merged.
        """

        # Copy the map.
        _map = map

        # This is a merged molecule.
        if self.isMerged():
            if is_lambda1:
                _map = { "charge" : "charge1" }
            else:
                _map = { "charge" : "charge0" }

        # Calculate the charge.
        charge = self._sire_molecule.evaluate().charge(_map).value()

        # Return the charge.
        return charge * _Units.Charge.electron_charge

    def translate(self, vector, map={}):
        """Translate the molecule.

           Positional arguments
           --------------------

           vector : list, tuple
               The translation vector (in Angstroms).


           Keyword argument
           ----------------

           map : dict
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
                        + "'BioSimSpace.Types.Length' types only!")
        else:
            raise TypeError("'vector' must be of type 'list' or 'tuple'")

        # Perform the translation.
        self._sire_molecule = self._sire_molecule.move().translate(_SireMaths.Vector(vec), map).commit()

    def toSystem(self):
        """Convert a single Molecule to a System."""
        return _System(self)

    def _getSireMolecule(self):
        """Return the full Sire Molecule object."""
        return self._sire_molecule

    def _makeCompatibleWith(self, molecule, map={}, overwrite=True,
            rename_atoms=False, verbose=False):
        """Make this molecule compatible with passed one, i.e. match atoms and
           add all additional properties.

           Positional arguments
           --------------------

           molecule : BioSimSpace._SireWrappers.Molecule
               The molecule to match with.


           Keyword arguments
           -----------------

           map : dict
               A map between property names and user supplied names.

           overwrite : bool
               Whether to overwrite any duplicate properties.

           rename_atoms : bool
               Whether to rename atoms if they have changed.

           verbose : bool
               Whether to report status updates to stdout.
        """

        # Validate input.

        if isinstance(molecule, _SireMol.Molecule):
            mol1 = molecule
        elif type(molecule) is Molecule:
            mol1 = molecule._sire_molecule
        else:
            raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule', or 'Sire.Mol._Mol.Molecule'")

        if type(map) is not dict:
            raise TypeError("'map' must be of type 'dict'")

        if type(overwrite) is not bool:
            raise TypeError("'overwrite' must be of type 'bool'")

        if type(rename_atoms) is not bool:
            raise TypeError("'rename_atoms' must be of type 'bool'")

        if type(verbose) is not bool:
            raise TypeError("'verbose' must be of type 'bool'")

        # Get the two Sire molecules.
        mol0 = self._sire_molecule

        # Store the number of atoms to match.
        num_atoms = mol0.nAtoms()

        # The new molecule must have at least as many atoms.
        if mol1.nAtoms() < num_atoms:
            raise _IncompatibleError("The passed molecule does not contain enough atoms!")

        # Whether the atoms have been renamed.
        is_renamed = False

        # Instantiate the default atom matcher (match by residue index and atom name).
        matcher = _SireMol.ResIdxAtomNameMatcher()

        # Match the atoms based on residue index and atom name.
        matches = matcher.match(mol0, mol1)

        # Have we matched all of the atoms?
        if len(matches) < num_atoms:
            # Atom names might have changed. Try to match by residue index
            # and coordinates.
            matcher = _SireMol.ResIdxAtomCoordMatcher()
            matches = matcher.match(mol0, mol1)

            # We need to rename the atoms.
            is_renamed = True

            # Have we matched all of the atoms?
            if len(matches) < num_atoms:
                raise _IncompatibleError("Failed to match all atoms!")

            # Are the atoms in the same order?
            is_reordered = matcher.changesOrder(mol0, mol1)

        else:
            # Are the atoms in the same order?
            is_reordered = matcher.changesOrder(mol0, mol1)

        if verbose:
            print("\nAtom matching successful.\nAtom indices %s reordered." % ("" if is_reordered else "not"))

        # Get a list of the property keys for each molecule.
        props0 = mol0.propertyKeys()
        props1 = mol1.propertyKeys()

        # Copy the property map.
        _map = map.copy()

        # See if any of the new properties are in the map, add them if not.
        for prop in props1:
            if not prop in _map:
                _map[prop] = prop

        # Make the molecule editable.
        edit_mol = mol0.edit()

        if "parameters" in _map:
            param = _map["parameters"]
        else:
            param = "parameters"

        # The atom order is the same, simply copy across properties as is.
        if not is_reordered:
            if verbose:
                print("\nSetting properties...")
            # Loop over all of the keys in the new molecule.
            for prop in props1:
                # Skip 'parameters' property, since it contains references to other parameters.
                if prop != param:
                    # This is a new property, or we are allowed to overwrite.
                    if (not mol0.hasProperty(_map[prop])) or overwrite:
                        if verbose:
                            print("  %s" % _map[prop])
                        try:
                            edit_mol = edit_mol.setProperty(_map[prop], mol1.property(prop))
                        except:
                            raise _IncompatibleError("Failed to set property '%s'" % _map[prop]) from None

        # The atom order is different, we need to map the atoms when setting properties.
        else:
            # Create a dictionary to flag whether a property has been seen.
            seen_prop = {}
            for prop in props1:
                seen_prop[prop] = False

            # First, set atom based properties.

            if verbose:
                print("\nSetting atom properties...")

            # Loop over all of the keys in the new molecule.
            for prop in props1:
                # This is a new property, or we are allowed to overwrite.
                if (not mol0.hasProperty(_map[prop])) or overwrite:
                    # Loop over all of the atom mapping pairs and set the property.
                    for idx0, idx1 in matches.items():
                        # Does the atom have this property?
                        # If so, add it to the matching atom in this molecule.
                        if mol1.atom(idx1).hasProperty(prop):
                            if verbose:
                                print("  %-20s %s --> %s" % (_map[prop], idx1, idx0))
                            try:
                                edit_mol = edit_mol.atom(idx0).setProperty(_map[prop], mol1.atom(idx1).property(prop)).molecule()
                                seen_prop[prop] = True
                            except:
                                raise _IncompatibleError("Failed to copy property '%s' from %s to %s."
                                    % (_map[prop], idx1, idx0)) from None

            # Now deal with all unseen properties. These will be non atom-based
            # properties, such as TwoAtomFunctions, StringProperty, etc.

            if verbose:
                print("\nSetting molecule properties...")

            # Loop over all of the unseen properties.
            for prop in seen_prop:
                if not seen_prop[prop]:
                    # Skip 'parameters' property, since it contains references to other parameters.
                    if prop != "parameters":
                        # This is a new property, or we are allowed to overwrite.
                        if (not mol0.hasProperty(_map[prop])) or overwrite:
                            if verbose:
                                print("  %s" % _map[prop])
                            # Try to match with atoms in the original molecule.
                            try:
                                edit_mol.setProperty(_map[prop], mol1.property(prop).makeCompatibleWith(mol0, matches))
                            # This is probably just a molecule property, such as a forcefield definition.
                            except AttributeError:
                                try:
                                    edit_mol.setProperty(_map[prop], mol1.property(prop))
                                except:
                                    raise _IncompatibleError("Failed to set property '%s'" % _map[prop]) from None

        # Finally, rename the atoms.

        if rename_atoms and is_renamed:
            if verbose:
                print("\nRenaming atoms...")

            for idx0, idx1 in matches.items():
                # Get the name of the atom in each molecule.
                name0 = mol0.atom(idx0).name()
                name1 = mol1.atom(idx1).name()

                if verbose:
                    print("  %s --> %s" % (name0, name1))

                # Try to rename the atom.
                try:
                    edit_mol = edit_mol.atom(idx0).rename(mol1.atom(idx1).name()).molecule()
                except:
                    raise _IncompatibleError("Failed to rename atom: %s --> %s" % (name0, name1)) from None

        # Commit the changes.
        self._sire_molecule = edit_mol.commit()

    def _convertFromMergedMolecule(self):
        """Convert from a merged molecule."""

        # We need to create two new Sire molecules and copy across the
        # atoms and properties from the two sub-molecules.

        # Try to get the list of atom indices in each molecule.
        try:
            idx0 = self._sire_molecule.property("atoms0")
            idx1 = self._sire_molecule.property("atoms1")
        except:
            raise _IncompatibleError("The merged molecule doesn't have the required properties!")

        # Create two selection objects to extract the atoms.
        sel0 = self._sire_molecule.selection()
        sel1 = self._sire_molecule.selection()
        sel0.deselectAll()
        sel1.deselectAll()

        try:
            # Loop over all atoms in the merged molecule and add the unique
            # atoms for each component.
            for idx in idx0:
                sel0.select(_SireMol.AtomIdx(idx))
            for idx in idx1:
                sel1.select(_SireMol.AtomIdx(idx))
        except:
            raise _IncompatibleError("Mismatched atom indices in the merged molecule!")

        # Use the selections to create two partial molecules for the atoms
        # at lamba = 0/1.
        par_mol0 = _SireMol.PartialMolecule(self._sire_molecule, sel0)
        par_mol1 = _SireMol.PartialMolecule(self._sire_molecule, sel1)

        # Extract the sub-molecules and make them editable.
        mol0 = par_mol0.extract().molecule().edit()
        mol1 = par_mol1.extract().molecule().edit()

        # Now rename all of the properties for each molecule and remove the
        # redundant values.

        try:
            # mol0
            for prop in mol0.propertyKeys():
                if prop[-1] == "0":
                    mol0.setProperty(prop[:-1], mol0.property(prop))
                    mol0.removeProperty(prop)
                elif prop[-1] == "1":
                    mol0.removeProperty(prop)
                elif prop == "is_perturbable":
                    mol0.removeProperty(prop)

            # mol1
            for prop in mol1.propertyKeys():
                if prop[-1] == "1":
                    mol1.setProperty(prop[:-1], mol1.property(prop))
                    mol1.removeProperty(prop)
                elif prop[-1] == "0":
                    mol1.removeProperty(prop)
                elif prop == "is_perturbable":
                    mol1.removeProperty(prop)
        except:
            raise _IncompatibleError("The merged molecule doesn't have the required properties!")

        # Finally, commit and save the molecules.
        self._molecule0 = Molecule(mol0.commit())
        self._molecule1 = Molecule(mol1.commit())

    def _fixCharge(self, map={}):
        """Make the molecular charge an integer value.

           Keyword arguments
           -----------------

           map : dict
               A dictionary that maps "properties" in this molecule to their
               user defined values. This allows the user to refer to properties
               with their own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Get the user defined charge property.
        if "charge" in map:
            prop = map["charge"]
        else:
            prop = "charge"

        # Calculate the charge.
        charge = self.charge().magnitude()

        # Calculate the difference from the nearest integer value.
        delta = round(charge) - charge

        # The difference is too small to care about.
        if charge + delta == _approx(charge):
            return

        # Normalise by the number of atoms.
        delta /= self.nAtoms()

        # Make the molecule editable.
        edit_mol = self._sire_molecule.edit()

        # Shift the charge of each atom in the molecule by delta.
        for atom in edit_mol.atoms():
            charge = edit_mol.atom(atom.index()).property(prop).value()
            edit_mol = edit_mol.atom(atom.index()) \
                               .setProperty(prop, (charge + delta)*_SireUnits.e_charge) \
                               .molecule()

        # Update the Sire molecule.
        self._sire_molecule = edit_mol.commit()

    def _merge(self, other, mapping, map0={}, map1={}):
        """Merge this molecule with 'other'.

           Positional arguments
           --------------------

           other : BioSimSpace._SireWrappers.Molecule
               The molecule to merge with.

           mapping : dict
               The mapping between matching atom indices in the two molecules.


           Keyword arguments
           -----------------

           map0 : dict
               A dictionary that maps "properties" in this molecule to their
               user defined values. This allows the user to refer to properties
               with their own naming scheme, e.g. { "charge" : "my-charge" }

           map1 : dict
               A dictionary that maps "properties" in ther other molecule to
               their user defined values.


            Returns
            -------

            merged : Sire.Mol.Molecule
                The merged molecule.
        """

        # Cannot merge an already merged molecule.
        if self._is_merged:
            raise IncompatibleError("This molecule has already been merged!")
        if other._is_merged:
            raise IncompatibleError("'other' has already been merged!")

        # Validate input.

        if type(other) is not Molecule:
            raise TypeError("'other' must be of type 'BioSimSpace._SireWrappers.Molecule'")

        if type(map0) is not dict:
            raise TypeError("'map0' must be of type 'dict'")

        if type(map1) is not dict:
            raise TypeError("'map1' must be of type 'dict'")

        if type(mapping) is not dict:
            raise TypeError("'mapping' must be of type 'dict'.")
        else:
            # Make sure all key/value pairs are of type AtomIdx.
            for idx0, idx1 in mapping.items():
                if type(idx0) is not _SireMol.AtomIdx or type(idx1) is not _SireMol.AtomIdx:
                    raise TypeError("key:value pairs in 'mapping' must be of type 'Sire.Mol.AtomIdx'")

        # Create a copy of this molecule.
        mol = Molecule(self)

        # Set the two molecule objects.
        molecule0 = mol._sire_molecule
        molecule1 = other._sire_molecule

        # Get the atom indices from the mapping.
        idx0 = mapping.keys()
        idx1 = mapping.values()

        # Create the reverse mapping: molecule1 --> molecule0
        inv_mapping = {v: k for k, v in mapping.items()}

        # Invert the user property mappings.
        inv_map0 = {v: k for k, v in map0.items()}
        inv_map1 = {v: k for k, v in map1.items()}

        # Make sure that the molecules have a "forcefield" property and that
        # the two force fields are compatible.

        # Get the user name for the "forcefield" property.
        ff0 = "forcefield"
        ff1 = "forcefield"
        if "forcefield" in inv_map0:
            ff0 = inv_map0["forcefield"]
        if "forcefield" in inv_map1:
            ff1 = inv_map1["forcefield"]

        # Force field information is missing.
        if not molecule0.hasProperty(ff0):
            raise _IncompatibleError("Cannot determine 'forcefield' of 'molecule0'!")
        if not molecule1.hasProperty(ff0):
            raise _IncompatibleError("Cannot determine 'forcefield' of 'molecule1'!")

        # The force fields are incompatible.
        if not molecule0.property(ff0).isCompatibleWith(molecule1.property(ff1)):
            raise _IncompatibleError("Cannot merge molecules with incompatible force fields!")

        # Create lists to store the atoms that are unique to each molecule,
        # along with their indices.
        atoms0 = []
        atoms1 = []
        atoms0_idx = []
        atoms1_idx = []

        # Loop over each molecule to find the unique atom indices.

        # molecule0
        for atom in molecule0.atoms():
            if atom.index() not in idx0:
                atoms0.append(atom)
                atoms0_idx.append(atom.index())

        # molecule1
        for atom in molecule1.atoms():
            if atom.index() not in idx1:
                atoms1.append(atom)
                atoms1_idx.append(atom.index())

        # Create lists of the actual property names in the molecules.
        props0 = []
        props1 = []

        # molecule0
        for prop in molecule0.propertyKeys():
            if prop in inv_map0:
                prop = inv_map0[prop]
            props0.append(prop)

        # molecule1
        for prop in molecule1.propertyKeys():
            if prop in inv_map1:
                prop = inv_map1[prop]
            props1.append(prop)

        # Determine the common properties between the two molecules.
        # These are the properties that can be perturbed.
        shared_props = list(set(props0).intersection(props1))
        del(props0)
        del(props1)

        # Create a new molecule to hold the merged molecule.
        molecule = _SireMol.Molecule("Merged_Molecule")

        # Add a single residue called LIG.
        res = molecule.edit().add(_SireMol.ResNum(1))
        res.rename(_SireMol.ResName("LIG"))

        # Create a single cut-group.
        cg = res.molecule().add(_SireMol.CGName("1"))

        # Counter for the number of atoms.
        num = 1

        # First add all of the atoms from molecule0.
        for atom in molecule0.atoms():
            # Add the atom.
            added = cg.add(atom.name())
            added.renumber(_SireMol.AtomNum(num))
            added.reparent(_SireMol.ResIdx(0))
            num += 1

        # Now add all of the atoms from molecule1 that aren't mapped from molecule0.
        for atom in atoms1:
            added = cg.add(atom.name())
            added.renumber(_SireMol.AtomNum(num))
            added.reparent(_SireMol.ResIdx(0))
            inv_mapping[atom.index()] = _SireMol.AtomIdx(num-1)
            num += 1

        # Commit the changes to the molecule.
        molecule = cg.molecule().commit()

        # Make the molecule editable.
        edit_mol = molecule.edit()

        # We now add properties to the merged molecule. The properties are used
        # to represent the molecule at two states along the alchemical pathway:
        #
        # lambda = 0:
        #   Here only molecule0 is active. The "charge" and "LJ" properties
        #   for atoms that are part of molecule1 are set to zero. Other force
        #   field properties, e.g. "bond", "angle", "dihedral", and "improper",
        #   are retained for the atoms in molecule1, although the indices of
        #   the atoms involved in the interactions must be re-mapped to their
        #   positions in the merged molecule. (Also, the interactions may now
        #   be between atoms of different type.) The properties are given the
        #   suffix "0", e.g. "charge0".
        #
        # lambda = 1:
        #   Here only molecule1 is active. We perform the same process as above,
        #   only we modify the properties of the atoms that are unique to
        #   molecule0. The properties are given the suffix "1", e.g. "charge1".
        #
        # Properties that aren't shared between the molecules (and thus can't
        # be merged) are set using their original names.

        ##############################
        # SET PROPERTIES AT LAMBDA = 0
        ##############################

        # Add the atom properties from molecule0.
        for atom in molecule0.atoms():
            # Loop over all atom properties.
            for prop in atom.propertyKeys():
                # Get the actual property name.
                if prop in inv_map0:
                    name = inv_map0[prop]
                else:
                    name = prop

                # This is a perturbable property. Rename to "property0", e.g. "charge0".
                if name in shared_props:
                    name = name + "0"

                # Add the property to the atom in the merged molecule.
                edit_mol = edit_mol.atom(atom.index()).setProperty(name, atom.property(prop)).molecule()

        # Add the atom properties from molecule1.
        for atom in atoms1:
            # Get the atom index in the merged molecule.
            idx = inv_mapping[atom.index()]

            # Loop over all atom properties.
            for prop in atom.propertyKeys():
                # Get the actual property name.
                if prop in inv_map1:
                    name = inv_map1[prop]
                else:
                    name = prop

                # Zero the "charge" and "LJ" property for atoms that are unique to molecule1.
                if name == "charge":
                    edit_mol = edit_mol.atom(idx).setProperty("charge0", 0*_SireUnits.e_charge).molecule()
                elif name == "LJ":
                    edit_mol = edit_mol.atom(idx).setProperty("LJ0", _SireMM.LJParameter()).molecule()
                else:
                    # This is a perturbable property. Rename to "property0", e.g. "charge0".
                    if name in shared_props:
                        name = name + "0"

                    # Add the property to the atom in the merged molecule.
                    edit_mol = edit_mol.atom(idx).setProperty(name, atom.property(prop)).molecule()

        # We now need to merge "bond", "angle", "dihedral", and "improper" parameters.
        # To do so, we extract the properties from molecule0, then add the additional
        # properties from molecule1, making sure to update the atom indices, and bond
        # atoms from molecule1 to the atoms to which they map in molecule0.

        # 1) bonds
        if "bond" in shared_props:
            # Get the user defined property names.
            prop0 = "bond"
            if prop0 in inv_map0:
                prop0 = inv_map0["bond"]

            prop1 = "bond"
            if prop1 in inv_map1:
                prop1 = inv_map1["bond"]

            # Get the "bond" property from the two molecules.
            bonds0 = molecule0.property(prop0)
            bonds1 = molecule1.property(prop1)

            # Get the molInfo object for molecule1.
            info = molecule1.info()

            # Create the new set of bonds.
            bonds = _SireMM.TwoAtomFunctions(edit_mol.info())

            # Add all of the bonds from molecule0.
            for bond in bonds0.potentials():
                bonds.set(bond.atom0(), bond.atom1(), bond.function())

            # Loop over all bonds in molecule1.
            for bond in bonds1.potentials():
                # This bond contains an atom that is unique to molecule1.
                if info.atomIdx(bond.atom0()) in atoms1_idx or \
                   info.atomIdx(bond.atom1()) in atoms1_idx:

                    # Extract the bond information.
                    atom0 = info.atomIdx(bond.atom0())
                    atom1 = info.atomIdx(bond.atom1())
                    exprn = bond.function()

                    # Map the atom indices to their position in the merged molecule.
                    atom0 = inv_mapping[atom0]
                    atom1 = inv_mapping[atom1]

                    # Set the new bond.
                    bonds.set(atom0, atom1, exprn)

            # Add the bonds to the merged molecule.
            edit_mol.setProperty("bond0", bonds)

        # 2) angles
        if "angle" in shared_props:
            # Get the user defined property names.
            prop0 = "angle"
            if prop0 in inv_map0:
                prop0 = inv_map0["angle"]

            prop1 = "angle"
            if prop1 in inv_map1:
                prop1 = inv_map1["angle"]

            # Get the "angle" property from the two molecules.
            angles0 = molecule0.property(prop0)
            angles1 = molecule1.property(prop1)

            # Get the molInfo object for molecule1.
            info = molecule1.info()

            # Create the new set of angles.
            angles = _SireMM.ThreeAtomFunctions(edit_mol.info())

            # Add all of the angles from molecule0.
            for angle in angles0.potentials():
                angles.set(angle.atom0(), angle.atom1(), angle.atom2(), angle.function())

            # Loop over all angles in molecule1.
            for angle in angles1.potentials():
                # This angle contains an atom that is unique to molecule1.
                if info.atomIdx(angle.atom0()) in atoms1_idx or \
                   info.atomIdx(angle.atom1()) in atoms1_idx or \
                   info.atomIdx(angle.atom2()) in atoms1_idx:

                       # Extract the angle information.
                       atom0 = info.atomIdx(angle.atom0())
                       atom1 = info.atomIdx(angle.atom1())
                       atom2 = info.atomIdx(angle.atom2())
                       exprn = angle.function()

                       # Map the atom indices to their position in the merged molecule.
                       atom0 = inv_mapping[atom0]
                       atom1 = inv_mapping[atom1]
                       atom2 = inv_mapping[atom2]

                       # Set the new angle.
                       angles.set(atom0, atom1, atom2, exprn)

            # Add the angles to the merged molecule.
            edit_mol.setProperty("angle0", angles)

        # 3) dihedrals
        if "dihedral" in shared_props:
            # Get the user defined property names.
            prop0 = "dihedral"
            if prop0 in inv_map0:
                prop0 = inv_map0["dihedral"]

            prop1 = "dihedral"
            if prop1 in inv_map1:
                prop1 = inv_map1["dihedral"]

            # Get the "dihedral" property from the two molecules.
            dihedrals0 = molecule0.property(prop0)
            dihedrals1 = molecule1.property(prop1)

            # Get the molInfo object for molecule1.
            info = molecule1.info()

            # Create the new set of dihedrals.
            dihedrals = _SireMM.FourAtomFunctions(edit_mol.info())

            # Add all of the dihedrals from molecule0.
            for dihedral in dihedrals0.potentials():
                dihedrals.set(dihedral.atom0(), dihedral.atom1(),
                    dihedral.atom2(), dihedral.atom3(), dihedral.function())

            # Loop over all dihedrals in molecule1.
            for dihedral in dihedrals1.potentials():
                # This dihedral contains an atom that is unique to molecule1.
                if info.atomIdx(dihedral.atom0()) in atoms1_idx or \
                   info.atomIdx(dihedral.atom1()) in atoms1_idx or \
                   info.atomIdx(dihedral.atom2()) in atoms1_idx or \
                   info.atomIdx(dihedral.atom3()) in atoms1_idx:

                       # Extract the dihedral information.
                       atom0 = info.atomIdx(dihedral.atom0())
                       atom1 = info.atomIdx(dihedral.atom1())
                       atom2 = info.atomIdx(dihedral.atom2())
                       atom3 = info.atomIdx(dihedral.atom3())
                       exprn = dihedral.function()

                       # Map the atom indices to their position in the merged molecule.
                       atom0 = inv_mapping[atom0]
                       atom1 = inv_mapping[atom1]
                       atom2 = inv_mapping[atom2]
                       atom3 = inv_mapping[atom3]

                       # Set the new dihedral.
                       dihedrals.set(atom0, atom1, atom2, atom3, exprn)

            # Add the dihedrals to the merged molecule.
            edit_mol.setProperty("dihedral0", dihedrals)

        # 4) impropers
        if "improper" in shared_props:
            # Get the user defined property names.
            prop0 = "improper"
            if prop0 in inv_map0:
                prop0 = inv_map0["improper"]

            prop1 = "improper"
            if prop1 in inv_map1:
                prop1 = inv_map1["improper"]

            # Get the "improper" property from the two molecules.
            impropers0 = molecule0.property(prop0)
            impropers1 = molecule1.property(prop1)

            # Get the molInfo object for molecule1.
            info = molecule1.info()

            # Create the new set of impropers.
            impropers = _SireMM.FourAtomFunctions(edit_mol.info())

            # Add all of the impropers from molecule0.
            for improper in impropers0.potentials():
                impropers.set(improper.atom0(), improper.atom1(),
                    improper.atom2(), improper.atom3(), improper.function())

            # Loop over all impropers in molecule1.
            for improper in impropers1.potentials():
                # This improper contains an atom that is unique to molecule1.
                if info.atomIdx(improper.atom0()) in atoms1_idx or \
                   info.atomIdx(improper.atom1()) in atoms1_idx or \
                   info.atomIdx(improper.atom2()) in atoms1_idx or \
                   info.atomIdx(improper.atom3()) in atoms1_idx:

                       # Extract the improper information.
                       atom0 = info.atomIdx(improper.atom0())
                       atom1 = info.atomIdx(improper.atom1())
                       atom2 = info.atomIdx(improper.atom2())
                       atom3 = info.atomIdx(improper.atom3())
                       exprn = improper.function()

                       # Map the atom indices to their position in the merged molecule.
                       atom0 = inv_mapping[atom0]
                       atom1 = inv_mapping[atom1]
                       atom2 = inv_mapping[atom2]
                       atom3 = inv_mapping[atom3]

                       # Set the new improper.
                       impropers.set(atom0, atom1, atom2, atom3, exprn)

            # Add the impropers to the merged molecule.
            edit_mol.setProperty("improper0", impropers)

        ##############################
        # SET PROPERTIES AT LAMBDA = 1
        ##############################

        # Add the atom properties from molecule1.
        for atom in molecule1.atoms():
            # Get the atom index in the merged molecule.
            idx = inv_mapping[atom.index()]

            # Loop over all atom properties.
            for prop in atom.propertyKeys():
                # Get the actual property name.
                if prop in inv_map1:
                    name = inv_map1[prop]
                else:
                    name = prop

                # This is a perturbable property. Rename to "property1", e.g. "charge1".
                if name in shared_props:
                    name = name + "1"

                # Add the property to the atom in the merged molecule.
                edit_mol = edit_mol.atom(idx).setProperty(name, atom.property(prop)).molecule()

        # Add the properties from atoms unique to molecule0.
        for atom in atoms0:
            # Loop over all atom properties.
            for prop in atom.propertyKeys():
                # Get the actual property name.
                if prop in inv_map0:
                    name = inv_map0[prop]
                else:
                    name = prop

                # Zero the "charge" and "LJ" property for atoms that are unique to molecule0.
                if name == "charge":
                    edit_mol = edit_mol.atom(atom.index()).setProperty("charge1", 0*_SireUnits.e_charge).molecule()
                elif name == "LJ":
                    edit_mol = edit_mol.atom(atom.index()).setProperty("LJ1", _SireMM.LJParameter()).molecule()
                else:
                    # This is a perturbable property. Rename to "property1", e.g. "charge1".
                    if name in shared_props:
                        name = name + "1"

                    # Add the property to the atom in the merged molecule.
                    edit_mol = edit_mol.atom(atom.index()).setProperty(name, atom.property(prop)).molecule()

        # We now need to merge "bond", "angle", "dihedral", and "improper" parameters.
        # To do so, we extract the properties from molecule1, then add the additional
        # properties from molecule0, making sure to update the atom indices, and bond
        # atoms from molecule0 to the atoms to which they map in molecule1.

        # 1) bonds
        if "bond" in shared_props:
            # Get the info objects for the two molecules.
            info0 = molecule0.info()
            info1 = molecule1.info()

            # Get the user defined property names.
            prop0 = "bond"
            if prop0 in inv_map0:
                prop0 = inv_map0["bond"]

            prop1 = "bond"
            if prop1 in inv_map1:
                prop1 = inv_map1["bond"]

            # Get the "bond" property from the two molecules.
            bonds0 = molecule0.property(prop0)
            bonds1 = molecule1.property(prop1)

            # Create the new set of bonds.
            bonds = _SireMM.TwoAtomFunctions(edit_mol.info())

            # Add all of the bonds from molecule1.
            for bond in bonds1.potentials():
                # Extract the bond information.
                atom0 = info1.atomIdx(bond.atom0())
                atom1 = info1.atomIdx(bond.atom1())
                exprn = bond.function()

                # Map the atom indices to their position in the merged molecule.
                atom0 = inv_mapping[atom0]
                atom1 = inv_mapping[atom1]

                # Set the new bond.
                bonds.set(atom0, atom1, exprn)

            # Loop over all bonds in molecule0
            for bond in bonds0.potentials():
                # This bond contains an atom that is unique to molecule0.
                if info0.atomIdx(bond.atom0()) in atoms0_idx or \
                   info0.atomIdx(bond.atom1()) in atoms0_idx:

                   # Extract the bond information.
                   atom0 = info0.atomIdx(bond.atom0())
                   atom1 = info0.atomIdx(bond.atom1())
                   exprn = bond.function()

                   # Set the new bond.
                   bonds.set(atom0, atom1, exprn)

            # Add the bonds to the merged molecule.
            edit_mol.setProperty("bond1", bonds)

        # 2) angles
        if "angle" in shared_props:
            # Get the info objects for the two molecules.
            info0 = molecule0.info()
            info1 = molecule1.info()

            # Get the user defined property names.
            prop0 = "angle"
            if prop0 in inv_map0:
                prop0 = inv_map0["angle"]

            prop1 = "angle"
            if prop1 in inv_map1:
                prop1 = inv_map1["angle"]

            # Get the "angle" property from the two molecules.
            angles0 = molecule0.property(prop0)
            angles1 = molecule1.property(prop1)

            # Create the new set of angles.
            angles = _SireMM.ThreeAtomFunctions(edit_mol.info())

            # Add all of the angles from molecule1.
            for angle in angles1.potentials():
                # Extract the angle information.
                atom0 = info1.atomIdx(angle.atom0())
                atom1 = info1.atomIdx(angle.atom1())
                atom2 = info1.atomIdx(angle.atom2())
                exprn = angle.function()

                # Map the atom indices to their position in the merged molecule.
                atom0 = inv_mapping[atom0]
                atom1 = inv_mapping[atom1]
                atom2 = inv_mapping[atom2]

                # Set the new angle.
                angles.set(atom0, atom1, atom2, exprn)

            # Loop over all angles in molecule0.
            for angle in angles0.potentials():
                # This angle contains an atom that is unique to molecule0.
                if info0.atomIdx(angle.atom0()) in atoms0_idx or \
                   info0.atomIdx(angle.atom1()) in atoms0_idx or \
                   info0.atomIdx(angle.atom2()) in atoms0_idx:

                       # Extract the angle information.
                       atom0 = info0.atomIdx(angle.atom0())
                       atom1 = info0.atomIdx(angle.atom1())
                       atom2 = info0.atomIdx(angle.atom2())
                       exprn = angle.function()

                       # Set the new angle.
                       angles.set(atom0, atom1, atom2, exprn)

            # Add the angles to the merged molecule.
            edit_mol.setProperty("angle1", angles)

        # 3) dihedrals
        if "dihedral" in shared_props:
            # Get the info objects for the two molecules.
            info0 = molecule0.info()
            info1 = molecule1.info()

            # Get the user defined property names.
            prop0 = "dihedral"
            if prop0 in inv_map0:
                prop0 = inv_map0["dihedral"]

            prop1 = "dihedral"
            if prop1 in inv_map1:
                prop1 = inv_map1["dihedral"]

            # Get the "dihedral" property from the two molecules.
            dihedrals0 = molecule0.property(prop0)
            dihedrals1 = molecule1.property(prop1)

            # Create the new set of dihedrals.
            dihedrals = _SireMM.FourAtomFunctions(edit_mol.info())

            # Add all of the dihedrals from molecule1.
            for dihedral in dihedrals1.potentials():
                # Extract the dihedral information.
                atom0 = info1.atomIdx(dihedral.atom0())
                atom1 = info1.atomIdx(dihedral.atom1())
                atom2 = info1.atomIdx(dihedral.atom2())
                atom3 = info1.atomIdx(dihedral.atom3())
                exprn = dihedral.function()

                # Map the atom indices to their position in the merged molecule.
                atom0 = inv_mapping[atom0]
                atom1 = inv_mapping[atom1]
                atom2 = inv_mapping[atom2]
                atom3 = inv_mapping[atom3]

                # Set the new dihedral.
                dihedrals.set(atom0, atom1, atom2, atom3, exprn)

            # Loop over all dihedrals in molecule0.
            for dihedral in dihedrals0.potentials():
                # This dihedral contains an atom that is unique to molecule0.
                if info0.atomIdx(dihedral.atom0()) in atoms0_idx or \
                   info0.atomIdx(dihedral.atom1()) in atoms0_idx or \
                   info0.atomIdx(dihedral.atom2()) in atoms0_idx or \
                   info0.atomIdx(dihedral.atom3()) in atoms0_idx:

                       # Extract the dihedral information.
                       atom0 = info0.atomIdx(dihedral.atom0())
                       atom1 = info0.atomIdx(dihedral.atom1())
                       atom2 = info0.atomIdx(dihedral.atom2())
                       atom3 = info0.atomIdx(dihedral.atom3())
                       exprn = dihedral.function()

                       # Set the new dihedral.
                       dihedrals.set(atom0, atom1, atom2, atom3, exprn)

            # Add the dihedrals to the merged molecule.
            edit_mol.setProperty("dihedral1", dihedrals)

        # 4) impropers
        if "improper" in shared_props:
            # Get the info objects for the two molecules.
            info0 = molecule0.info()
            info1 = molecule1.info()

            # Get the user defined property names.
            prop0 = "improper"
            if prop0 in inv_map0:
                prop0 = inv_map0["improper"]

            prop1 = "improper"
            if prop1 in inv_map1:
                prop1 = inv_map1["improper"]

            # Get the "improper" property from the two molecules.
            impropers0 = molecule0.property(prop0)
            impropers1 = molecule1.property(prop1)

            # Create the new set of impropers.
            impropers = _SireMM.FourAtomFunctions(edit_mol.info())

            # Add all of the impropers from molecule1.
            for improper in impropers1.potentials():
                # Extract the improper information.
                atom0 = info1.atomIdx(improper.atom0())
                atom1 = info1.atomIdx(improper.atom1())
                atom2 = info1.atomIdx(improper.atom2())
                atom3 = info1.atomIdx(improper.atom3())
                exprn = improper.function()

                # Map the atom indices to their position in the merged molecule.
                atom0 = inv_mapping[atom0]
                atom1 = inv_mapping[atom1]
                atom2 = inv_mapping[atom2]
                atom3 = inv_mapping[atom3]

                # Set the new improper.
                impropers.set(atom0, atom1, atom2, atom3, exprn)

            # Loop over all impropers in molecule0.
            for improper in impropers0.potentials():
                # This improper contains an atom that is unique to molecule0.
                if info0.atomIdx(improper.atom0()) in atoms0_idx or \
                   info0.atomIdx(improper.atom1()) in atoms0_idx or \
                   info0.atomIdx(improper.atom2()) in atoms0_idx or \
                   info0.atomIdx(improper.atom3()) in atoms0_idx:

                       # Extract the improper information.
                       atom0 = info0.atomIdx(improper.atom0())
                       atom1 = info0.atomIdx(improper.atom1())
                       atom2 = info0.atomIdx(improper.atom2())
                       atom3 = info0.atomIdx(improper.atom3())
                       exprn = improper.function()

                       # Set the new improper.
                       impropers.set(atom0, atom1, atom2, atom3, exprn)

            # Add the impropers to the merged molecule.
            edit_mol.setProperty("improper1", impropers)

        # The number of potentials should be consistent for the "bond0"
        # and "bond1" properties.
        if edit_mol.property("bond0").nFunctions() != edit_mol.property("bond1").nFunctions():
            raise RuntimeError("Inconsistent number of bonds in merged molecule!")

        # Create the connectivity object
        conn = _SireMol.Connectivity(edit_mol.info()).edit()

        # Connect the bonded atoms. Connectivity is the same at lambda = 0
        # and lambda = 1.
        for bond in edit_mol.property("bond0").potentials():
            conn.connect(bond.atom0(), bond.atom1())
        conn = conn.commit()

        # Get the connectivity of the two molecules.
        c0 = molecule0.property("connectivity")
        c1 = molecule1.property("connectivity")

        # Check that the merge hasn't modified the connectivity.

        # molecule0
        for x in range(0, molecule0.nAtoms()):
            for y in range(x+1, molecule0.nAtoms()):
                if c0.connectionType(_SireMol.AtomIdx(x), _SireMol.AtomIdx(y)) != \
                   conn.connectionType(_SireMol.AtomIdx(x), _SireMol.AtomIdx(y)):
                       raise _IncompatibleError("Merge has changed the molecular connectivity! "
                           + "Check your atom mapping.")

        # molecule1
        for x in range(0, molecule1.nAtoms()):
            # Convert to an AtomIdx.
            idx = _SireMol.AtomIdx(x)

            # Map the index to its position in the merged molecule.
            idx = inv_mapping[idx]

            for y in range(x+1, molecule1.nAtoms()):
                # Convert to an AtomIdx.
                idy = _SireMol.AtomIdx(y)

                # Map the index to its position in the merged molecule.
                idy = inv_mapping[idy]

                if c1.connectionType(_SireMol.AtomIdx(x), _SireMol.AtomIdx(y)) != \
                   conn.connectionType(idx, idy):
                       raise _IncompatibleError("Merge has changed the molecular connectivity! "
                           + "Check your atom mapping.")

        # Set the "connectivity" property.
        edit_mol.setProperty("connectivity", conn)

        # Create the CLJNBPairs matrices.
        ff = molecule0.property(ff0)

        clj_nb_pairs0 = _SireMM.CLJNBPairs(edit_mol.info(),
            _SireMM.CLJScaleFactor(0, 0))

        # Loop over all atoms unique to molecule0.
        for idx0 in atoms0_idx:
            # Loop over all atoms unique to molecule1.
            for idx1 in atoms1_idx:
                # Map the index to its position in the merged molecule.
                idx1 = inv_mapping[idx1]

                # Work out the connection type between the atoms.
                conn_type = conn.connectionType(idx0, idx1)

                # The atoms aren't bonded.
                if conn_type == 0:
                    clj_scale_factor = _SireMM.CLJScaleFactor(1, 1)
                    clj_nb_pairs0.set(idx0, idx1, clj_scale_factor)

                # The atoms are part of a dihedral.
                elif conn_type == 4:
                    clj_scale_factor = _SireMM.CLJScaleFactor(ff.electrostatic14ScaleFactor(),
                                                              ff.vdw14ScaleFactor())
                    clj_nb_pairs0.set(idx0, idx1, clj_scale_factor)

        # Copy the intrascale matrix.
        clj_nb_pairs1 = clj_nb_pairs0.__deepcopy__()

        # Get the user defined "intrascale" property names.
        prop0 = "intrascale"
        if prop0 in inv_map0:
            prop0 = inv_map0["intrascale"]

        prop1 = "intrascale"
        if prop1 in inv_map1:
            prop1 = inv_map1["intrascale"]

        # Get the "intrascale" property from the two molecules.
        intrascale0 = molecule0.property(prop0)
        intrascale1 = molecule1.property(prop1)

        # Copy the intrascale from molecule1 into clj_nb_pairs0.

        # Perform a triangular loop over atoms from molecule1.
        for x in range(0, molecule1.nAtoms()):
            # Convert to an AtomIdx.
            idx = _SireMol.AtomIdx(x)

            # Map the index to its position in the merged molecule.
            idx = inv_mapping[idx]

            for y in range(x+1, molecule1.nAtoms()):
                # Convert to an AtomIdx.
                idy = _SireMol.AtomIdx(y)

                # Map the index to its position in the merged molecule.
                idy = inv_mapping[idy]

                # Get the intrascale value.
                intra = intrascale1.get(_SireMol.AtomIdx(x), _SireMol.AtomIdx(y))

                # Only set if there is a non-zero value.
                # Set using the re-mapped atom indices.
                if not intra.coulomb() == 0:
                    clj_nb_pairs0.set(idx, idy, intra)

        # Now copy in all intrascale values from molecule0 into both
        # clj_nb_pairs matrices.

        # Perform a triangular loop over atoms from molecule0.
        for x in range(0, molecule0.nAtoms()):
            for y in range(x+1, molecule0.nAtoms()):
                # Get the intrascale value.
                intra = intrascale0.get(_SireMol.AtomIdx(x), _SireMol.AtomIdx(y))

                # Set the value in the new matrix, overwriting existing value.
                clj_nb_pairs0.set(_SireMol.AtomIdx(x), _SireMol.AtomIdx(y), intra)

                # Only set if there is a non-zero value.
                if not intra.coulomb() == 0:
                    clj_nb_pairs1.set(_SireMol.AtomIdx(x), _SireMol.AtomIdx(y), intra)

        # Finally, copy the intrascale from molecule1 into clj_nb_pairs1.

        # Perform a triangular loop over atoms from molecule1.
        for x in range(0, molecule1.nAtoms()):
            # Convert to an AtomIdx.
            idx = _SireMol.AtomIdx(x)

            # Map the index to its position in the merged molecule.
            idx = inv_mapping[idx]

            for y in range(x+1, molecule1.nAtoms()):
                # Convert to an AtomIdx.
                idy = _SireMol.AtomIdx(y)

                # Map the index to its position in the merged molecule.
                idy = inv_mapping[idy]

                # Get the intrascale value.
                intra = intrascale1.get(_SireMol.AtomIdx(x), _SireMol.AtomIdx(y))

                # Set the value in the new matrix, overwriting existing value.
                clj_nb_pairs1.set(idx, idy, intra)

        # Set the "intrascale" properties.
        edit_mol.setProperty("intrascale0", clj_nb_pairs0)
        edit_mol.setProperty("intrascale1", clj_nb_pairs1)

        # Set the "forcefield" properties.
        edit_mol.setProperty("forcefield0", molecule0.property(ff0))
        edit_mol.setProperty("forcefield1", molecule1.property(ff1))

        # Flag that this molecule is perturbable.
        edit_mol.setProperty("is_perturbable", _SireBase.wrap(True))

        # Store the indices of the atoms that are unique to each molecule.
        atom_list0 = _SireBase.IntegerArrayProperty()
        atom_list1 = _SireBase.IntegerArrayProperty()
        for atom in molecule0.atoms():
            atom_list0.append(atom.index().value())
        for atom in molecule1.atoms():
            atom_list1.append(inv_mapping[atom.index()].value())
        edit_mol.setProperty("atoms0", atom_list0)
        edit_mol.setProperty("atoms1", atom_list1)

        # Update the Sire molecule object of the new molecule.
        mol._sire_molecule = edit_mol.commit()

        # Flag that the molecule has been merged.
        mol._is_merged = True

        # Store the components of the merged molecule.
        mol._molecule0 = Molecule(molecule0)
        mol._molecule1 = Molecule(molecule1)

        # Return the new molecule.
        return mol

    def _getAABox(self, map={}):
        """Get the axis-aligned bounding box for the molecule.

           Keyword arguments
           -----------------

           map : dict
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

        # Extract the atomic coordinates and append them to the vector.
        try:
            if "coordinates" in map:
                prop = map["coordinates"]
            else:
                prop = "coordinates"
            coord.extend(self._sire_molecule.property(prop).toVector())

        except UserWarning:
            raise UserWarning("Molecule has no coordinate property.") from None

        # Return the AABox for the coordinates.
        return _SireVol.AABox(coord)

# Import at bottom of module to avoid circular dependency.
from ._system import System as _System
