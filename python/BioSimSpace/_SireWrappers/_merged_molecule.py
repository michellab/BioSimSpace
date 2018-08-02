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
single- and dual-topology free energy calculations.

Author: Lester Hedges <lester.hedges@gmail.com>
"""

import Sire.Maths as _SireMaths
import Sire.MM as _SireMM
import Sire.Mol as _SireMol
import Sire.Units as _SireUnits

__all__ = ["MergedMolecule"]

class MergedMolecule():
    """A container class for storing merged molecules."""

    def __init__(self, molecule0, molecule1, mapping, map0={}, map1={}):
        """Constructor.

           Positional arguments
           --------------------

           molecule0 : Sire.Mol.Molecule
               The initial molecule.

           molecule1 : Sire.Mol.Molecule
               The final molecule.

           mapping : dict
               The mapping between matching atom indices in the two molecules.


           Keyword arguments
           -----------------

           map0 : dict
               A dictionary that maps "properties" in molecule0 to their user
               defined values. This allows the user to refer to properties
               with their own naming scheme, e.g. { "charge" : "my-charge" }

           map1 : dict
               A dictionary that maps "properties" in molecule1 to their user
               defined values.
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

        if type(map0) is not dict:
            raise TypeError("'map0' must be of type 'dict'")

        if type(map1) is not dict:
            raise TypeError("'map1' must be of type 'dict'")

        # Create the merged molecule.
        self._sire_molecule = self._merge(self._sire_molecule0, self._sire_molecule1,
                                          mapping, map0, map1)

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

    def _getSireMolecule(self):
        """Return the full Sire Molecule object."""
        return self._sire_molecule

    def _merge(self, molecule0, molecule1, mapping, map0={}, map1={}):
        """Create the merged molecule.

           Positional arguments
           --------------------

           molecule0 : Sire.Mol.Molecule
               The initial molecule.

           molecule1 : Sire.Mol.Molecule
               The final molecule.

           mapping : dict
               The mapping between matching atom indices in the two molecules.


           Keyword arguments
           -----------------

           map0 : dict
               A dictionary that maps "properties" in molecule0 to their user
               defined values. This allows the user to refer to properties
               with their own naming scheme, e.g. { "charge" : "my-charge" }

           map1 : dict
               A dictionary that maps "properties" in molecule1 to their user
               defined values.


            Returns
            -------

            merged : Sire.Mol.Molecule
                The merged molecule.
        """

        # Get the atom indices from the mapping.
        idx0 = mapping.keys()
        idx1 = mapping.values()

        # Create the reverse mapping: molecule1 --> molecule0
        inv_mapping = {v: k for k, v in mapping.items()}

        # Invert the user property mappings.
        inv_map0 = {v: k for k, v in map0.items()}
        inv_map1 = {v: k for k, v in map1.items()}

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
        molecule = _SireMol.Molecule("Merged Molecule")

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

        # Create a dictionary to map between the indices of the unique atoms in
        # molecule1 to their index within the new merged molecule.
        new_idx = {}

        # Now add all of the atoms from molecule1 that aren't mapped from molecule0.
        for atom in atoms1:
            added = cg.add(atom.name())
            added.renumber(_SireMol.AtomNum(num))
            added.reparent(_SireMol.ResIdx(0))
            new_idx[atom.index()] = _SireMol.AtomIdx(num-1)
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
            idx = new_idx[atom.index()]

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
                    lj = _SireMM.LJParameter(0*_SireUnits.angstrom, 0*_SireUnits.kcal_per_mol)
                    edit_mol = edit_mol.atom(idx).setProperty("LJ0", lj).molecule()
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

            # Create a list to store the indices of bonds in molecule1 that
            # involve atoms that are part of the merged molecule at lambda = 0.
            bond_idxs = []

            # Get the molInfo object for molecule1.
            info = molecule1.info()

            # Loop over all bonds in molecule1.
            for idx, bond in enumerate(bonds1.potentials()):
                # This bond contains an atom that is unique to molecule1.
                if info.atomIdx(bond.atom0()) in atoms1_idx or \
                   info.atomIdx(bond.atom1()) in atoms1_idx:
                       bond_idxs.append(idx)

            # Store the bond potentials for molecule1.
            potentials = bonds1.potentials()

            # Create the new set of bonds.
            bonds = _SireMM.TwoAtomFunctions(edit_mol.info())

            # Add all of the bonds from molecule0.
            for bond in bonds0.potentials():
                bonds.set(bond.atom0(), bond.atom1(), bond.function())

            # Loop over all of the additional bonds from molecule1.
            for idx in bond_idxs:
                # Extract the bond information.
                atom0 = info.atomIdx(potentials[idx].atom0())
                atom1 = info.atomIdx(potentials[idx].atom1())
                exprn = potentials[idx].function()

                # Map the atom indices to their position in the merged molecule.

                if atom0 in inv_mapping:
                    atom0 = inv_mapping[atom0]
                else:
                    atom0 = new_idx[atom0]

                if atom1 in inv_mapping:
                    atom1 = inv_mapping[atom1]
                else:
                    atom1 = new_idx[atom1]

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

            # Create a list to store the indices of angles in molecule1 that
            # involve atoms that are part of the merged molecule at lambda = 0.
            angle_idxs = []

            # Get the molInfo object for molecule1.
            info = molecule1.info()

            # Loop over all angles in molecule1.
            for idx, angle in enumerate(angles1.potentials()):
                # This angle contains an atom that is unique to molecule1.
                if info.atomIdx(angle.atom0()) in atoms1_idx or \
                   info.atomIdx(angle.atom1()) in atoms1_idx or \
                   info.atomIdx(angle.atom2()) in atoms1_idx:
                       angle_idxs.append(idx)

            # Store the angle potentials for molecule1.
            potentials = angles1.potentials()

            # Create the new set of angles.
            angles = _SireMM.ThreeAtomFunctions(edit_mol.info())

            # Add all of the angles from molecule0.
            for angle in angles0.potentials():
                angles.set(angle.atom0(), angle.atom1(), angle.atom2(), angle.function())

            # Loop over all of the additional angles from molecule1.
            for idx in angle_idxs:
                # Extract the angle information.
                atom0 = info.atomIdx(potentials[idx].atom0())
                atom1 = info.atomIdx(potentials[idx].atom1())
                atom2 = info.atomIdx(potentials[idx].atom2())
                exprn = potentials[idx].function()

                # Map the atom indices to their position in the merged molecule.

                if atom0 in inv_mapping:
                    atom0 = inv_mapping[atom0]
                else:
                    atom0 = new_idx[atom0]

                if atom1 in inv_mapping:
                    atom1 = inv_mapping[atom1]
                else:
                    atom1 = new_idx[atom1]

                if atom2 in inv_mapping:
                    atom2 = inv_mapping[atom2]
                else:
                    atom2 = new_idx[atom2]

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

            # Create a list to store the indices of dihedrals in molecule1 that
            # involve atoms that are part of the merged molecule at lambda = 0.
            dihedral_idxs = []

            # Get the molInfo object for molecule1.
            info = molecule1.info()

            # Loop over all dihedrals in molecule1.
            for idx, dihedral in enumerate(dihedrals1.potentials()):
                # This dihedral contains an atom that is unique to molecule1.
                if info.atomIdx(dihedral.atom0()) in atoms1_idx or \
                   info.atomIdx(dihedral.atom1()) in atoms1_idx or \
                   info.atomIdx(dihedral.atom2()) in atoms1_idx or \
                   info.atomIdx(dihedral.atom3()) in atoms1_idx:
                       dihedral_idxs.append(idx)

            # Store the dihedral potentials for molecule1.
            potentials = dihedrals1.potentials()

            # Create the new set of dihedrals.
            dihedrals = _SireMM.FourAtomFunctions(edit_mol.info())

            # Add all of the dihedrals from molecule0.
            for dihedral in dihedrals0.potentials():
                dihedrals.set(dihedral.atom0(), dihedral.atom1(),
                    dihedral.atom2(), dihedral.atom3(), dihedral.function())

            # Loop over all of the additional dihedrals from molecule1.
            for idx in dihedral_idxs:
                # Extract the dihedral information.
                atom0 = info.atomIdx(potentials[idx].atom0())
                atom1 = info.atomIdx(potentials[idx].atom1())
                atom2 = info.atomIdx(potentials[idx].atom2())
                atom3 = info.atomIdx(potentials[idx].atom3())
                exprn = potentials[idx].function()

                # Map the atom indices to their position in the merged molecule.

                if atom0 in inv_mapping:
                    atom0 = inv_mapping[atom0]
                else:
                    atom0 = new_idx[atom0]

                if atom1 in inv_mapping:
                    atom1 = inv_mapping[atom1]
                else:
                    atom1 = new_idx[atom1]

                if atom2 in inv_mapping:
                    atom2 = inv_mapping[atom2]
                else:
                    atom2 = new_idx[atom2]

                if atom3 in inv_mapping:
                    atom3 = inv_mapping[atom3]
                else:
                    atom3 = new_idx[atom3]

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

            # Create a list to store the indices of impropers in molecule1 that
            # involve atoms that are part of the merged molecule at lambda = 0.
            improper_idxs = []

            # Get the molInfo object for molecule1.
            info = molecule1.info()

            # Loop over all impropers in molecule1.
            for idx, improper in enumerate(impropers1.potentials()):
                # This improper contains an atom that is unique to molecule1.
                if info.atomIdx(improper.atom0()) in atoms1_idx or \
                   info.atomIdx(improper.atom1()) in atoms1_idx or \
                   info.atomIdx(improper.atom2()) in atoms1_idx or \
                   info.atomIdx(improper.atom3()) in atoms1_idx:
                       improper_idxs.append(idx)

            # Store the improper potentials for molecule1.
            potentials = impropers1.potentials()

            # Create the new set of impropers.
            impropers = _SireMM.FourAtomFunctions(edit_mol.info())

            # Add all of the impropers from molecule0.
            for improper in impropers0.potentials():
                impropers.set(improper.atom0(), improper.atom1(),
                    improper.atom2(), improper.atom3(), improper.function())

            # Loop over all of the additional impropers from molecule1.
            for idx in improper_idxs:
                # Extract the improper information.
                atom0 = info.atomIdx(potentials[idx].atom0())
                atom1 = info.atomIdx(potentials[idx].atom1())
                atom2 = info.atomIdx(potentials[idx].atom2())
                atom3 = info.atomIdx(potentials[idx].atom3())
                exprn = potentials[idx].function()

                # Map the atom indices to their position in the merged molecule.

                if atom0 in inv_mapping:
                    atom0 = inv_mapping[atom0]
                else:
                    atom0 = new_idx[atom0]

                if atom1 in inv_mapping:
                    atom1 = inv_mapping[atom1]
                else:
                    atom1 = new_idx[atom1]

                if atom2 in inv_mapping:
                    atom2 = inv_mapping[atom2]
                else:
                    atom2 = new_idx[atom2]

                if atom3 in inv_mapping:
                    atom3 = inv_mapping[atom3]
                else:
                    atom3 = new_idx[atom3]

                # Set the new improper.
                impropers.set(atom0, atom1, atom2, atom3, exprn)

            # Add the impropers to the merged molecule.
            edit_mol.setProperty("improper0", impropers)

        ##############################
        # SET PROPERTIES AT LAMBDA = 1
        ##############################

        # Add the atom properties from molecule1.
        for atom in molecule1.atoms():
            # This atom has a match in molecule0.
            if atom.index() in inv_mapping:
                idx = inv_mapping[atom.index()]
            # This atom is unique to molecule, get its index in the merged molecule.
            else:
                idx = new_idx[atom.index()]

            # Loop over all atom properties.
            for prop in atom.propertyKeys():
                # Get the actual property name.
                if prop in inv_map0:
                    name = inv_map0[prop]
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
                if prop in inv_map1:
                    name = inv_map1[prop]
                else:
                    name = prop

                # Zero the "charge" and "LJ" property for atoms that are unique to molecule0.
                if name == "charge":
                    edit_mol = edit_mol.atom(atom.index()).setProperty("charge1", 0*_SireUnits.e_charge).molecule()
                elif name == "LJ":
                    lj = _SireMM.LJParameter(0*_SireUnits.angstrom, 0*_SireUnits.kcal_per_mol)
                    edit_mol = edit_mol.atom(atom.index()).setProperty("LJ1", lj).molecule()
                else:
                    # This is a perturbable property. Rename to "property1", e.g. "charge1".
                    if name in shared_props:
                        name = name + "0"

                    # Add the property to the atom in the merged molecule.
                    edit_mol = edit_mol.atom(atom.index()).setProperty(name, atom.property(prop)).molecule()

        # We now need to merge "bond", "angle", "dihedral", and "improper" parameters.
        # To do so, we extract the properties from molecule1, then add the additional
        # properties from molecule0, making sure to update the atom indices, and bond
        # atoms from molecule0 to the atoms to which they map in molecule1.

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

            # Create a list to store the indices of bonds in molecule0 that
            # involve atoms that are part of the merged molecule at lambda = 1.
            bond_idxs = []

            # Get the molInfo object for molecule0.
            info = molecule0.info()

            # Loop over all bonds in molecule0
            for idx, bond in enumerate(bonds0.potentials()):
                # This bond contains an atom that is unique to molecule0.
                if info.atomIdx(bond.atom0()) in atoms0_idx or \
                   info.atomIdx(bond.atom1()) in atoms0_idx:
                       bond_idxs.append(idx)

            # Store the bond potentials for molecule0.
            potentials = bonds0.potentials()

            # Create the new set of bonds.
            bonds = _SireMM.TwoAtomFunctions(edit_mol.info())

            # Add all of the bonds from molecule1.
            for bond in bonds1.potentials():
                # Extract the bond information.
                atom0 = info.atomIdx(bond.atom0())
                atom1 = info.atomIdx(bond.atom1())

                # Map the atom indices to their position in the merged molecule.

                if atom0 in inv_mapping:
                    atom0 = inv_mapping[atom0]
                else:
                    atom0 = new_idx[atom0]

                if atom1 in inv_mapping:
                    atom1 = inv_mapping[atom1]
                else:
                    atom1 = new_idx[atom1]

                # Set the new bond.
                bonds.set(atom0, atom1, exprn)

            # Loop over all of the additional bonds from molecule0.
            for idx in bond_idxs:
                # Extract the bond information.
                atom0 = info.atomIdx(potentials[idx].atom0())
                atom1 = info.atomIdx(potentials[idx].atom1())
                exprn = potentials[idx].function()

                # Set the new bond.
                bonds.set(atom0, atom1, exprn)

            # Add the bonds to the merged molecule.
            edit_mol.setProperty("bond1", bonds)

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

            # Create a list to store the indices of angles in molecule0 that
            # involve atoms that are part of the merged molecule at lambda = 1.
            angle_idxs = []

            # Get the molInfo object for molecule0.
            info = molecule0.info()

            # Loop over all angles in molecule0.
            for idx, angle in enumerate(angles0.potentials()):
                # This angle contains an atom that is unique to molecule0.
                if info.atomIdx(angle.atom0()) in atoms0_idx or \
                   info.atomIdx(angle.atom1()) in atoms0_idx or \
                   info.atomIdx(angle.atom2()) in atoms0_idx:
                       angle_idxs.append(idx)

            # Store the angle potentials for molecule0.
            potentials = angles0.potentials()

            # Create the new set of angles.
            angles = _SireMM.ThreeAtomFunctions(edit_mol.info())

            # Add all of the angles from molecule1.
            for angle in angles1.potentials():
                # Extract the angle information.
                atom0 = info.atomIdx(angle.atom0())
                atom1 = info.atomIdx(angle.atom1())
                atom2 = info.atomIdx(angle.atom2())
                exprn = angle.function()

                # Map the atom indices to their position in the merged molecule.

                if atom0 in inv_mapping:
                    atom0 = inv_mapping[atom0]
                else:
                    atom0 = new_idx[atom0]

                if atom1 in inv_mapping:
                    atom1 = inv_mapping[atom1]
                else:
                    atom1 = new_idx[atom1]

                if atom2 in inv_mapping:
                    atom2 = inv_mapping[atom2]
                else:
                    atom2 = new_idx[atom2]

                # Set the new angle.
                angles.set(atom0, atom1, atom2, exprn)

            # Loop over all of the additional angles from molecule0.
            for idx in angle_idxs:
                # Extract the angle information.
                atom0 = info.atomIdx(potentials[idx].atom0())
                atom1 = info.atomIdx(potentials[idx].atom1())
                atom2 = info.atomIdx(potentials[idx].atom2())
                exprn = potentials[idx].function()

                # Set the new angle.
                angles.set(atom0, atom1, atom2, exprn)

            # Add the angles to the merged molecule.
            edit_mol.setProperty("angle1", angles)

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

            # Create a list to store the indices of dihedrals in molecule0 that
            # involve atoms that are part of the merged molecule at lambda = 1.
            dihedral_idxs = []

            # Get the molInfo object for molecule0
            info = molecule0.info()

            # Loop over all dihedrals in molecule0.
            for idx, dihedral in enumerate(dihedrals0.potentials()):
                # This dihedral contains an atom that is unique to molecule0.
                if info.atomIdx(dihedral.atom0()) in atoms0_idx or \
                   info.atomIdx(dihedral.atom1()) in atoms0_idx or \
                   info.atomIdx(dihedral.atom2()) in atoms0_idx or \
                   info.atomIdx(dihedral.atom3()) in atoms0_idx:
                       dihedral_idxs.append(idx)

            # Store the dihedral potentials for molecule0.
            potentials = dihedrals0.potentials()

            # Create the new set of dihedrals.
            dihedrals = _SireMM.FourAtomFunctions(edit_mol.info())

            # Add all of the dihedrals from molecule1.
            for dihedral in dihedrals1.potentials():
                # Extract the dihedral information.
                atom0 = info.atomIdx(dihedral.atom0())
                atom1 = info.atomIdx(dihedral.atom1())
                atom2 = info.atomIdx(dihedral.atom2())
                atom3 = info.atomIdx(dihedral.atom3())
                exprn = dihedral.function()

                # Map the atom indices to their position in the merged molecule.

                if atom0 in inv_mapping:
                    atom0 = inv_mapping[atom0]
                else:
                    atom0 = new_idx[atom0]

                if atom1 in inv_mapping:
                    atom1 = inv_mapping[atom1]
                else:
                    atom1 = new_idx[atom1]

                if atom2 in inv_mapping:
                    atom2 = inv_mapping[atom2]
                else:
                    atom2 = new_idx[atom2]

                if atom3 in inv_mapping:
                    atom3 = inv_mapping[atom3]
                else:
                    atom3 = new_idx[atom3]

                # Set the new dihedral.
                dihedrals.set(atom0, atom1, atom2, atom3, exprn)

            # Loop over all of the additional dihedrals from molecule0.
            for idx in dihedral_idxs:
                # Extract the dihedral information.
                atom0 = info.atomIdx(potentials[idx].atom0())
                atom1 = info.atomIdx(potentials[idx].atom1())
                atom2 = info.atomIdx(potentials[idx].atom2())
                atom3 = info.atomIdx(potentials[idx].atom3())
                exprn = potentials[idx].function()

                # Set the new dihedral.
                dihedrals.set(atom0, atom1, atom2, atom3, exprn)

            # Add the dihedrals to the merged molecule.
            edit_mol.setProperty("dihedral1", dihedrals)

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

            # Create a list to store the indices of impropers in molecule0 that
            # involve atoms that are part of the merged molecule at lambda = 1.
            improper_idxs = []

            # Get the molInfo object for molecule0.
            info = molecule0.info()

            # Loop over all impropers in molecule0.
            for idx, improper in enumerate(impropers0.potentials()):
                # This improper contains an atom that is unique to molecule0.
                if info.atomIdx(improper.atom0()) in atoms0_idx or \
                   info.atomIdx(improper.atom1()) in atoms0_idx or \
                   info.atomIdx(improper.atom2()) in atoms0_idx or \
                   info.atomIdx(improper.atom3()) in atoms0_idx:
                       improper_idxs.append(idx)

            # Store the improper potentials for molecule0.
            potentials = impropers0.potentials()

            # Create the new set of impropers.
            impropers = _SireMM.FourAtomFunctions(edit_mol.info())

            # Add all of the impropers from molecule1.
            for improper in impropers1.potentials():
                # Extract the improper information.
                atom0 = info.atomIdx(improper.atom0())
                atom1 = info.atomIdx(improper.atom1())
                atom2 = info.atomIdx(improper.atom2())
                atom3 = info.atomIdx(improper.atom3())
                exprn = improper.function()

                # Map the atom indices to their position in the merged molecule.

                if atom0 in inv_mapping:
                    atom0 = inv_mapping[atom0]
                else:
                    atom0 = new_idx[atom0]

                if atom1 in inv_mapping:
                    atom1 = inv_mapping[atom1]
                else:
                    atom1 = new_idx[atom1]

                if atom2 in inv_mapping:
                    atom2 = inv_mapping[atom2]
                else:
                    atom2 = new_idx[atom2]

                if atom3 in inv_mapping:
                    atom3 = inv_mapping[atom3]
                else:
                    atom3 = new_idx[atom3]

                # Set the new improper.
                impropers.set(atom0, atom1, atom2, atom3, exprn)

            # Loop over all of the additional impropers from molecule0.
            for idx in improper_idxs:
                # Extract the improper information.
                atom0 = info.atomIdx(potentials[idx].atom0())
                atom1 = info.atomIdx(potentials[idx].atom1())
                atom2 = info.atomIdx(potentials[idx].atom2())
                atom3 = info.atomIdx(potentials[idx].atom3())
                exprn = potentials[idx].function()

                # Set the new improper.
                impropers.set(atom0, atom1, atom2, atom3, exprn)

            # Add the impropers to the merged molecule.
            edit_mol.setProperty("improper1", impropers)

        # Return the merged molecule.
        return edit_mol.commit()

# Import at bottom of module to avoid circular dependency.
from ._molecule import Molecule as _Molecule
