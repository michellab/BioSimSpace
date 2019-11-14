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
A thin wrapper around Sire.Mol.Molecule. This is an internal package and should
not be directly exposed to the user.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["Molecule"]

from pytest import approx as _approx

import os.path as _path
import random as _random
import string as _string

from Sire import Base as _SireBase
from Sire import CAS as _SireCAS
from Sire import MM as _SireMM
from Sire import Mol as _SireMol
from Sire import Units as _SireUnits

from BioSimSpace import _isVerbose
from BioSimSpace._Exceptions import IncompatibleError as _IncompatibleError
from BioSimSpace.Types import Length as _Length

from ._sire_wrapper import SireWrapper as _SireWrapper

class Molecule(_SireWrapper):
    """A container class for storing a molecule."""

    def __init__(self, molecule):
        """Constructor.

           Parameters
           ----------

           molecule : Sire.Mol.Molecule, :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
               A Sire or BioSimSpace Molecule object.
        """

        # Set the force field variable. This records the force field with which
        # the molecule has been parameterised, i.e. by BSS.Parameters.
        self._forcefield = None

        # Set the components of the merged molecule to None.
        self._molecule0 = None
        self._molecule1 = None

        # Check that the molecule is valid.

        # A Sire Molecule object.
        if type(molecule) is _SireMol.Molecule:
            super().__init__(molecule)
            if self._sire_object.hasProperty("is_perturbable"):
                self._convertFromMergedMolecule()

        # Another BioSimSpace Molecule object.
        elif type(molecule) is Molecule:
            super().__init__(molecule._sire_object)
            if molecule._molecule0 is not None:
                self._molecule0 = Molecule(molecule._molecule0)
            if molecule._molecule1 is not None:
                self._molecule1 = Molecule(molecule._molecule1)
            self._forcefield = molecule._forcefield
            self._is_merged = molecule._is_merged

        # Invalid type.
        else:
            raise TypeError("'molecule' must be of type 'Sire.Mol.Molecule' "
                            "or 'BioSimSpace._SireWrappers.Molecule'.")

        # Flag that this object holds multiple atoms.
        self._is_multi_atom = True

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

        # Whether other is a container of molecules.
        is_sire_container = False

        # Validate the input.

        molecules = [self]

        # A System object.
        if type(other) is _System:
            system = _System(other)
            system.addMolecules(self)
            return system

        # A single Molecule object.
        elif type(other) is Molecule:
            molecules.append(other)

        # A Molecules object.
        elif type(other) is _Molecules:
            is_sire_container = True

        # A list of Molecule objects.
        elif type(other) is list and all(isinstance(x, Molecule) for x in other):
            molecules.extend(other)

        # Unsupported.
        else:
            raise TypeError("'other' must be of type 'BioSimSpace._SireWrappers.System', "
                            "'BioSimSpace._SireWrappers.Molecule', 'BioSimSpace._SireWrappers.Molecules' "
                            "or a list of 'BioSimSpace._SireWrappers.Molecule' types")

        # Add this molecule to the container and return.
        if is_sire_container:
            return other + self

        # Create a new molecule group to store the molecules, then create a
        # container and return.
        else:
            molgrp = _SireMol.MoleculeGroup("all")

            for molecule in molecules:
                molgrp.add(molecule._sire_object)

            return _Molecules(molgrp.molecules())

    def number(self):
        """Return the number of the molecule. Each molecule has a unique
           identification number.

           Returns
           -------

           mol_num : int
               The unique number of the molecule.
        """
        return self._sire_object.number().value()

    def getResidues(self):
        """Return a list containing the residues in the molecule.

           Returns
           -------

           residues : [:class:`Residue <BioSimSpace._SireWrappers.Residue>`]
               The list of residues in the molecule.
        """
        residues = []
        for residue in self._sire_object.residues():
            residues.append(_Residue(residue))
        return residues

    def getAtoms(self):
        """Return a list containing the atoms in the molecule.

           Returns
           -------

           atoms : [:class:`Atom <BioSimSpace._SireWrappers.Atom>`]
               The list of atoms in the molecule.
        """
        atoms = []
        for atom in self._sire_object.atoms():
            atoms.append(_Atom(atom))
        return atoms

    def molecule0(self):
        """Return the component of the merged molecule at lambda = 0.

           Returns
           -------

           molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
               The component of the merged molecule at lambda = 0.
               Returns None if this isn't a merged molecule.
        """
        return self._molecule0

    def molecule1(self):
        """Return the component of the merged molecule at lambda = 1.

           Returns
           -------

           molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
               The component of the merged molecule at lambda = 1.
               Returns None if this isn't a merged molecule.
        """
        return self._molecule1

    def nAtoms(self):
        """Return the number of atoms in the molecule.

           Returns
           -------

           num_atoms : int
               The number of atoms in the molecule.
        """
        return self._sire_object.nAtoms()

    def nResidues(self):
        """Return the number of residues in the molecule.

           Returns
           -------

           num_residues : int
               The number of residues in the molecule.
        """
        return self._sire_object.nResidues()

    def nChains(self):
        """Return the number of chains in the molecule.

           Returns
           -------

           num_chains : int
               The number of chains in the molecule.
        """
        return self._sire_object.nChains()

    def isMerged(self):
        """Whether this molecule has been merged with another.

           Returns
           -------

           is_merged : bool
               Whether the molecule has been merged.
        """
        return self._is_merged

    def isWater(self):
        """Whether this is a water molecule.

           Returns
           -------

           is_water : bool
               Whether this is a water molecule.
        """
        # Water models have 5 or less atoms.
        if self.nAtoms() > 5:
            return False

        # Tally counters for the number of H and O atoms.
        num_hydrogen = 0
        num_oxygen = 0

        # Loop over all atoms in the molecule.
        for atom in self._sire_object.atoms():

            # First try using the "element" property of the atom.
            try:
                # Hydrogen.
                if atom.property("element") == _SireMol.Element("H"):
                    num_hydrogen += 1
                # Oxygen.
                elif atom.property("element") == _SireMol.Element("O"):
                    num_oxygen += 1

            # Otherwise, try to infer the element from the atom name.
            except:
                # Strip all digits from the name.
                name = "".join([x for x in atom.name().value() if not x.isdigit()])

                # Remove any whitespace.
                name = name.replace(" ", "")

                # Try to infer the element.
                element = _SireMol.Element.biologicalElement(name)

                # Hydrogen.
                if element == _SireMol.Element("H"):
                    num_hydrogen += 1
                # Oxygen.
                elif element == _SireMol.Element("O"):
                    num_oxygen += 1

        # A water molecule has two Hydrogens and one Oxygen.
        if num_hydrogen == 2 and num_oxygen == 1:
            return True
        else:
            return False

    def toSystem(self):
        """Convert a single Molecule to a System.

           Returns
           -------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
        """
        return _System(self)

    def search(self, query):
        """Search the molecule for atoms and residues. Search results will be
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

           >>> result = molecule.search("resname ALA")

           Search for all oxygen or hydrogen atoms.

           >>> result = molecule.search("element oxygen or element hydrogen")

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

    def _getPropertyMap0(self):
        """Generate a property map for the lambda = 0 state of the merged molecule."""

        property_map = {}

        if self._is_merged:
            for prop in self._sire_object.propertyKeys():
                if prop[-1] == "0":
                    property_map[prop[:-1]] = prop

        return property_map

    def _getPropertyMap1(self):
        """Generate a property map for the lambda = 1 state of the merged molecule."""

        property_map = {}

        if self._is_merged:
            for prop in self._sire_object.propertyKeys():
                if prop[-1] == "1":
                    property_map[prop[:-1]] = prop

        return property_map

    def _makeCompatibleWith(self, molecule, property_map={}, overwrite=True,
            rename_atoms=False, verbose=False):
        """Make this molecule compatible with passed one, i.e. match atoms and
           add all additional properties.

           Parameters
           ----------

           molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
               The molecule to match with.

           property_map : dict
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
            mol1 = molecule._sire_object
        else:
            raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule', or 'Sire.Mol.Molecule'")

        if type(property_map) is not dict:
            raise TypeError("'property_map' must be of type 'dict'")

        if type(overwrite) is not bool:
            raise TypeError("'overwrite' must be of type 'bool'")

        if type(rename_atoms) is not bool:
            raise TypeError("'rename_atoms' must be of type 'bool'")

        if type(verbose) is not bool:
            raise TypeError("'verbose' must be of type 'bool'")

        # Get the two Sire molecules.
        mol0 = self._sire_object

        # Store the number of atoms to match.
        num_atoms = mol0.nAtoms()

        # The new molecule must have at least as many atoms.
        if mol1.nAtoms() < num_atoms:
            raise _IncompatibleError("The passed molecule is incompatible with the original! "
                                     "self.nAtoms() = %d, other.nAtoms() = %d" % (num_atoms, mol1.nAtoms()))

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
        _property_map = property_map.copy()

        # See if any of the new properties are in the map, add them if not.
        for prop in props0:
            if not prop in _property_map:
                _property_map[prop] = prop
        for prop in props1:
            if not prop in _property_map:
                _property_map[prop] = prop

        # Make the molecule editable.
        edit_mol = mol0.edit()

        if "parameters" in _property_map:
            param = _property_map["parameters"]
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
                    if (not mol0.hasProperty(_property_map[prop])) or overwrite:
                        if verbose:
                            print("  %s" % _property_map[prop])
                        try:
                            edit_mol = edit_mol.setProperty(_property_map[prop], mol1.property(prop))
                        except Exception as e:
                            msg = "Failed to set property '%s'" % _property_map[prop]
                            if _isVerbose():
                                raise _IncompatibleError(msg) from e
                            else:
                                raise _IncompatibleError(msg) from None

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
                if (not mol0.hasProperty(_property_map[prop])) or overwrite:
                    # Loop over all of the atom mapping pairs and set the property.
                    for idx0, idx1 in matches.items():
                        # Does the atom have this property?
                        # If so, add it to the matching atom in this molecule.
                        if mol1.atom(idx1).hasProperty(prop):
                            if verbose:
                                print("  %-20s %s --> %s" % (_property_map[prop], idx1, idx0))
                            try:
                                edit_mol = edit_mol.atom(idx0).setProperty(_property_map[prop], mol1.atom(idx1).property(prop)).molecule()
                                seen_prop[prop] = True
                            except Exception as e:
                                msg = "Failed to copy property '%s' from %s to %s." % (_property_map[prop], idx1, idx0)
                                if _isVerbose():
                                    raise _IncompatibleError(msg) from e
                                else:
                                    raise _IncompatibleError(msg) from None

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
                        if (not mol0.hasProperty(_property_map[prop])) or overwrite:
                            if verbose:
                                print("  %s" % _property_map[prop])

                            # Get the property from the parameterised molecule.
                            propty = mol1.property(_property_map[prop])

                            # Try making it compatible with the original molecule.
                            if hasattr(propty, "makeCompatibleWith"):
                                try:
                                    propty = propty.makeCompatibleWith(mol0, matches)
                                except Exception as e:
                                    msg = "Incompatible property: %s" % _property_map[prop]
                                    if _isVerbose():
                                        raise _IncompatibleError(msg) from e
                                    else:
                                        raise _IncompatibleError(msg) from None

                            # Now try to set the property.
                            edit_mol.setProperty(_property_map[prop], propty)

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
                except Exception as e:
                    msg = "Failed to rename atom: %s --> %s" % (name0, name1)
                    if _isVerbose():
                        raise _IncompatibleError(msg) from e
                    else:
                        raise _IncompatibleError(msg) from None

        # Commit the changes.
        self._sire_object = edit_mol.commit()

    def _convertFromMergedMolecule(self):
        """Convert from a merged molecule."""

        # Extract the components of the merged molecule.
        try:
            mol0 = self._sire_object.property("molecule0")
            mol1 = self._sire_object.property("molecule1")
        except:
            raise _IncompatibleError("The merged molecule doesn't have the required properties!")

        # Store the components.
        self._molecule0 = Molecule(mol0)
        self._molecule1 = Molecule(mol1)

        # Flag that the molecule is merged.
        self._is_merged = True

    def _fixCharge(self, property_map={}):
        """Make the molecular charge an integer value.

           Parameters
           ----------

           property_map : dict
               A dictionary that maps "properties" in this molecule to their
               user defined values. This allows the user to refer to properties
               with their own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Get the user defined charge property.
        prop = property_map.get("charge", "charge")

        if not self._sire_object.hasProperty(prop):
            raise _IncompatibleError("Molecule does not have charge property: '%s'." % prop)

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
        edit_mol = self._sire_object.edit()

        # Shift the charge of each atom in the molecule by delta.
        # Make sure to invert the sign of the charge since it is in
        # units of -1 |e|.
        for atom in edit_mol.atoms():
            charge = edit_mol.atom(atom.index()).property(prop).value()
            charge = -(charge + delta)
            edit_mol = edit_mol.atom(atom.index()) \
                               .setProperty(prop, charge * _SireUnits.e_charge) \
                               .molecule()

        # Update the Sire molecule.
        self._sire_object = edit_mol.commit()

    def _fromPertFile(self, filename):
        """Create a merged molecule from a perturbation file.

           Parameters
           ----------

           filename: str
               The location of the perturbation file.
        """

        raise NotImplementedError("Perturbation file reader is not yet implemented.")

        if not _path.isfile(filename):
            raise IOError("Perturbation file doesn't exist: '%s'" % filename)

    def _toPertFile(self, filename="MORPH.pert", zero_dummy_dihedrals=False,
            zero_dummy_impropers=False, print_all_atoms=False, property_map={}):
        """Write the merged molecule to a perturbation file.

           Parameters
           ----------

           filename : str
               The name of the perturbation file.

           zero_dummy_dihedrals : bool
               Whether to zero the barrier height for dihedrals involving
               dummy atoms.

           zero_dummy_impropers : bool
               Whether to zero the barrier height for impropers involving
               dummy atoms.

           print_all_atoms : bool
               Whether to print all atom records to the pert file, not just
               the atoms that are perturbed.

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }

           Returns
           -------

           molecule : Sire.Mol.Molecule
               The molecule with properties corresponding to the lamda = 0 state.
        """

        if not self._is_merged:
            raise _IncompatibleError("This isn't a merged molecule. Cannot write perturbation file!")

        if not self._sire_object.property("forcefield0").isAmberStyle():
            raise _IncompatibleError("Can only write perturbation files for AMBER style force fields.")

        if type(zero_dummy_dihedrals) is not bool:
            raise TypeError("'zero_dummy_dihedrals' must be of type 'bool'")

        if type(zero_dummy_impropers) is not bool:
            raise TypeError("'zero_dummy_impropers' must be of type 'bool'")

        if type(print_all_atoms) is not bool:
            raise TypeError("'print_all_atoms' must be of type 'bool'")

        if type(property_map) is not dict:
            raise TypeError("'property_map' must be of type 'dict'")

        # Extract and copy the Sire molecule.
        mol = self._sire_object.__deepcopy__()

        # First work out the indices of atoms that are perturbed.
        pert_idxs = []

        # Perturbed atoms change one of the following properties:
        # "ambertype", "LJ", or "charge".
        for atom in mol.atoms():
            if atom.property("ambertype0") != atom.property("ambertype1") or \
               atom.property("LJ0") != atom.property("LJ1")               or \
               atom.property("charge0") != atom.property("charge1"):
                pert_idxs.append(atom.index())

        # The pert file uses atom names for identification purposes. This means
        # that the names must be unique. As such we need to count the number of
        # atoms with a particular name, then append an index to their name.

        # A dictionary to track the atom names.
        atom_names = {}

        # Loop over all atoms in the molecule.
        for atom in mol.atoms():
            if atom.name() in atom_names:
                atom_names[atom.name()] += 1
            else:
                atom_names[atom.name()] = 1

        # Create a set from the atoms names seen so far.
        names = set(atom_names.keys())

        # If there are duplicate names, then we need to rename the atoms.
        if sum(atom_names.values()) > len(names):

            # Make the molecule editable.
            edit_mol = mol.edit()

            # Create a dictionary to flag whether we've seen each atom name.
            is_seen = { name : False for name in names }

            # Tally counter for the number of dummy atoms.
            num_dummy = 1

            # Loop over all atoms.
            for atom in mol.atoms():
                # Store the original atom.
                name = atom.name()

                # If this is a dummy atom, then rename it as "DU##", where ## is a
                # two-digit number padded with a leading zero.
                if atom.property("element0") == _SireMol.Element("X"):
                    # Create the new atom name.
                    new_name = "DU%02d" % num_dummy

                    # Convert to an AtomName and rename the atom.
                    new_name = _SireMol.AtomName(new_name)
                    edit_mol = edit_mol.atom(atom.index()).rename(new_name).molecule()

                    # Update the number of dummy atoms that have been named.
                    num_dummy += 1

                    # Since ligands typically have less than 100 atoms, the following
                    # exception shouldn't be triggered. We can't support perturbations
                    # with 100 or more dummy atoms in the lambda = 0 state because of
                    # AMBER fixed width atom naming restrictions (4 character width).
                    # We could give dummies a unique name in the same way that non-dummy
                    # atoms are handled (see else) block below, but instead we'll raise
                    # an exception.
                    if num_dummy == 100:
                        raise RuntimeError("Dummy atom naming limit exceeded! (100 atoms)")

                    # Append to the list of seen names.
                    names.add(new_name)

                else:
                    # There is more than one atom with this name, and this is the second
                    # time we've come across it.
                    if atom_names[name] > 1 and is_seen[name]:
                        # Create the base of the new name.
                        new_name = name.value()

                        # Create a random suffix.
                        suffix = _random_suffix(new_name)

                        # Zero the number of attempted renamings.
                        num_attempts = 0

                        # If this name already exists, keep trying until we get a unique name.
                        while new_name + suffix in names:
                            suffix = _random_suffix(new_name)
                            num_attempts += 1

                            # Abort if we've tried more than 100 times.
                            if num_attempts == 100:
                                raise RuntimeError("Error while writing SOMD pert file. "
                                                   "Unable to generate a unique suffix for "
                                                   "atom name: '%s'" % new_name)

                        # Append the suffix to the name and store in the set of seen names.
                        new_name = new_name + suffix
                        names.add(new_name)

                        # Convert to an AtomName and rename the atom.
                        new_name = _SireMol.AtomName(new_name)
                        edit_mol = edit_mol.atom(atom.index()).rename(new_name).molecule()

                    # Record that we've seen this atom name.
                    is_seen[name] = True

            # Store the updated molecule.
            mol = edit_mol.commit()

        # Now write the perturbation file.

        with open(filename, "w") as file:
            # Get the info object for the molecule.
            info = mol.info()

            # Write the version header.
            file.write("version 1\n")

            # Start molecule record.
            file.write("molecule LIG\n")

            # 1) Atoms.

            # Print all atom records.
            if print_all_atoms:
                for atom in mol.atoms():
                    # Start atom record.
                    file.write("    atom\n")

                    # Get the initial/final Lennard-Jones properties.
                    LJ0 = atom.property("LJ0");
                    LJ1 = atom.property("LJ1");

                    # Atom data.
                    file.write("        name           %s\n" % atom.name().value())
                    file.write("        initial_type   %s\n" % atom.property("ambertype0"))
                    file.write("        final_type     %s\n" % atom.property("ambertype1"))
                    file.write("        initial_LJ     %.5f %.5f\n" % (LJ0.sigma().value(), LJ0.epsilon().value()))
                    file.write("        final_LJ       %.5f %.5f\n" % (LJ1.sigma().value(), LJ1.epsilon().value()))
                    file.write("        initial_charge %.5f\n" % atom.property("charge0").value())
                    file.write("        final_charge   %.5f\n" % atom.property("charge1").value())

                    # End atom record.
                    file.write("    endatom\n")

            # Only print records for the atoms that are perturbed.
            else:
                for idx in pert_idxs:
                    # Get the perturbed atom.
                    atom = mol.atom(idx)

                    # Start atom record.
                    file.write("    atom\n")

                    # Get the initial/final Lennard-Jones properties.
                    LJ0 = atom.property("LJ0");
                    LJ1 = atom.property("LJ1");

                    # Atom data.
                    file.write("        name           %s\n" % atom.name().value())
                    file.write("        initial_type   %s\n" % atom.property("ambertype0"))
                    file.write("        final_type     %s\n" % atom.property("ambertype1"))
                    file.write("        initial_LJ     %.5f %.5f\n" % (LJ0.sigma().value(), LJ0.epsilon().value()))
                    file.write("        final_LJ       %.5f %.5f\n" % (LJ1.sigma().value(), LJ1.epsilon().value()))
                    file.write("        initial_charge %.5f\n" % atom.property("charge0").value())
                    file.write("        final_charge   %.5f\n" % atom.property("charge1").value())

                    # End atom record.
                    file.write("    endatom\n")

            # 2) Bonds.

            # Extract the bonds at lambda = 0 and 1.
            bonds0 = mol.property("bond0").potentials()
            bonds1 = mol.property("bond1").potentials()

            # Dictionaries to store the BondIDs at lambda = 0 and 1.
            bonds0_idx = {}
            bonds1_idx = {}

            # Loop over all bonds at lambda = 0.
            for idx, bond in enumerate(bonds0):
                # Get the AtomIdx for the atoms in the bond.
                idx0 = info.atomIdx(bond.atom0())
                idx1 = info.atomIdx(bond.atom1())

                # Create the BondID.
                bond_id = _SireMol.BondID(idx0, idx1)

                # Add to the list of ids.
                bonds0_idx[bond_id] = idx

            # Loop over all bonds at lambda = 1.
            for idx, bond in enumerate(bonds1):
                # Get the AtomIdx for the atoms in the bond.
                idx0 = info.atomIdx(bond.atom0())
                idx1 = info.atomIdx(bond.atom1())

                # Create the AngleID.
                bond_id = _SireMol.BondID(idx0, idx1)

                # Add to the list of ids.
                if bond_id.mirror() in bonds0_idx:
                    bonds1_idx[bond_id.mirror()] = idx
                else:
                    bonds1_idx[bond_id] = idx

            # Now work out the BondIDs that are unique at lambda = 0 and 1
            # as well as those that are shared.
            bonds0_unique_idx = {}
            bonds1_unique_idx = {}
            bonds_shared_idx = {}

            # lambda = 0.
            for idx in bonds0_idx.keys():
                if idx not in bonds1_idx.keys():
                    bonds0_unique_idx[idx] = bonds0_idx[idx]
                else:
                    bonds_shared_idx[idx] = (bonds0_idx[idx], bonds1_idx[idx])

            # lambda = 1.
            for idx in bonds1_idx.keys():
                if idx not in bonds0_idx.keys():
                    bonds1_unique_idx[idx] = bonds1_idx[idx]
                elif idx not in bonds_shared_idx.keys():
                    bonds_shared_idx[idx] = (bonds0_idx[idx], angles1_idx[idx])

            # First create records for the bonds that are unique to lambda = 0 and 1.

            # lambda = 0.
            for idx in bonds0_unique_idx.values():
                # Get the bond potential.
                bond = bonds0[idx]

                # Get the AtomIdx for the atoms in the bond.
                idx0 = info.atomIdx(bond.atom0())
                idx1 = info.atomIdx(bond.atom1())

                # Cast the function as an AmberBond.
                amber_bond = _SireMM.AmberBond(bond.function(), _SireCAS.Symbol("r"))

                # Start angle record.
                file.write("    bond\n")

                # Bond data.
                file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                file.write("        initial_force  %.5f\n" % amber_bond.k())
                file.write("        initial_equil  %.5f\n" % amber_bond.r0())
                file.write("        final_force    %.5f\n" % 0.0)
                file.write("        final_equil    %.5f\n" % amber_bond.r0())

                # End bond record.
                file.write("    endbond\n")

            # lambda = 1.
            for idx in bonds1_unique_idx.values():
                # Get the bond potential.
                bond = bonds1[idx]

                # Get the AtomIdx for the atoms in the bond.
                idx0 = info.atomIdx(bond.atom0())
                idx1 = info.atomIdx(bond.atom1())

                # Cast the function as an AmberBond.
                amber_bond = _SireMM.AmberBond(bond.function(), _SireCAS.Symbol("r"))

                # Start angle record.
                file.write("    bond\n")

                # Bond data.
                file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                file.write("        initial_force  %.5f\n" % 0.0)
                file.write("        initial_equil  %.5f\n" % amber_bond.r0())
                file.write("        final_force    %.5f\n" % amber_bond.k())
                file.write("        final_equil    %.5f\n" % amber_bond.r0())

                # End bond record.
                file.write("    endbond\n")

            # Now add records for the shared bonds.
            for idx0, idx1 in bonds_shared_idx.values():
                # Get the bond potentials.
                bond0 = bonds0[idx0]
                bond1 = bonds1[idx1]

                # Get the AtomIdx for the atoms in the bond.
                idx0 = info.atomIdx(bond0.atom0())
                idx1 = info.atomIdx(bond0.atom1())

                # Check that an atom in the bond is perturbed.
                if _has_pert_atom([idx0, idx1], pert_idxs):

                    # Cast the bonds as AmberBonds.
                    amber_bond0 = _SireMM.AmberBond(bond0.function(), _SireCAS.Symbol("r"))
                    amber_bond1 = _SireMM.AmberBond(bond1.function(), _SireCAS.Symbol("r"))

                    # Check whether a dummy atoms are present in the lambda = 0
                    # and lambda = 1 states.
                    initial_dummy = _has_dummy(mol, [idx0, idx1])
                    final_dummy = _has_dummy(mol, [idx0, idx1], True)

                    # Cannot have a bond with a dummy in both states.
                    if initial_dummy and final_dummy:
                        raise _IncompatibleError("Dummy atoms are present in both the initial "
                                                 "and final bond?")

                    # Set the bond parameters of the dummy state to those of the non-dummy end state.
                    if initial_dummy or final_dummy:
                        has_dummy = True
                        if initial_dummy:
                            amber_bond0 = amber_bond1
                        else:
                            amber_bond1 = amber_bond0
                    else:
                        has_dummy = False

                    # Only write record if the bond parameters change.
                    if has_dummy or amber_bond0 != amber_bond1:

                        # Start bond record.
                        file.write("    bond\n")

                        # Angle data.
                        file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                        file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                        file.write("        initial_force  %.5f\n" % amber_bond0.k())
                        file.write("        initial_equil  %.5f\n" % amber_bond0.r0())
                        file.write("        final_force    %.5f\n" % amber_bond1.k())
                        file.write("        final_equil    %.5f\n" % amber_bond1.r0())

                        # End bond record.
                        file.write("    endbond\n")

            # 3) Angles.

            # Extract the angles at lambda = 0 and 1.
            angles0 = mol.property("angle0").potentials()
            angles1 = mol.property("angle1").potentials()

            # Dictionaries to store the AngleIDs at lambda = 0 and 1.
            angles0_idx = {}
            angles1_idx = {}

            # Loop over all angles at lambda = 0.
            for idx, angle in enumerate(angles0):
                # Get the AtomIdx for the atoms in the angle.
                idx0 = info.atomIdx(angle.atom0())
                idx1 = info.atomIdx(angle.atom1())
                idx2 = info.atomIdx(angle.atom2())

                # Create the AngleID.
                angle_id = _SireMol.AngleID(idx0, idx1, idx2)

                # Add to the list of ids.
                angles0_idx[angle_id] = idx

            # Loop over all angles at lambda = 1.
            for idx, angle in enumerate(angles1):
                # Get the AtomIdx for the atoms in the angle.
                idx0 = info.atomIdx(angle.atom0())
                idx1 = info.atomIdx(angle.atom1())
                idx2 = info.atomIdx(angle.atom2())

                # Create the AngleID.
                angle_id = _SireMol.AngleID(idx0, idx1, idx2)

                # Add to the list of ids.
                if angle_id.mirror() in angles0_idx:
                    angles1_idx[angle_id.mirror()] = idx
                else:
                    angles1_idx[angle_id] = idx

            # Now work out the AngleIDs that are unique at lambda = 0 and 1
            # as well as those that are shared.
            angles0_unique_idx = {}
            angles1_unique_idx = {}
            angles_shared_idx = {}

            # lambda = 0.
            for idx in angles0_idx.keys():
                if idx not in angles1_idx.keys():
                    angles0_unique_idx[idx] = angles0_idx[idx]
                else:
                    angles_shared_idx[idx] = (angles0_idx[idx], angles1_idx[idx])

            # lambda = 1.
            for idx in angles1_idx.keys():
                if idx not in angles0_idx.keys():
                    angles1_unique_idx[idx] = angles1_idx[idx]
                elif idx not in angles_shared_idx.keys():
                    angles_shared_idx[idx] = (angles0_idx[idx], angles1_idx[idx])

            # First create records for the angles that are unique to lambda = 0 and 1.

            # lambda = 0.
            for idx in angles0_unique_idx.values():
                # Get the angle potential.
                angle = angles0[idx]

                # Get the AtomIdx for the atoms in the angle.
                idx0 = info.atomIdx(angle.atom0())
                idx1 = info.atomIdx(angle.atom1())
                idx2 = info.atomIdx(angle.atom2())

                # Cast the function as an AmberAngle.
                amber_angle = _SireMM.AmberAngle(angle.function(), _SireCAS.Symbol("theta"))

                # Start angle record.
                file.write("    angle\n")

                # Angle data.
                file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                file.write("        atom2          %s\n" % mol.atom(idx2).name().value())
                file.write("        initial_force  %.5f\n" % amber_angle.k())
                file.write("        initial_equil  %.5f\n" % amber_angle.theta0())
                file.write("        final_force    %.5f\n" % 0.0)
                file.write("        final_equil    %.5f\n" % amber_angle.theta0())

                # End angle record.
                file.write("    endangle\n")

            # lambda = 1.
            for idx in angles1_unique_idx.values():
                # Get the angle potential.
                angle = angles1[idx]

                # Get the AtomIdx for the atoms in the angle.
                idx0 = info.atomIdx(angle.atom0())
                idx1 = info.atomIdx(angle.atom1())
                idx2 = info.atomIdx(angle.atom2())

                # Cast the function as an AmberAngle.
                amber_angle = _SireMM.AmberAngle(angle.function(), _SireCAS.Symbol("theta"))

                # Start angle record.
                file.write("    angle\n")

                # Angle data.
                file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                file.write("        atom2          %s\n" % mol.atom(idx2).name().value())
                file.write("        initial_force  %.5f\n" % 0.0)
                file.write("        initial_equil  %.5f\n" % amber_angle.theta0())
                file.write("        final_force    %.5f\n" % amber_angle.k())
                file.write("        final_equil    %.5f\n" % amber_angle.theta0())

                # End angle record.
                file.write("    endangle\n")

            # Now add records for the shared angles.
            for idx0, idx1 in angles_shared_idx.values():
                # Get the angle potentials.
                angle0 = angles0[idx0]
                angle1 = angles1[idx1]

                # Get the AtomIdx for the atoms in the angle.
                idx0 = info.atomIdx(angle0.atom0())
                idx1 = info.atomIdx(angle0.atom1())
                idx2 = info.atomIdx(angle0.atom2())

                # Check that an atom in the angle is perturbed.
                if _has_pert_atom([idx0, idx1, idx2], pert_idxs):

                    # Cast the functions as AmberAngles.
                    amber_angle0 = _SireMM.AmberAngle(angle0.function(), _SireCAS.Symbol("theta"))
                    amber_angle1 = _SireMM.AmberAngle(angle1.function(), _SireCAS.Symbol("theta"))

                    # Check whether a dummy atoms are present in the lambda = 0
                    # and lambda = 1 states.
                    initial_dummy = _has_dummy(mol, [idx0, idx1, idx2])
                    final_dummy = _has_dummy(mol, [idx0, idx1, idx2], True)

                    # Set the angle parameters of the dummy state to those of the non-dummy end state.
                    if initial_dummy and final_dummy:
                        has_dummy = True
                        amber_angle0 = _SireMM.AmberAngle()
                        amber_angle1 = _SireMM.AmberAngle()
                    elif initial_dummy or final_dummy:
                        has_dummy = True
                        if initial_dummy:
                            amber_angle0 = amber_angle1
                        else:
                            amber_angle1 = amber_angle0
                    else:
                        has_dummy = False

                    # Only write record if the angle parameters change.
                    if has_dummy or amber_angle0 != amber_angle1:

                        # Start angle record.
                        file.write("    angle\n")

                        # Angle data.
                        file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                        file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                        file.write("        atom2          %s\n" % mol.atom(idx2).name().value())
                        file.write("        initial_force  %.5f\n" % amber_angle0.k())
                        file.write("        initial_equil  %.5f\n" % amber_angle0.theta0())
                        file.write("        final_force    %.5f\n" % amber_angle1.k())
                        file.write("        final_equil    %.5f\n" % amber_angle1.theta0())

                        # End angle record.
                        file.write("    endangle\n")

            # 4) Dihedrals.

            # Extract the dihedrals at lambda = 0 and 1.
            dihedrals0 = mol.property("dihedral0").potentials()
            dihedrals1 = mol.property("dihedral1").potentials()

            # Dictionaries to store the DihedralIDs at lambda = 0 and 1.
            dihedrals0_idx = {}
            dihedrals1_idx = {}

            # Loop over all dihedrals at lambda = 0.
            for idx, dihedral in enumerate(dihedrals0):
                # Get the AtomIdx for the atoms in the dihedral.
                idx0 = info.atomIdx(dihedral.atom0())
                idx1 = info.atomIdx(dihedral.atom1())
                idx2 = info.atomIdx(dihedral.atom2())
                idx3 = info.atomIdx(dihedral.atom3())

                # Create the DihedralID.
                dihedral_id = _SireMol.DihedralID(idx0, idx1, idx2, idx3)

                # Add to the list of ids.
                dihedrals0_idx[dihedral_id] = idx

            # Loop over all dihedrals at lambda = 1.
            for idx, dihedral in enumerate(dihedrals1):
                # Get the AtomIdx for the atoms in the dihedral.
                idx0 = info.atomIdx(dihedral.atom0())
                idx1 = info.atomIdx(dihedral.atom1())
                idx2 = info.atomIdx(dihedral.atom2())
                idx3 = info.atomIdx(dihedral.atom3())

                # Create the DihedralID.
                dihedral_id = _SireMol.DihedralID(idx0, idx1, idx2, idx3)

                # Add to the list of ids.
                if dihedral_id.mirror() in dihedrals0_idx:
                    dihedrals1_idx[dihedral_id.mirror()] = idx
                else:
                    dihedrals1_idx[dihedral_id] = idx

            # Now work out the DihedralIDs that are unique at lambda = 0 and 1
            # as well as those that are shared.
            dihedrals0_unique_idx = {}
            dihedrals1_unique_idx = {}
            dihedrals_shared_idx = {}

            # lambda = 0.
            for idx in dihedrals0_idx.keys():
                if idx not in dihedrals1_idx.keys():
                    dihedrals0_unique_idx[idx] = dihedrals0_idx[idx]
                else:
                    dihedrals_shared_idx[idx] = (dihedrals0_idx[idx], dihedrals1_idx[idx])

            # lambda = 1.
            for idx in dihedrals1_idx.keys():
                if idx not in dihedrals0_idx.keys():
                    dihedrals1_unique_idx[idx] = dihedrals1_idx[idx]
                elif idx not in dihedrals_shared_idx.keys():
                    dihedrals_shared_idx[idx] = (dihedrals0_idx[idx], dihedrals1_idx[idx])

            # First create records for the dihedrals that are unique to lambda = 0 and 1.

            # lambda = 0.
            for idx in dihedrals0_unique_idx.values():
                # Get the dihedral potential.
                dihedral = dihedrals0[idx]

                # Get the AtomIdx for the atoms in the dihedral.
                idx0 = info.atomIdx(dihedral.atom0())
                idx1 = info.atomIdx(dihedral.atom1())
                idx2 = info.atomIdx(dihedral.atom2())
                idx3 = info.atomIdx(dihedral.atom3())

                # Cast the function as an AmberDihedral.
                amber_dihedral = _SireMM.AmberDihedral(dihedral.function(), _SireCAS.Symbol("phi"))

                # Start dihedral record.
                file.write("    dihedral\n")

                # Dihedral data.
                file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                file.write("        atom2          %s\n" % mol.atom(idx2).name().value())
                file.write("        atom3          %s\n" % mol.atom(idx3).name().value())
                file.write("        initial_form  ")
                for term in amber_dihedral.terms():
                    file.write(" %5.4f %.1f %7.6f" % (term.k(), term.periodicity(), term.phase()))
                file.write("\n")
                file.write("        final form    ")
                for term in amber_dihedral.terms():
                    file.write(" %5.4f %.1f %7.6f" % (0.0, term.periodicity(), term.phase()))
                file.write("\n")

                # End dihedral record.
                file.write("    enddihedral\n")

            # lambda = 1.
            for idx in dihedrals1_unique_idx.values():
                # Get the dihedral potential.
                dihedral = dihedrals1[idx]

                # Get the AtomIdx for the atoms in the dihedral.
                idx0 = info.atomIdx(dihedral.atom0())
                idx1 = info.atomIdx(dihedral.atom1())
                idx2 = info.atomIdx(dihedral.atom2())
                idx3 = info.atomIdx(dihedral.atom3())

                # Cast the function as an AmberDihedral.
                amber_dihedral = _SireMM.AmberDihedral(dihedral.function(), _SireCAS.Symbol("phi"))

                # Start dihedral record.
                file.write("    dihedral\n")

                # Dihedral data.
                file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                file.write("        atom2          %s\n" % mol.atom(idx2).name().value())
                file.write("        atom3          %s\n" % mol.atom(idx3).name().value())
                file.write("        initial_form  ")
                for term in amber_dihedral.terms():
                    file.write(" %5.4f %.1f %7.6f" % (0.0, term.periodicity(), term.phase()))
                file.write("\n")
                file.write("        final_form    ")
                for term in amber_dihedral.terms():
                    file.write(" %5.4f %.1f %7.6f" % (term.k(), term.periodicity(), term.phase()))
                file.write("\n")

                # End dihedral record.
                file.write("    enddihedral\n")

            # Now add records for the shared dihedrals.
            for idx0, idx1 in dihedrals_shared_idx.values():
                # Get the dihedral potentials.
                dihedral0 = dihedrals0[idx0]
                dihedral1 = dihedrals1[idx1]

                # Get the AtomIdx for the atoms in the dihedral.
                idx0 = info.atomIdx(dihedral0.atom0())
                idx1 = info.atomIdx(dihedral0.atom1())
                idx2 = info.atomIdx(dihedral0.atom2())
                idx3 = info.atomIdx(dihedral0.atom3())

                # Check that an atom in the dihedral is perturbed.
                if _has_pert_atom([idx0, idx1, idx2, idx3], pert_idxs):

                    # Cast the functions as AmberDihedrals.
                    amber_dihedral0 = _SireMM.AmberDihedral(dihedral0.function(), _SireCAS.Symbol("phi"))
                    amber_dihedral1 = _SireMM.AmberDihedral(dihedral1.function(), _SireCAS.Symbol("phi"))

                    # Whether to zero the barrier height of the initial state dihedral.
                    zero_k = False

                    # Whether to force writing the dihedral to the perturbation file.
                    force_write = False

                    # Whether any atom in each end state is a dummy.
                    has_dummy_initial = _has_dummy(mol, [idx0, idx1, idx2, idx3])
                    has_dummy_final = _has_dummy(mol, [idx0, idx1, idx2, idx3], True)

                    # Whether all atoms in each state are dummies.
                    all_dummy_initial = all(_is_dummy(mol, [idx0, idx1, idx2, idx3]))
                    all_dummy_final = all(_is_dummy(mol, [idx0, idx1, idx2, idx3], True))

                    # Dummies are present in both end states, use null potentials.
                    if has_dummy_initial and has_dummy_final:
                        amber_dihedral0 = _SireMM.AmberDihedral()
                        amber_dihedral1 = _SireMM.AmberDihedral()

                    # Dummies in the initial state.
                    elif has_dummy_initial:
                        if all_dummy_initial and not zero_dummy_dihedrals:
                            # Use the potential at lambda = 1 and write to the pert file.
                            amber_dihedral0 = amber_dihedral1
                            force_write = True
                        else:
                            zero_k = True

                    # Dummies in the final state.
                    elif has_dummy_final:
                        if all_dummy_final and not zero_dummy_dihedrals:
                            # Use the potential at lambda = 0 and write to the pert file.
                            amber_dihedral1 = amber_dihedral0
                            force_write = True
                        else:
                            zero_k = True

                    # Only write record if the dihedral parameters change.
                    if zero_k or force_write or amber_dihedral0 != amber_dihedral1:

                        # Initialise a null dihedral.
                        null_dihedral = _SireMM.AmberDihedral()

                        # Don't write the dihedral record if both potentials are null.
                        if not (amber_dihedral0 == null_dihedral and amber_dihedral1 == null_dihedral):

                            # Start dihedral record.
                            file.write("    dihedral\n")

                            # Dihedral data.
                            file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                            file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                            file.write("        atom2          %s\n" % mol.atom(idx2).name().value())
                            file.write("        atom3          %s\n" % mol.atom(idx3).name().value())
                            file.write("        initial_form  ")
                            for term in amber_dihedral0.terms():
                                if zero_k and has_dummy_initial:
                                    k = 0.0
                                else:
                                    k = term.k()
                                file.write(" %5.4f %.1f %7.6f" % (k, term.periodicity(), term.phase()))
                            file.write("\n")
                            file.write("        final_form    ")
                            for term in amber_dihedral1.terms():
                                if zero_k and has_dummy_final:
                                    k = 0.0
                                else:
                                    k = term.k()
                                file.write(" %5.4f %.1f %7.6f" % (k, term.periodicity(), term.phase()))
                            file.write("\n")

                            # End dihedral record.
                            file.write("    enddihedral\n")

            # 5) Impropers.

            # Extract the impropers at lambda = 0 and 1.
            impropers0 = mol.property("improper0").potentials()
            impropers1 = mol.property("improper1").potentials()

            # Dictionaries to store the ImproperIDs at lambda = 0 and 1.
            impropers0_idx = {}
            impropers1_idx = {}

            # Loop over all impropers at lambda = 0.
            for idx, improper in enumerate(impropers0):
                # Get the AtomIdx for the atoms in the improper.
                idx0 = info.atomIdx(improper.atom0())
                idx1 = info.atomIdx(improper.atom1())
                idx2 = info.atomIdx(improper.atom2())
                idx3 = info.atomIdx(improper.atom3())

                # Create the ImproperID.
                improper_id = _SireMol.ImproperID(idx0, idx1, idx2, idx3)

                # Add to the list of ids.
                impropers0_idx[improper_id] = idx

            # Loop over all impropers at lambda = 1.
            for idx, improper in enumerate(impropers1):
                # Get the AtomIdx for the atoms in the improper.
                idx0 = info.atomIdx(improper.atom0())
                idx1 = info.atomIdx(improper.atom1())
                idx2 = info.atomIdx(improper.atom2())
                idx3 = info.atomIdx(improper.atom3())

                # Create the ImproperID.
                improper_id = _SireMol.ImproperID(idx0, idx1, idx2, idx3)

                # Add to the list of ids.
                # You cannot mirror an improper!
                impropers1_idx[improper_id] = idx

            # Now work out the ImproperIDs that are unique at lambda = 0 and 1
            # as well as those that are shared. Note that the ordering of
            # impropers is inconsistent between molecular topology formats so
            # we test all permutations of atom ordering to find matches. This
            # is achieved using the ImproperID.equivalent() method.
            impropers0_unique_idx = {}
            impropers1_unique_idx = {}
            impropers_shared_idx = {}

            # lambda = 0.
            for idx0 in impropers0_idx.keys():
                is_shared = False
                for idx1 in impropers1_idx.keys():
                    if idx0.equivalent(idx1):
                        impropers_shared_idx[idx0] = (impropers0_idx[idx0], impropers1_idx[idx1])
                        is_shared = True
                        break
                if not is_shared:
                    impropers0_unique_idx[idx0] = impropers0_idx[idx0]

            # lambda = 1.
            for idx1 in impropers1_idx.keys():
                is_shared = False
                for idx0 in impropers0_idx.keys():
                    if idx1.equivalent(idx0):
                        # Don't store duplicates.
                        if not idx0 in impropers_shared_idx.keys():
                            impropers_shared_idx[idx1] = (impropers0_idx[idx0], impropers1_idx[idx1])
                        is_shared = True
                        break
                if not is_shared:
                    impropers1_unique_idx[idx1] = impropers1_idx[idx1]

            # First create records for the impropers that are unique to lambda = 0 and 1.

            # lambda = 0.
            for idx in impropers0_unique_idx.values():
                # Get the improper potential.
                improper = impropers0[idx]

                # Get the AtomIdx for the atoms in the improper.
                idx0 = info.atomIdx(improper.atom0())
                idx1 = info.atomIdx(improper.atom1())
                idx2 = info.atomIdx(improper.atom2())
                idx3 = info.atomIdx(improper.atom3())

                # Cast the function as an AmberDihedral.
                amber_dihedral = _SireMM.AmberDihedral(improper.function(), _SireCAS.Symbol("phi"))

                # Start improper record.
                file.write("    improper\n")

                # Dihedral data.
                file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                file.write("        atom2          %s\n" % mol.atom(idx2).name().value())
                file.write("        atom3          %s\n" % mol.atom(idx3).name().value())
                file.write("        initial_form  ")
                for term in amber_dihedral.terms():
                    file.write(" %5.4f %.1f %7.6f" % (term.k(), term.periodicity(), term.phase()))
                file.write("\n")
                file.write("        final form    ")
                for term in amber_dihedral.terms():
                    file.write(" %5.4f %.1f %7.6f" % (0.0, term.periodicity(), term.phase()))
                file.write("\n")

                # End improper record.
                file.write("    endimproper\n")

            # lambda = 1.
            for idx in impropers1_unique_idx.values():
                # Get the improper potential.
                improper = impropers1[idx]

                # Get the AtomIdx for the atoms in the dihedral.
                idx0 = info.atomIdx(improper.atom0())
                idx1 = info.atomIdx(improper.atom1())
                idx2 = info.atomIdx(improper.atom2())
                idx3 = info.atomIdx(improper.atom3())

                # Cast the function as an AmberDihedral.
                amber_dihedral = _SireMM.AmberDihedral(improper.function(), _SireCAS.Symbol("phi"))

                # Start improper record.
                file.write("    improper\n")

                # Dihedral data.
                file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                file.write("        atom2          %s\n" % mol.atom(idx2).name().value())
                file.write("        atom3          %s\n" % mol.atom(idx3).name().value())
                file.write("        initial_form  ")
                for term in amber_dihedral0.terms():
                    file.write(" %5.4f %.1f %7.6f" % (0.0, term.periodicity(), term.phase()))
                file.write("\n")
                file.write("        final_form    ")
                for term in amber_dihedral1.terms():
                    file.write(" %5.4f %.1f %7.6f" % (term.k(), term.periodicity(), term.phase()))
                file.write("\n")

                # End improper record.
                file.write("    endimproper\n")

            # Now add records for the shared impropers.
            for idx0, idx1 in impropers_shared_idx.values():
                # Get the improper potentials.
                improper0 = impropers0[idx0]
                improper1 = impropers1[idx1]

                # Get the AtomIdx for the atoms in the improper.
                idx0 = info.atomIdx(improper0.atom0())
                idx1 = info.atomIdx(improper0.atom1())
                idx2 = info.atomIdx(improper0.atom2())
                idx3 = info.atomIdx(improper0.atom3())

                # Check that an atom in the improper is perturbed.
                if _has_pert_atom([idx0, idx1, idx2, idx3], pert_idxs):

                    # Cast the functions as AmberDihedrals.
                    amber_dihedral0 = _SireMM.AmberDihedral(improper0.function(), _SireCAS.Symbol("phi"))
                    amber_dihedral1 = _SireMM.AmberDihedral(improper1.function(), _SireCAS.Symbol("phi"))

                    # Whether to zero the barrier height of the initial/final improper.
                    zero_k = False

                    # Whether to force writing the improper to the perturbation file.
                    force_write = False

                    # Whether any atom in each end state is a dummy.
                    has_dummy_initial = _has_dummy(mol, [idx0, idx1, idx2, idx3])
                    has_dummy_final = _has_dummy(mol, [idx0, idx1, idx2, idx3], True)

                    # Whether all atoms in each state are dummies.
                    all_dummy_initial = all(_is_dummy(mol, [idx0, idx1, idx2, idx3]))
                    all_dummy_final = all(_is_dummy(mol, [idx0, idx1, idx2, idx3], True))

                    # Dummies are present in both end states, use null potentials.
                    if has_dummy_initial and has_dummy_final:
                        amber_dihedral0 = _SireMM.AmberDihedral()
                        amber_dihedral1 = _SireMM.AmberDihedral()

                    # Dummies in the initial state.
                    elif has_dummy_initial:
                        if all_dummy_initial and not zero_dummy_impropers:
                            # Use the potential at lambda = 1 and write to the pert file.
                            amber_dihedral0 = amber_dihedral1
                            force_write = True
                        else:
                            zero_k = True

                    # Dummies in the final state.
                    elif has_dummy_final:
                        if all_dummy_final and not zero_dummy_impropers:
                            # Use the potential at lambda = 0 and write to the pert file.
                            amber_dihedral1 = amber_dihedral0
                            force_write = True
                        else:
                            zero_k = True

                    # Only write record if the improper parameters change.
                    if zero_k or force_write or amber_dihedral0 != amber_dihedral1:

                        # Initialise a null dihedral.
                        null_dihedral = _SireMM.AmberDihedral()

                        # Don't write the improper record if both potentials are null.
                        if not (amber_dihedral0 == null_dihedral and amber_dihedral1 == null_dihedral):

                            # Start improper record.
                            file.write("    improper\n")

                            # Improper data.
                            file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                            file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                            file.write("        atom2          %s\n" % mol.atom(idx2).name().value())
                            file.write("        atom3          %s\n" % mol.atom(idx3).name().value())
                            file.write("        initial_form  ")
                            for term in amber_dihedral0.terms():
                                if zero_k and has_dummy_initial:
                                    k = 0.0
                                else:
                                    k = term.k()
                                file.write(" %5.4f %.1f %7.6f" % (k, term.periodicity(), term.phase()))
                            file.write("\n")
                            file.write("        final_form    ")
                            for term in amber_dihedral1.terms():
                                if zero_k and has_dummy_final:
                                    k = 0.0
                                else:
                                    k = term.k()
                                file.write(" %5.4f %.1f %7.6f" % (k, term.periodicity(), term.phase()))
                            file.write("\n")

                            # End improper record.
                            file.write("    endimproper\n")

            # End molecule record.
            file.write("endmolecule\n")

        # Finally, convert the molecule to the lambda = 0 state.

        # Make the molecule editable.
        mol = mol.edit()

        # Remove the perturbable molecule flag.
        mol = mol.removeProperty("is_perturbable").molecule()

        # Special handling for the mass and element properties. Perturbed atoms
        # take the mass and atomic number from the maximum of both states,
        # not the lambda = 0 state.
        if mol.hasProperty("mass0") and mol.hasProperty("element0"):
            # See if the mass or element properties exists in the user map.
            new_mass_prop = property_map.get("mass", "mass")
            new_element_prop = property_map.get("element", "element")

            for idx in range(0, mol.nAtoms()):
                # Convert to an AtomIdx.
                idx = _SireMol.AtomIdx(idx)

                # Extract the elements of the end states.
                element0 = mol.atom(idx).property("element0")
                element1 = mol.atom(idx).property("element1")

                # The end states are different elements.
                if element0 != element1:
                    # Extract the mass of the end states.
                    mass0 = mol.atom(idx).property("mass0")
                    mass1 = mol.atom(idx).property("mass1")

                    # Choose the heaviest mass.
                    if mass0.value() > mass1.value():
                        mass = mass0
                    else:
                        mass = mass1

                    # Choose the element with the most protons.
                    if element0.nProtons() > element1.nProtons():
                        element = element0
                    else:
                        element = element1

                    # Set the updated properties.
                    mol = mol.atom(idx).setProperty(new_mass_prop, mass).molecule()
                    mol = mol.atom(idx).setProperty(new_element_prop, element).molecule()

                else:
                    # Use the properties at lambda = 0.
                    mass = mol.atom(idx).property("mass0")
                    mol = mol.atom(idx).setProperty(new_mass_prop, mass).molecule()
                    mol = mol.atom(idx).setProperty(new_element_prop, element0).molecule()

            # Delete redundant properties.
            mol = mol.removeProperty("mass0").molecule()
            mol = mol.removeProperty("mass1").molecule()
            mol = mol.removeProperty("element0").molecule()
            mol = mol.removeProperty("element1").molecule()

        # Rename all properties in the molecule: "prop0" --> "prop".
        # Delete all properties named "prop0" and "prop1".
        for prop in mol.propertyKeys():
            if prop[-1] == "0" and prop != "mass0" and prop != "element0":
                # See if this property exists in the user map.
                new_prop = property_map.get(prop[:-1], prop[:-1])

                # Copy the property using the updated name.
                mol = mol.setProperty(new_prop, mol.property(prop)).molecule()

                # Delete redundant properties.
                mol = mol.removeProperty(prop).molecule()
                mol = mol.removeProperty(prop[:-1] + "1").molecule()

        # Return the updated molecule.
        return mol.commit()

    def _toRegularMolecule(self, property_map={}, is_lambda1=False):
        """Internal function to convert a merged molecule to a regular molecule.

           Parameters
           ----------

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }

           is_lambda1 : bool
               Whether to use the molecule at the lambda = 1 end state.
               By default, the state at lambda = 0 is used.

           Returns
           -------

           molecule : BioSimSpace._SireWrappers.Molecule
               The molecule at the chosen end state.
        """

        if type(is_lambda1) is not bool:
            raise TypeError("'is_lambda1' must be of type 'bool'")

        if is_lambda1:
            lam = "1"
        else:
            lam = "0"

        if not self._is_merged:
            return Molecule(self._sire_object)

        # Extract and copy the Sire molecule.
        mol = self._sire_object.__deepcopy__()

        # Make the molecule editable.
        mol = mol.edit()

        # Remove the perturbable molecule flag.
        mol = mol.removeProperty("is_perturbable").molecule()

        # Rename all properties in the molecule for the corrsponding end state,
        # e.g.: "prop0" --> "prop". Then delete all properties named "prop0"
        # and "prop1".
        for prop in mol.propertyKeys():
            if prop[-1] == lam:
                # See if this property exists in the user map.
                new_prop = property_map.get(prop[:-1], prop[:-1])

                # Copy the property using the updated name.
                mol = mol.setProperty(new_prop, mol.property(prop)).molecule()

                # Delete redundant properties.
                mol = mol.removeProperty(prop[:-1] + "0").molecule()
                mol = mol.removeProperty(prop[:-1] + "1").molecule()

        # Return the updated molecule.
        return Molecule(mol.commit())

    def _merge(self, other, mapping, allow_ring_breaking=False,
            allow_ring_size_change=False, property_map0={}, property_map1={}):
        """Merge this molecule with 'other'.

           Parameters
           ----------

           other : BioSimSpace._SireWrappers.Molecule
               The molecule to merge with.

           mapping : dict
               The mapping between matching atom indices in the two molecules.

           allow_ring_breaking : bool
               Whether to allow the opening/closing of rings during a merge.

           allow_ring_size_change : bool
               Whether to allow changes in ring size.

           property_map0 : dict
               A dictionary that maps "properties" in this molecule to their
               user defined values. This allows the user to refer to properties
               with their own naming scheme, e.g. { "charge" : "my-charge" }

           property_map1 : dict
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

        if type(property_map0) is not dict:
            raise TypeError("'property_map0' must be of type 'dict'")

        if type(property_map1) is not dict:
            raise TypeError("'property_map1' must be of type 'dict'")

        if type(allow_ring_breaking) is not bool:
            raise TypeError("'allow_ring_breaking' must be of type 'bool'")

        if type(allow_ring_size_change) is not bool:
            raise TypeError("'allow_ring_size_change' must be of type 'bool'")

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
        molecule0 = mol._sire_object
        molecule1 = other._sire_object

        # Get the atom indices from the mapping.
        idx0 = mapping.keys()
        idx1 = mapping.values()

        # Create the reverse mapping: molecule1 --> molecule0
        inv_mapping = {v: k for k, v in mapping.items()}

        # Invert the user property mappings.
        inv_property_map0 = {v: k for k, v in property_map0.items()}
        inv_property_map1 = {v: k for k, v in property_map1.items()}

        # Make sure that the molecules have a "forcefield" property and that
        # the two force fields are compatible.

        # Get the user name for the "forcefield" property.
        ff0 = inv_property_map0.get("forcefield", "forcefield")
        ff1 = inv_property_map1.get("forcefield", "forcefield")

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
            if prop in inv_property_map0:
                prop = inv_property_map0[prop]
            props0.append(prop)

        # molecule1
        for prop in molecule1.propertyKeys():
            if prop in inv_property_map1:
                prop = inv_property_map1[prop]
            props1.append(prop)

        # Determine the common properties between the two molecules.
        # These are the properties that can be perturbed.
        shared_props = list(set(props0).intersection(props1))
        del props0
        del props1

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
            # Add an "name0" property.
            edit_mol = edit_mol.atom(atom.index()) \
                               .setProperty("name0", atom.name().value()).molecule()

            # Loop over all atom properties.
            for prop in atom.propertyKeys():
                # Get the actual property name.
                name = inv_property_map0.get(prop, prop)

                # This is a perturbable property. Rename to "property0", e.g. "charge0".
                if name in shared_props:
                    name = name + "0"

                # Add the property to the atom in the merged molecule.
                edit_mol = edit_mol.atom(atom.index()).setProperty(name, atom.property(prop)).molecule()

        # Add the atom properties from molecule1.
        for atom in atoms1:
            # Get the atom index in the merged molecule.
            idx = inv_mapping[atom.index()]

            # Add an "name0" property.
            edit_mol = edit_mol.atom(idx).setProperty("name0", atom.name().value()).molecule()

            # Loop over all atom properties.
            for prop in atom.propertyKeys():
                # Get the actual property name.
                name = inv_property_map1.get(prop, prop)

                # Zero the "charge" and "LJ" property for atoms that are unique to molecule1.
                if name == "charge":
                    edit_mol = edit_mol.atom(idx).setProperty("charge0", 0*_SireUnits.e_charge).molecule()
                elif name == "LJ":
                    edit_mol = edit_mol.atom(idx).setProperty("LJ0", _SireMM.LJParameter()).molecule()
                elif name == "ambertype":
                    edit_mol = edit_mol.atom(idx).setProperty("ambertype0", "du").molecule()
                elif name == "element":
                    edit_mol = edit_mol.atom(idx).setProperty("element0", _SireMol.Element(0)).molecule()
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
            prop0 = inv_property_map0.get("bond", "bond")
            prop1 = inv_property_map1.get("bond", "bond")

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
            prop0 = inv_property_map0.get("angle", "angle")
            prop1 = inv_property_map1.get("angle", "angle")

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
            prop0 = inv_property_map0.get("dihedral", "dihedral")
            prop1 = inv_property_map1.get("dihedral", "dihedral")

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
            prop0 = inv_property_map0.get("improper", "improper")
            prop1 = inv_property_map1.get("improper", "improper")

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

            # Add an "name1" property.
            edit_mol = edit_mol.atom(idx).setProperty("name1", atom.name().value()).molecule()

            # Loop over all atom properties.
            for prop in atom.propertyKeys():
                # Get the actual property name.
                name = inv_property_map1.get(prop, prop)

                # This is a perturbable property. Rename to "property1", e.g. "charge1".
                if name in shared_props:
                    name = name + "1"

                # Add the property to the atom in the merged molecule.
                edit_mol = edit_mol.atom(idx).setProperty(name, atom.property(prop)).molecule()

        # Add the properties from atoms unique to molecule0.
        for atom in atoms0:
            # Add an "name1" property.
            edit_mol = edit_mol.atom(atom.index()) \
                               .setProperty("name1", atom.name().value()).molecule()

            # Loop over all atom properties.
            for prop in atom.propertyKeys():
                # Get the actual property name.
                name = inv_property_map0.get(prop, prop)

                # Zero the "charge" and "LJ" property for atoms that are unique to molecule0.
                if name == "charge":
                    edit_mol = edit_mol.atom(atom.index()).setProperty("charge1", 0*_SireUnits.e_charge).molecule()
                elif name == "LJ":
                    edit_mol = edit_mol.atom(atom.index()).setProperty("LJ1", _SireMM.LJParameter()).molecule()
                elif name == "ambertype":
                    edit_mol = edit_mol.atom(atom.index()).setProperty("ambertype1", "du").molecule()
                elif name == "element":
                    edit_mol = edit_mol.atom(atom.index()).setProperty("element1", _SireMol.Element(0)).molecule()
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
            prop0 = inv_property_map0.get("bond", "bond")
            prop1 = inv_property_map1.get("bond", "bond")

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
            prop0 = inv_property_map0.get("angle", "angle")
            prop1 = inv_property_map1.get("angle", "angle")

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
            prop0 = inv_property_map0.get("dihedral", "dihedral")
            prop1 = inv_property_map1.get("dihedral", "dihedral")

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
            prop0 = inv_property_map0.get("improper", "improper")
            prop1 = inv_property_map1.get("improper", "improper")

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
        # and "bond1" properties, unless a ring is broken or changes size.
        if not (allow_ring_breaking or allow_ring_size_change):
            if edit_mol.property("bond0").nFunctions() != edit_mol.property("bond1").nFunctions():
                raise _IncompatibleError("Inconsistent number of bonds in merged molecule! "
                                         "A ring may have broken, or changed size. If you want to "
                                         "allow this perturbation, try using the 'allow_ring_breaking' "
                                         "or 'allow_ring_size_change' options.")

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
            # Convert to an AtomIdx.
            idx = _SireMol.AtomIdx(x)

            for y in range(x+1, molecule0.nAtoms()):
                # Convert to an AtomIdx.
                idy = _SireMol.AtomIdx(y)

                # Was a ring openend/closed?
                is_ring_broken =  _is_ring_broken(c0, conn, idx, idy, idx, idy)

                # A ring was broken and it is not allowed.
                if is_ring_broken and not allow_ring_breaking:
                    raise _IncompatibleError("The merge has opened/closed a ring. To allow this "
                                             "perturbation, set the 'allow_ring_breaking' option "
                                             "to 'True'.")

                # Did a ring change size?
                is_ring_size_change =  _is_ring_size_changed(c0, conn, idx, idy, idx, idy)

                # A ring changed size and it is not allowed.
                if not is_ring_broken and is_ring_size_change and not allow_ring_size_change:
                    raise _IncompatibleError("The merge has changed the size of a ring. To allow this "
                                             "perturbation, set the 'allow_ring_size_change' option "
                                             "to 'True'. Be aware that this perturbation may not work "
                                             "and a transition through an intermediate state may be "
                                             "preferable.")

                # The connectivity has changed.
                if c0.connectionType(idx, idy) != conn.connectionType(idx, idy):

                    # The connectivity changed for an unknown reason.
                    if not (is_ring_broken or is_ring_size_change):
                        raise _IncompatibleError("The merge has changed the molecular connectivity "
                                                "but a ring didn't open/close or change size. "
                                                "Check your atom mapping.")
        # molecule1
        for x in range(0, molecule1.nAtoms()):
            # Convert to an AtomIdx.
            idx = _SireMol.AtomIdx(x)

            # Map the index to its position in the merged molecule.
            idx_map = inv_mapping[idx]

            for y in range(x+1, molecule1.nAtoms()):
                # Convert to an AtomIdx.
                idy = _SireMol.AtomIdx(y)

                # Map the index to its position in the merged molecule.
                idy_map = inv_mapping[idy]

                # Was a ring openend/closed?
                is_ring_broken =  _is_ring_broken(c1, conn, idx, idy, idx_map, idy_map)

                # A ring was broken and it is not allowed.
                if is_ring_broken and not allow_ring_breaking:
                    raise _IncompatibleError("The merge has opened/closed a ring. To allow this "
                                             "perturbation, set the 'allow_ring_breaking' option "
                                             "to 'True'.")

                # Did a ring change size?
                is_ring_size_change =  _is_ring_size_changed(c1, conn, idx, idy, idx_map, idy_map)

                # A ring changed size and it is not allowed.
                if not is_ring_broken and is_ring_size_change and not allow_ring_size_change:
                    raise _IncompatibleError("The merge has changed the size of a ring. To allow this "
                                             "perturbation, set the 'allow_ring_size_change' option "
                                             "to 'True'. Be aware that this perturbation may not work "
                                             "and a transition through an intermediate state may be "
                                             "preferable.")

                # The connectivity has changed.
                if c1.connectionType(idx, idy) != conn.connectionType(idx_map, idy_map):

                    # The connectivity changed for an unknown reason.
                    if not (is_ring_broken or is_ring_size_change):
                        raise _IncompatibleError("The merge has changed the molecular connectivity "
                                                "but a ring didn't open/close or change size. "
                                                "Check your atom mapping.")

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
        prop0 = inv_property_map0.get("intrascale", "intrascale")
        prop1 = inv_property_map1.get("intrascale", "intrascale")

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

        # Store the two molecular components.
        edit_mol.setProperty("molecule0", molecule0)
        edit_mol.setProperty("molecule1", molecule1)

        # Set the "intrascale" properties.
        edit_mol.setProperty("intrascale0", clj_nb_pairs0)
        edit_mol.setProperty("intrascale1", clj_nb_pairs1)

        # Set the "forcefield" properties.
        edit_mol.setProperty("forcefield0", molecule0.property(ff0))
        edit_mol.setProperty("forcefield1", molecule1.property(ff1))

        # Flag that this molecule is perturbable.
        edit_mol.setProperty("is_perturbable", _SireBase.wrap(True))

        # Update the Sire molecule object of the new molecule.
        mol._sire_object = edit_mol.commit()

        # Flag that the molecule has been merged.
        mol._is_merged = True

        # Store the components of the merged molecule.
        mol._molecule0 = Molecule(molecule0)
        mol._molecule1 = Molecule(molecule1)

        # Return the new molecule.
        return mol

def _has_pert_atom(idxs, pert_idxs):
    """Internal function to check whether a potential contains perturbed atoms.

       Parameters
       ----------

       idxs : [AtomIdx]
           A list of atom indices involved in the potential.

       pert_idxs : [AtomIdx]
           A list of atom indices that are perturbed.

       Returns
       -------

       has_pert_atom : bool
           Whether the potential includes a perturbed atom.
    """

    for idx in idxs:
        if idx in pert_idxs:
            return True

    return False

def _has_dummy(mol, idxs, is_lambda1=False):
    """Internal function to check whether any atom is a dummy.

       Parameters
       ----------

       mol : Sire.Mol.Molecule
           The molecule.

       idxs : [AtomIdx]
           A list of atom indices.

       is_lambda1 : bool
           Whether to check the lambda = 1 state.

       Returns
       -------

       has_dummy : bool
           Whether a dummy atom is present.
    """

    # Set the element property associated with the end state.
    if is_lambda1:
        prop = "element1"
    else:
        prop = "element0"

    dummy = _SireMol.Element(0)

    # Check whether an of the atoms is a dummy.
    for idx in idxs:
        if mol.atom(idx).property(prop) == dummy:
            return True

    return False

def _is_dummy(mol, idxs, is_lambda1=False):
    """Internal function to return whether each atom is a dummy.

       Parameters
       ----------

       mol : Sire.Mol.Molecule
           The molecule.

       idxs : [AtomIdx]
           A list of atom indices.

       is_lambda1 : bool
           Whether to check the lambda = 1 state.

       Returns
       -------

       is_dummy : [bool]
           Whether each atom is a dummy.
    """

    # Set the element property associated with the end state.
    if is_lambda1:
        prop = "element1"
    else:
        prop = "element0"

    # Store a dummy element.
    dummy = _SireMol.Element(0)

    # Initialise a list to store the state of each atom.
    is_dummy = []

    # Check whether each of the atoms is a dummy.
    for idx in idxs:
        is_dummy.append(mol.atom(idx).property(prop) == dummy)

    return is_dummy

def _random_suffix(basename, size=4, chars=_string.ascii_uppercase + _string.digits):
    """Internal helper function to generate a random atom name suffix to avoid
       naming clashes.

       Adapted from:
       https://stackoverflow.com/questions/2257441/random-string-generation-with-upper-case-letters-and-digits-in-python

       Parameters
       ----------

       basename : str
           The base string to which a suffix will be appended.

       size : int
           The maximum width of the string, i.e. len(basename + suffix).

       chars : str
           The set of characters to include in the suffix.

       Returns
       -------

       suffix : str
           The randomly generated suffix.
    """

    basename_size = len(basename)
    if basename_size >= size:
        raise ValueError("Cannot generate suffix for basename '%s'. " % basename
                       + "AMBER atom names can only be 4 characters wide.")
    return "".join(_random.choice(chars) for _ in range(size-basename_size))

def _is_ring_broken(conn0, conn1, idx0, idy0, idx1, idy1):
    """Internal function to test whether a perturbation changes the connectivity
       around two atoms such that a ring is broken.

       Parameters
       ----------

       conn0 : Sire.Mol.Connectivity
           The connectivity object for the first end state.

       conn1 : Sire.Mol.Connectivity
           The connectivity object for the second end state.

       idx0 : Sire.Mol.AtomIdx
           The index of the first atom in the first state.

       idy0 : Sire.Mol.AtomIdx
           The index of the second atom in the first state.

       idx1 : Sire.Mol.AtomIdx
           The index of the first atom in the second state.

       idy1 : Sire.Mol.AtomIdx
           The index of the second atom in the second state.
    """

    # Have we opened/closed a ring? This means that both atoms are part of a
    # ring in one end state (either in it, or on it), whereas at least one
    # are the result of changes in ring size, where atoms remain in or on a
    # ring in both end states.

    # Whether each atom is in a ring in both end states.
    in_ring_idx0 = conn0.inRing(idx0)
    in_ring_idy0 = conn0.inRing(idy0)
    in_ring_idx1 = conn1.inRing(idx1)
    in_ring_idy1 = conn1.inRing(idy1)

    # Whether each atom is on a ring in both end states.
    on_ring_idx0 = _onRing(idx0, conn0)
    on_ring_idy0 = _onRing(idy0, conn0)
    on_ring_idx1 = _onRing(idx1, conn1)
    on_ring_idy1 = _onRing(idy1, conn1)

    # Both atoms are in a ring in one end state and at least one isn't in the other.
    if (in_ring_idx0 & in_ring_idy0) ^ (in_ring_idx1 & in_ring_idy1):
        return True

    # Both atoms are on a ring in one end state and at least one isn't in the other.
    if ((on_ring_idx0 & on_ring_idy0 & (conn0.connectionType(idx0, idy0) == 4))
        ^ (on_ring_idx1 & on_ring_idy1 & (conn1.connectionType(idx1, idy1) == 4))):
        return True

    # Both atoms are in or on a ring in one state and at least one isn't in the other.
    if (((in_ring_idx0 | on_ring_idx0) & (in_ring_idy0 | on_ring_idy0) & (conn0.connectionType(idx0, idy0) == 3)) ^
        ((in_ring_idx1 | on_ring_idx1) & (in_ring_idy1 | on_ring_idy1) & (conn1.connectionType(idx1, idy1) == 3))):
        iscn0 = set(conn0.connectionsTo(idx0)).intersection(set(conn0.connectionsTo(idy0)))
        if (len(iscn0) != 1):
            return True
        common_idx = iscn0.pop()
        in_ring_bond0 = (conn0.inRing(idx0, common_idx) | conn0.inRing(idy0, common_idx))
        iscn1 = set(conn1.connectionsTo(idx1)).intersection(set(conn1.connectionsTo(idy1)))
        if (len(iscn1) != 1):
            return True
        common_idx = iscn1.pop()
        in_ring_bond1 = (conn1.inRing(idx1, common_idx) | conn1.inRing(idy1, common_idx))
        if (in_ring_bond0 ^ in_ring_bond1):
            return True

    # If we get this far, then a ring wasn't broken.
    return False

def _is_ring_size_changed(conn0, conn1, idx0, idy0, idx1, idy1, max_ring_size=12):
    """Internal function to test whether a perturbation changes the connectivity
       around two atoms such that a ring changes size.

       Parameters
       ----------

       conn0 : Sire.Mol.Connectivity
           The connectivity object for the first end state.

       conn1 : Sire.Mol.Connectivity
           The connectivity object for the second end state.

       idx0 : Sire.Mol.AtomIdx
           The index of the first atom in the first state.

       idy0 : Sire.Mol.AtomIdx
           The index of the second atom in the first state.

       idx1 : Sire.Mol.AtomIdx
           The index of the first atom in the second state.

       idy1 : Sire.Mol.AtomIdx
           The index of the second atom in the second state.

       max_ring_size : int
           The maximum size of what is considered to be a ring.
    """

    # Have a ring changed size? If so, then the minimum path size between
    # two atoms will have changed.

    # Work out the paths connecting the atoms in the two end states.
    paths0 = conn0.findPaths(idx0, idy0, max_ring_size)
    paths1 = conn1.findPaths(idx1, idy1, max_ring_size)

    # Initalise the ring size in each end state.
    ring0 = None
    ring1 = None

    # Determine the minimum path in the lambda = 0 state.
    if len(paths0) > 1:
        path_lengths0 = []
        for path in paths0:
            path_lengths0.append(len(path))
        ring0 = min(path_lengths0)

    if ring0 is None:
        return False

    # Determine the minimum path in the lambda = 1 state.
    if len(paths1) > 1:
        path_lengths1 = []
        for path in paths1:
            path_lengths1.append(len(path))
        ring1 = min(path_lengths1)

    # Return whether the ring has changed size.
    if ring1:
        return ring0 != ring1
    else:
        return False

def _onRing(idx, conn):
    """Internal function to test whether an atom is adjacent to a ring.

       Parameters
       ----------

       idx : Sire.Mol.AtomIdx
           The index of the atom

       conn : Sire.Mol.Connectivity
           The connectivity object.

       Returns
       -------

       is_on_ring : bool
           Whether the atom is adjacent to a ring.
    """

    # Loop over all atoms connected to this atom.
    for x in conn.connectionsTo(idx):
        # The neighbour is in a ring.
        if conn.inRing(x) and (not conn.inRing(x, idx)):
            return True

    # If we get this far, then the atom is not adjacent to a ring.
    return False

# Import at bottom of module to avoid circular dependency.
from ._atom import Atom as _Atom
from ._molecules import Molecules as _Molecules
from ._residue import Residue as _Residue
from ._search_result import SearchResult as _SearchResult
from ._system import System as _System
