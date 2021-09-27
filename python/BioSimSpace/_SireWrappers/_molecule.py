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
A thin wrapper around Sire.Mol.Molecule. This is an internal package and should
not be directly exposed to the user.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Molecule"]

from pytest import approx as _approx

import os.path as _path

from Sire import Base as _SireBase
from Sire import IO as _SireIO
from Sire import MM as _SireMM
from Sire import Mol as _SireMol
from Sire import System as _SireSystem
from Sire import Units as _SireUnits

from BioSimSpace import _isVerbose
from BioSimSpace._Exceptions import IncompatibleError as _IncompatibleError
from BioSimSpace.Types import Coordinate as _Coordinate
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

        # Set the water model used to parameterise any structural ions
        # in the molecule.
        self._ion_water_model = None

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
            self._is_perturbable = molecule._is_perturbable

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
            for molecule in other:
                molecules.append(molecule)

        # Unsupported.
        else:
            raise TypeError("'other' must be of type 'BioSimSpace._SireWrappers.System', "
                            "'BioSimSpace._SireWrappers.Molecule', 'BioSimSpace._SireWrappers.Molecules' "
                            "or a list of 'BioSimSpace._SireWrappers.Molecule' types")

        # Add this molecule to the container and return.
        if is_sire_container:
            return other + self

        # Create a new Molecules container.
        else:
            return _Molecules(molecules)

    def copy(self):
        """Return a copy of this Molecule.

           Returns
           -------

           molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
               A copy of the object.
        """
        # Copy the Sire object.
        mol = self._sire_object.__deepcopy__()

        # Give the molecule a unique number.
        mol = mol.edit() \
                 .renumber(_SireMol.MolNum.getUniqueNumber()) \
                 .commit().molecule()

        return Molecule(mol)

    def number(self):
        """Return the number of the molecule. Each molecule has a unique
           identification number.

           Returns
           -------

           mol_num : int
               The unique number of the molecule.
        """
        return self._sire_object.number().value()

    def coordinates(self, property_map={}):
        """Return the coordinates of the atoms in the molecule.

           Parameters
           ----------

           property_map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }

           Returns
           -------

           [coordinates] : [class:`Coordinate <BioSimSpace.Types.Coordinate>`]
               The coordinates of the atoms in the molecule.
        """
        prop = property_map.get("coordinates", "coordinates")

        # Get the "coordinates" property from the molecule.
        try:
            sire_coord = self._sire_object.property(prop).toVector()
            coordinates = []
            for coord in sire_coord:
                coordinates.append(_Coordinate(_Length(coord[0], "Angstrom"),
                                               _Length(coord[1], "Angstrom"),
                                               _Length(coord[2], "Angstrom")))
        except:
            return None

        # Return the coordinates.
        return coordinates

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

    def isPerturbable(self):
        """Whether this molecule is perturbable, i.e. it can be used in a
           free-energy perturbation simulation.

           Returns
           -------

           is_perturbable : bool
               Whether the molecule is perturbable.
        """
        return self._is_perturbable

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

    def isAmberWater(self):
        """Whether this is an AMBER format water molecule.

           Returns
           -------

           is_amber_water : bool
               Whether this molecule is an AMBER format water.
        """

        # First check that this is a water molecule.
        if not self.isWater():
            return False

        # Now check the residue name.
        if self.getResidues()[0].name() != "WAT":
            return False

        # Now check the atom names.

        # Get the number of atoms.
        num_atoms = self.nAtoms()

        # SPC/E or TIP3P.
        if num_atoms == 3:
            atom_names = ["O", "H1", "H2"]
            for atom in self.getAtoms():
                if atom.name() not in atom_names:
                    return False

        # TIP4P.
        elif num_atoms == 4:
            atom_names = ["O", "H1", "H2", "EPW"]
            for atom in self.getAtoms():
                if atom.name() not in atom_names:
                    return False

        # TIP5P.
        elif num_atoms == 5:
            atom_names = ["O", "H1", "H2", "EP1", "EP2"]
            for atom in self.getAtoms():
                if atom.name() not in atom_names:
                    return False

        # If we get this far, then it is an AMBER format water.
        return True

    def isGromacsWater(self):
        """Whether this is a GROMACS format water molecule.

           Returns
           -------

           is_gromacs_water : bool
               Whether this molecule is a GROMACS format water.
        """

        # First check that this is a water molecule.
        if not self.isWater():
            return False

        # Now check the residue name.
        if self.getResidues()[0].name() != "SOL":
            return False

        # Get the number of atoms.
        num_atoms = self.nAtoms()

        # SPC/E or TIP3P.
        if num_atoms == 3:
            atom_names = ["OW", "HW1", "HW2"]
            for atom in self.getAtoms():
                if atom.name() not in atom_names:
                    return False

        # TIP4P.
        elif num_atoms == 4:
            atom_names = ["OW", "HW1", "HW2", "MW"]
            for atom in self.getAtoms():
                if atom.name() not in atom_names:
                    return False

        # TIP5P.
        elif num_atoms == 5:
            atom_names = ["OW", "HW1", "LP1", "LP2"]
            for atom in self.getAtoms():
                if atom.name() not in atom_names:
                    return False

        # If we get this far, then it is a GROMACS format water.
        return True

    def toSystem(self):
        """Convert a single Molecule to a System.

           Returns
           -------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
        """
        return _System(self)

    def search(self, query, property_map={}):
        """Search the molecule for atoms and residues. Search results will be
           reduced to their minimal representation, i.e. a residue containing
           a single atom will be returned as a atom.

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

        if type(property_map) is not dict:
            raise TypeError("'property_map' must be of type 'dict'")

        # Initialise a list to hold the search results.
        results = []

        try:
            # Query the Sire molecule.
            search_result = _SireMol.Select(query)(self._sire_object, property_map)

        except Exception as e:
            msg = "'Invalid search query: %r" % query
            if _isVerbose():
                raise ValueError(msg) from e
            else:
                raise ValueError(msg) from None

        return _SearchResult(search_result)

    def makeCompatibleWith(self, molecule, property_map={}, overwrite=True,
            rename_atoms=False, verbose=False):
        """Make this molecule compatible with passed one, i.e. match atoms and
           add all additional properties from the passed molecule while preserving
           the topology and naming/numbering convention of this molecule.

           Parameters
           ----------

           molecule : Sire.Mol.Molecule, :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                      Sire.System.System, :class:`System <Sire._SireWrappers.System>`
               The molecule, or system of molecules, to match with.

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

        is_system = False
        if isinstance(molecule, _SireMol.Molecule):
            mol1 = molecule
        elif type(molecule) is Molecule:
            mol1 = molecule._sire_object
        elif isinstance(molecule, _SireSystem.System):
            mol1 = molecule
            is_system = True
        elif type(molecule) is _System:
            mol1 = molecule._sire_object
            is_system = True
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
        num_atoms0 = mol0.nAtoms()

        # Work out the number of atoms in mol1.
        if is_system:
            num_atoms1 = _System(mol1).nAtoms()
        else:
            num_atoms1 = mol1.nAtoms()

        # The new molecule must have the same number of atoms.
        if num_atoms1 != num_atoms0:
            if is_system:
                raise _IncompatibleError("The passed system is incompatible with the original! "
                                         "self.nAtoms() = %d, other.nAtoms() = %d" % (num_atoms0, num_atoms1))
            else:
                raise _IncompatibleError("The passed molecule is incompatible with the original! "
                                         "self.nAtoms() = %d, other.nAtoms() = %d" % (num_atoms0, num_atoms1))

        # Whether the atoms have been renamed.
        is_renamed = False

        if not is_system:
            # Instantiate the default atom matcher (match by residue index and atom name).
            matcher = _SireMol.ResIdxAtomNameMatcher()

            # Match the atoms based on residue index and atom name.
            matches = matcher.match(mol0, mol1)

            # Have we matched all of the atoms?
            if len(matches) < num_atoms0:
                # Atom names might have changed. Try to match by residue index
                # and coordinates.
                matcher = _SireMol.ResIdxAtomCoordMatcher()
                matches = matcher.match(mol0, mol1)

                # We need to rename the atoms.
                is_renamed = True

                # Have we matched all of the atoms?
                if len(matches) < num_atoms0:
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

                # Create a dictionary mapping the atom matches from mol1 to mol0.
                inv_matches = {}
                for key, value in matches.items():
                    inv_matches[value] = key

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
                                        propty = propty.makeCompatibleWith(mol0, inv_matches)
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

        # Atoms from molecule in the passed system (mol1) need to be matched against atoms
        # from the corresponding residue in mol0 and the properties aggregated.
        else:
            # Initialise a list to hold the matches for each molecule in mol1.
            matches = []

            # Tally counter for the total number of matches.
            num_matches = 0

            # Initialise the offset.
            offset = 0

            # Get the molecule numbers in the system.
            mol_nums = mol1.molNums()

            # Loop over all molecules in mol1.
            for num in mol_nums:
                # Extract the numbered molecule from the system mol1.
                mol = mol1[num]

                # Initialise the matcher.
                matcher = _SireMol.ResIdxAtomCoordMatcher(_SireMol.ResIdx(offset))

                # Get the matches for this molecule and append to the list.
                match = matcher.match(mol0, mol)
                matches.append(match)
                num_matches += len(match)

                # Increment the offset.
                offset += mol.nResidues()

            # Have we matched all of the atoms?
            if num_matches < num_atoms0:
                raise _IncompatibleError("Failed to match all atoms!")

            # Make the molecule editable.
            edit_mol = mol0.edit()

            # Create objects to hold all of the potential terms.
            bonds      = _SireMM.TwoAtomFunctions(edit_mol.info())
            angles     = _SireMM.ThreeAtomFunctions(edit_mol.info())
            dihedrals  = _SireMM.FourAtomFunctions(edit_mol.info())
            impropers  = _SireMM.FourAtomFunctions(edit_mol.info())

            # Next we need to work out what properties can be set at the atom level,
            # and which are associated with the molecule as a whole.
            mol_props = []
            atom_props = []

            # Each atom should have the same set of properties so we can check the
            # first atom in each molecule.
            for num in mol_nums:
                mol = mol1[num]

                # Get the molecule and atom properties.
                props_mol = mol.propertyKeys()
                props_atom = mol.atoms()[0].propertyKeys()

                # Check the atomic properties and add any new ones to the list.
                for prop in props_atom:
                    prop = property_map.get(prop, prop)
                    if prop not in atom_props:
                        atom_props.append(prop)

                # Check the molecular properties and add any new ones to the list.
                for prop in props_mol:
                    prop = property_map.get(prop, prop)
                    if prop not in mol_props and prop not in atom_props:
                        mol_props.append(prop)

            # Create a list of excluded molecular properties. These are ones that
            # must be re-mapped manually and set by hand.
            excluded_props = [property_map.get("bond", "bond"),
                              property_map.get("angle", "angle"),
                              property_map.get("dihedral", "dihedral"),
                              property_map.get("improper", "improper"),
                              property_map.get("connectivity", "connectivity"),
                              property_map.get("intrascale", "intrascale"),
                              property_map.get("parameters", "parameters")]

            # Loop over all atoms within each molecule. We set the allowed atomic
            # properties and build the molecular ones as we go.
            for idx, num in enumerate(mol_nums):
                # Extract the molecule and its associated info object.
                mol = mol1[num]
                info = mol.info()

                # An inverse mapping for indices in mol1 to those in mol0.
                inv_matches = {}

                # Loop over all matching atom pairs for this molecule.
                for idx0, idx1 in matches[idx].items():
                    # Add indices to the inverse mapping.
                    inv_matches[idx1] = idx0

                    # Check the atom names and see if they need updating.
                    if rename_atoms:
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

                    # Loop over all atom properties.
                    for prop in atom_props:
                        # This is a new property, or we're allowed to overwrite.
                        if (not mol0.hasProperty(prop)) or overwrite:
                            if verbose:
                                print("  %-20s %s --> %s" % (prop, idx1, idx0))
                            try:
                                edit_mol = edit_mol.atom(idx0).setProperty(prop, mol.atom(idx1).property(prop)).molecule()
                            except Exception as e:
                                msg = "Failed to copy property '%s' from %s to %s." % (prop, idx1, idx0)
                                if _isVerbose():
                                    raise _IncompatibleError(msg) from e
                                else:
                                    raise _IncompatibleError(msg) from None

                # Now deal with the molecular properties.
                for prop in mol_props:
                    # Get the property name from the user mapping.
                    prop = property_map.get(prop, prop)

                    # This is a new property, or we're allowed to overwrite, and it's not excluded.
                    if ((not mol0.hasProperty(prop)) or overwrite) and prop not in excluded_props:
                        if verbose:
                            print("  %s" % prop)

                        # Get the property from the parameterised molecule.
                        propty = mol.property(prop)

                        # Try making the property compatible with the original molecule.
                        if hasattr(propty, "makeCompatibleWith"):
                            try:
                                propty = propty.makeCompatibleWith(mol0, inv_matches)
                            except Exception as e:
                                msg = "Incompatible property: %s" % prop
                                if _isVerbose():
                                    raise _IncompatibleError(msg) from e
                                else:
                                    raise _IncompatibleError(msg) from None

                        # Now try to set the property.
                        edit_mol.setProperty(prop, propty)

                # Now re-map and build the properties for each of the potential terms.

                # Bonds.
                prop = property_map.get("bond", "bond")
                if mol.hasProperty(prop):
                    if verbose:
                        print("  %s" % prop)
                    for bond in mol.property(prop).potentials():
                        # Extract the bond information.
                        atom0 = info.atomIdx(bond.atom0())
                        atom1 = info.atomIdx(bond.atom1())
                        exprn = bond.function()

                        # Map the atom indices to their position in the merged molecule.
                        atom0 = inv_matches[atom0]
                        atom1 = inv_matches[atom1]

                        # Set the new bond.
                        bonds.set(atom0, atom1, exprn)

                # Angles.
                prop = property_map.get("angle", "angle")
                if mol.hasProperty(prop):
                    if verbose:
                        print("  %s" % prop)
                    for angle in mol.property(prop).potentials():
                        # Extract the angle information.
                        atom0 = info.atomIdx(angle.atom0())
                        atom1 = info.atomIdx(angle.atom1())
                        atom2 = info.atomIdx(angle.atom2())
                        exprn = angle.function()

                        # Map the atom indices to their position in the merged molecule.
                        atom0 = inv_matches[atom0]
                        atom1 = inv_matches[atom1]
                        atom2 = inv_matches[atom2]

                        # Set the new angle.
                        angles.set(atom0, atom1, atom2, exprn)

                # Dihedrals.
                prop = property_map.get("dihedral", "dihedral")
                if mol.hasProperty(prop):
                    if verbose:
                        print("  %s" % prop)
                    for dihedral in mol.property(prop).potentials():
                        # Extract the dihedral information.
                        atom0 = info.atomIdx(dihedral.atom0())
                        atom1 = info.atomIdx(dihedral.atom1())
                        atom2 = info.atomIdx(dihedral.atom2())
                        atom3 = info.atomIdx(dihedral.atom3())
                        exprn = dihedral.function()

                        # Map the atom indices to their position in the merged molecule.
                        atom0 = inv_matches[atom0]
                        atom1 = inv_matches[atom1]
                        atom2 = inv_matches[atom2]
                        atom3 = inv_matches[atom3]

                        # Set the new dihedral.
                        dihedrals.set(atom0, atom1, atom2, atom3, exprn)

                # Impropers.
                prop = property_map.get("improper", "improper")
                if mol.hasProperty(prop):
                    if verbose:
                        print("  %s" % prop)
                    for improper in mol.property(prop).potentials():
                        # Extract the improper information.
                        atom0 = info.atomIdx(improper.atom0())
                        atom1 = info.atomIdx(improper.atom1())
                        atom2 = info.atomIdx(improper.atom2())
                        atom3 = info.atomIdx(improper.atom3())
                        exprn = improper.function()

                        # Map the atom indices to their position in the merged molecule.
                        atom0 = inv_matches[atom0]
                        atom1 = inv_matches[atom1]
                        atom2 = inv_matches[atom2]
                        atom3 = inv_matches[atom3]

                        # Set the new improper.
                        impropers.set(atom0, atom1, atom2, atom3, exprn)

            # Set properties for the molecular potential.

            # Bonds.
            if bonds.nFunctions() > 0:
                prop = property_map.get("bond", "bond")
                if verbose:
                    print("  %s" % prop)
                try:
                    edit_mol.setProperty(prop, bonds)
                except Exception as e:
                    msg = "Incompatible property: %s" % prop
                    if _isVerbose():
                        raise _IncompatibleError(msg) from e
                    else:
                        raise _IncompatibleError(msg) from None

            # Angles.
            if angles.nFunctions() > 0:
                prop = property_map.get("angle", "angle")
                if verbose:
                    print("  %s" % prop)
                try:
                    edit_mol.setProperty(prop, angles)
                except Exception as e:
                    msg = "Incompatible property: %s" % prop
                    if _isVerbose():
                        raise _IncompatibleError(msg) from e
                    else:
                        raise _IncompatibleError(msg) from None

            # Dihedrals.
            if dihedrals.nFunctions() > 0:
                prop = property_map.get("dihedral", "dihedral")
                if verbose:
                    print("  %s" % prop)
                try:
                    edit_mol.setProperty(prop, dihedrals)
                except Exception as e:
                    msg = "Incompatible property: %s" % prop
                    if _isVerbose():
                        raise _IncompatibleError(msg) from e
                    else:
                        raise _IncompatibleError(msg) from None

            # Impropers.
            if impropers.nFunctions() > 0:
                prop = property_map.get("improper", "improper")
                if verbose:
                    print("  %s" % prop)
                try:
                    edit_mol.setProperty(prop, impropers)
                except Exception as e:
                    msg = "Incompatible property: %s" % prop
                    if _isVerbose():
                        raise _IncompatibleError(msg) from e
                    else:
                        raise _IncompatibleError(msg) from None

            # Now generate the molecular connectivity.
            if bonds.nFunctions() > 0:
                prop = property_map.get("connectivity", "connectivity")
                if verbose:
                    print("  %s" % prop)
                conn = _SireMol.Connectivity(edit_mol.info()).edit()

                # Connect the bonded atoms. Connectivity is the same at lambda = 0
                # and lambda = 1.
                for bond in bonds.potentials():
                    conn.connect(bond.atom0(), bond.atom1())
                conn = conn.commit()

                try:
                    edit_mol.setProperty(prop, conn)
                except Exception as e:
                    msg = "Incompatible property: %s" % prop
                    if _isVerbose():
                        raise _IncompatibleError(msg) from e
                    else:
                        raise _IncompatibleError(msg) from None

            # Next we construct the intrascale matrix for the non-bonded
            # interactions. It is prohibitively slow to do this on-the-fly so
            # we use the GroTop parser to re-construct it for us, then copy it
            # back into the original system.

            prop = property_map.get("intrascale", "intrascale")
            if (not mol0.hasProperty(prop)) or overwrite:
                if verbose:
                    print("  %s" % prop)
                # Delete any existing intrascale property from the molecule.
                if mol0.hasProperty(prop):
                    edit_mol.removeProperty(prop)
                mol = edit_mol.commit()
                # Convert to a "GROMACS system" using the GroTop parser.
                gro_system = _SireIO.GroTop(Molecule(mol).toSystem()._sire_object,
                             _SireBase.PropertyMap(property_map)).toSystem()
                # Extract the only molecule in the system.
                gro_mol = gro_system[_SireMol.MolIdx(0)]
                edit_mol = mol.edit()
                try:
                    edit_mol.setProperty(prop, gro_mol.property(prop))
                except Exception as e:
                    msg = "Incompatible property: %s" % prop
                    if _isVerbose():
                        raise _IncompatibleError(msg) from e
                    else:
                        raise _IncompatibleError(msg) from None

            # Finally, commit the changes to the internal object.
            self._sire_object = edit_mol.commit()

    def _getPropertyMap0(self):
        """Generate a property map for the lambda = 0 state of the merged molecule."""

        property_map = {}

        if self._is_perturbable:
            for prop in self._sire_object.propertyKeys():
                if prop[-1] == "0":
                    property_map[prop[:-1]] = prop

        return property_map

    def _getPropertyMap1(self):
        """Generate a property map for the lambda = 1 state of the merged molecule."""

        property_map = {}

        if self._is_perturbable:
            for prop in self._sire_object.propertyKeys():
                if prop[-1] == "1":
                    property_map[prop[:-1]] = prop

        return property_map

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

        # Flag that the molecule is perturbable.
        self._is_perturbable = True

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
        charge = self.charge(property_map=property_map).magnitude()

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

    def _toRegularMolecule(self, property_map={},
            is_lambda1=False, convert_amber_dummies=False):
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

           convert_amber_dummies : bool
               Whether to convert dummies to the correct AMBER formatting for
               non-FEP simulations. This will replace the "du" ambertype
               and "Xx" element with the properties from the other end state.

           Returns
           -------

           molecule : BioSimSpace._SireWrappers.Molecule
               The molecule at the chosen end state.
        """

        if type(is_lambda1) is not bool:
            raise TypeError("'is_lambda1' must be of type 'bool'")

        if type(convert_amber_dummies) is not bool:
            raise TypeError("'convert_amber_dummies' must be of type 'bool'")

        if is_lambda1:
            lam = "1"
        else:
            lam = "0"

        if not self._is_perturbable:
            return Molecule(self._sire_object)

        # Extract and copy the Sire molecule.
        mol = self._sire_object.__deepcopy__()

        # Make the molecule editable.
        mol = mol.edit()

        # Remove the perturbable molecule flag.
        mol = mol.removeProperty("is_perturbable").molecule()

        # Flag that the molecule was perturbable, so that dummies should be
        # treated as "normal" atoms.
        mol = mol.setProperty("was_perturbable", _SireBase.wrap(True)).molecule()

        # Rename all properties in the molecule for the corresponding end state,
        # e.g.: "prop0" --> "prop". Then delete all properties named "prop0"
        # and "prop1".
        for prop in mol.propertyKeys():
            if prop[-1] == lam:
                # See if this property exists in the user map.
                new_prop = property_map.get(prop[:-1], prop[:-1])

                # Copy the property using the updated name.
                mol = mol.setProperty(new_prop, mol.property(prop)).molecule()

                # Store the amber types in the opposie end state.
                if prop[:-1] == "ambertype":
                    if lam == "0":
                        amber_types = mol.property("ambertype1").toVector()
                    else:
                        amber_types = mol.property("ambertype0").toVector()

                elif prop[:-1] == "element":
                    if lam == "0":
                        elements = mol.property("element1").toVector()
                    else:
                        elements = mol.property("element0").toVector()

                else:
                    # Delete redundant properties.
                    mol = mol.removeProperty(prop[:-1] + "1").molecule()
                    mol = mol.removeProperty(prop[:-1] + "0").molecule()

        # Convert ambertype and element property of dummies to those of the
        # other end state.
        if convert_amber_dummies:
            amber_type = property_map.get("ambertype", "ambertype")
            element = property_map.get("element", "element")
            if mol.hasProperty(amber_type) and mol.hasProperty(element):
                # Search for any dummy atoms.
                search = mol.search("element Xx")

                # Replace the ambertype.
                for dummy in search:
                    mol = mol.atom(dummy.index()) \
                             .setProperty(amber_type, amber_types[dummy.index().value()]).molecule()
                    mol = mol.atom(dummy.index()) \
                             .setProperty(element, elements[dummy.index().value()]).molecule()

                # Delete redundant properties.
                mol = mol.removeProperty("ambertype0").molecule()
                mol = mol.removeProperty("ambertype1").molecule()
                mol = mol.removeProperty("element0").molecule()
                mol = mol.removeProperty("element1").molecule()

        # Return the updated molecule.
        return Molecule(mol.commit())

# Import at bottom of module to avoid circular dependency.
from ._atom import Atom as _Atom
from ._molecules import Molecules as _Molecules
from ._residue import Residue as _Residue
from ._search_result import SearchResult as _SearchResult
from ._system import System as _System
