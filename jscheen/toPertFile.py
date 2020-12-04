#!/bin/python

# tmp for dev: 
from bss_mol_helper_fns import _has_pert_atom, _has_dummy, \
                _is_dummy, _random_suffix, _is_ring_broken, \
                _is_ring_size_changed, _onRing, _random_suffix
import sys

from pytest import approx as _approx

import os.path as _path
import random as _random
import string as _string

from Sire import Base as _SireBase
from Sire import CAS as _SireCAS
from Sire import IO as _SireIO
from Sire import MM as _SireMM
from Sire import Mol as _SireMol
from Sire import System as _SireSystem
from Sire import Units as _SireUnits

from BioSimSpace import _isVerbose
from BioSimSpace._Exceptions import IncompatibleError as _IncompatibleError
from BioSimSpace.Types import Coordinate as _Coordinate
from BioSimSpace.Types import Length as _Length

# from ._sire_wrapper import SireWrapper as _SireWrapper

def _writePertAtoms(mol, file, pert_type, print_all_atoms=True):
    """Write the merged molecule atoms to a perturbation file.

    Parameters
    ----------

    mol : Sire.Mol.Molecule
    The molecule with properties corresponding to the lamda = 0 state.

    file : object
    Python file object to write Pert lines to.

    pert_type : str
    The type of perturbation to perform.

    Returns
    -------

    None 
    """
    # for testing, remove later:
    print_all_atoms=True

    def atom_sorting_criteria(atom):
        LJ0 = atom.property("LJ0");
        LJ1 = atom.property("LJ1");
        return (atom.name().value(),
                atom.property("ambertype0"),
                atom.property("ambertype1"),
                LJ0.sigma().value(),
                LJ1.sigma().value(),
                LJ0.epsilon().value(),
                LJ1.epsilon().value(),
                atom.property("charge0").value(),
                atom.property("charge1").value())

    if pert_type == "standard":
        # Print all atom records.
        if print_all_atoms:
            for atom in sorted(mol.atoms(), key=lambda atom: atom_sorting_criteria(atom)):
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
            for idx in sorted(pert_idxs, key=lambda idx: atom_sorting_criteria(mol.atom(idx))):
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


    else:        
        # Given multistep protocol, assume print all atom records.
        for atom in sorted(mol.atoms(), key=lambda atom: atom_sorting_criteria(atom)):
            # Start atom record.
            file.write("    atom\n")

            # Get the initial/final Lennard-Jones properties.
            LJ0 = atom.property("LJ0");
            LJ1 = atom.property("LJ1");


            
            # set LJ/charge based on requested perturbed term.
            #######################################################################
            if pert_type == "discharge_soft":
                if (atom.property("element0") == _SireMol.Element("X") or 
                    atom.property("element1") == _SireMol.Element("X")):
                    # In this step, only remove charges from soft-core perturbations.
                    LJ0_value = LJ1_value = LJ0.sigma().value(), LJ0.epsilon().value()

                    charge0_value = atom.property("charge0").value()
                    charge1_value = 0.0
                else:
                    # If only hard atoms in perturbation, hold parameters.
                    LJ0_value = LJ1_value = LJ0.sigma().value(), LJ0.epsilon().value()
                    charge0_value = charge1_value = atom.property("charge0").value()

            #######################################################################
            elif pert_type == "vanish_soft":
                if (atom.property("element0") == _SireMol.Element("X") or 
                    atom.property("element1") == _SireMol.Element("X")):
                    # In this step, only remove LJ from soft-core perturbations.
                    LJ0_value = LJ0.sigma().value(), LJ0.epsilon().value()
                    LJ1_value = (0.0, 0.0)

                    # soft discharge was previous step, so assume 0.0.
                    charge0_value = charge1_value = 0.0
                else:
                    # If only hard atoms in perturbation, hold parameters.
                    LJ0_value = LJ1_value = LJ0.sigma().value(), LJ0.epsilon().value()
                    charge0_value = charge1_value = atom.property("charge0").value()

            #######################################################################
            elif pert_type == "change_bonds":
                # this is handled in _writePertBonds(); states should be set to final
                # states of previous step.
                if (atom.property("element0") == _SireMol.Element("X") or 
                    atom.property("element1") == _SireMol.Element("X")):
                    # In previous steps, soft-core transformations were discharged and vanished.
                    LJ0_value = LJ1_value = (0.0, 0.0)
                    charge0_value = charge1_value = 0.0
                else:
                    # If only hard atoms in perturbation, hold parameters as in previous steps.
                    LJ0_value = LJ1_value = LJ0.sigma().value(), LJ0.epsilon().value()
                    charge0_value = charge1_value = atom.property("charge0").value()

            #######################################################################
            elif pert_type == "change_hard":
                if (atom.property("element0") == _SireMol.Element("X") or 
                    atom.property("element1") == _SireMol.Element("X")):
                    # In this step, soft-core perturbations are untouched.
                    LJ0_value = LJ1_value = (0.0, 0.0)
                    charge0_value = charge1_value = 0.0
                else:
                    # If only hard atoms in perturbation, change all parameters.
                    LJ0_value = LJ0.sigma().value(), LJ0.epsilon().value()
                    LJ1_value = LJ1.sigma().value(), LJ1.epsilon().value()
                    charge0_value = atom.property("charge0").value()
                    charge1_value = atom.property("charge1").value()

            #######################################################################
            elif pert_type == "grow_soft":
                if (atom.property("element0") == _SireMol.Element("X") or 
                    atom.property("element1") == _SireMol.Element("X")):
                    # In this step, soft-core perturbations are grown from 0.
                    LJ0_value = (0.0, 0.0)
                    LJ1_value = LJ1.sigma().value(), LJ1.epsilon().value()
                    charge0_value = 0.0
                    charge1_value = 0.0
                else:
                    # If only hard atoms in perturbation, parameters are already changed.
                    LJ0_value = LJ1_value = LJ1.sigma().value(), LJ1.epsilon().value()
                    charge0_value = charge1_value = atom.property("charge1").value()

            #######################################################################
            elif pert_type == "charge_soft":
                if (atom.property("element0") == _SireMol.Element("X") or 
                    atom.property("element1") == _SireMol.Element("X")):
                    # In this step, soft-core perturbations are charged from 0.
                    LJ0_value = LJ1.sigma().value(), LJ1.epsilon().value()
                    LJ1_value = LJ1.sigma().value(), LJ1.epsilon().value()
                    charge0_value = 0.0
                    charge1_value = atom.property("charge1").value()
                else:
                    # If only hard atoms in perturbation, parameters are already changed.
                    LJ0_value = LJ1_value = LJ1.sigma().value(), LJ1.epsilon().value()
                    charge0_value = charge1_value = atom.property("charge1").value()

            # Write atom data.
            file.write("        name           %s\n" % atom.name().value())
            file.write("        initial_type   %s\n" % atom.property("ambertype0"))
            file.write("        final_type     %s\n" % atom.property("ambertype1"))
            file.write("        initial_LJ     %.5f %.5f\n" % (LJ0_value))
            file.write("        final_LJ       %.5f %.5f\n" % (LJ1_value))
            file.write("        initial_charge %.5f\n" % charge0_value)
            file.write("        final_charge   %.5f\n" % charge1_value)

            # End atom record.
            file.write("    endatom\n")




def _toPertFile(molecule, filename="MORPH.pert", zero_dummy_dihedrals=False,
        zero_dummy_impropers=False, print_all_atoms=False, property_map={},
        pert_type=False):
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

       pert_type : str
           String that specifies whether to use the 'standard' or 'multistep'
           approach in perturbing the LJ, bond and angle terms.

       Returns
       -------

       molecule : Sire.Mol.Molecule
           The molecule with properties corresponding to the lamda = 0 state.
    """

    if not molecule._is_merged:
        raise _IncompatibleError("This isn't a merged molecule. Cannot write perturbation file!")

    if not molecule._sire_object.property("forcefield0").isAmberStyle():
        raise _IncompatibleError("Can only write perturbation files for AMBER style force fields.")

    if type(zero_dummy_dihedrals) is not bool:
        raise TypeError("'zero_dummy_dihedrals' must be of type 'bool'")

    if type(zero_dummy_impropers) is not bool:
        raise TypeError("'zero_dummy_impropers' must be of type 'bool'")

    if type(print_all_atoms) is not bool:
        raise TypeError("'print_all_atoms' must be of type 'bool'")

    if type(property_map) is not dict:
        raise TypeError("'property_map' must be of type 'dict'")
    if pert_type:
      if type(pert_type) is not str:
          raise TypeError("'pert_type' must be of type 'str'")

      allowed_pert_types = [  "standard",
                                "discharge_soft",
                                "vanish_soft",
                                "change_bonds",
                                "change_hard",
                                "grow_soft",
                                "charge_soft"]
      if pert_type not in allowed_pert_types:
        raise ValueError("'pert_type' must be any of: "+str(allowed_pert_types), pert_type)
      # change pertfile name:
      filename = "SOMD_INPUTS/"+pert_type+".pert"




    # Extract and copy the Sire molecule.
    mol = molecule._sire_object.__deepcopy__()

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

        _writePertAtoms(mol, file, pert_type=pert_type)


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

        def sort_bonds(bonds, idx):
            # Get the bond potential.
            bond = bonds[idx]

            # Get the AtomIdx for the atoms in the bond.
            idx0 = info.atomIdx(bond.atom0())
            idx1 = info.atomIdx(bond.atom1())

            return (mol.atom(idx0).name().value(),
                    mol.atom(idx1).name().value())

        # lambda = 0.
        for idx in sorted(bonds0_unique_idx.values(),
                          key=lambda idx: sort_bonds(bonds0, idx)):
            # Get the bond potential.
            bond = bonds0[idx]

            # Get the AtomIdx for the atoms in the bond.
            idx0 = info.atomIdx(bond.atom0())
            idx1 = info.atomIdx(bond.atom1())

            # Cast the function as an AmberBond.
            amber_bond = _SireMM.AmberBond(bond.function(), _SireCAS.Symbol("r"))

            # Start bond record.
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
        for idx in sorted(bonds1_unique_idx.values(),
                          key=lambda idx: sort_bonds(bonds1, idx)):
            # Get the bond potential.
            bond = bonds1[idx]

            # Get the AtomIdx for the atoms in the bond.
            idx0 = info.atomIdx(bond.atom0())
            idx1 = info.atomIdx(bond.atom1())

            # Cast the function as an AmberBond.
            amber_bond = _SireMM.AmberBond(bond.function(), _SireCAS.Symbol("r"))

            # Start bond record.
            file.write("    bond\n")
            if pert_type in [
                "discharge_soft",
                "vanish_soft"
                            ]:
                # Bond data is unchanged.
                file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                file.write("        initial_force  %.5f\n" % 0.0)
                file.write("        initial_equil  %.5f\n" % amber_bond.r0())
                file.write("        final_force    %.5f\n" % 0.0)
                file.write("        final_equil    %.5f\n" % amber_bond.r0())

            elif pert_type == "change_bonds" or pert_type == "standard":
                # bonds are perturbed.
                file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                file.write("        initial_force  %.5f\n" % 0.0)
                file.write("        initial_equil  %.5f\n" % amber_bond.r0())
                file.write("        final_force    %.5f\n" % amber_bond.k())
                file.write("        final_equil    %.5f\n" % amber_bond.r0())                

            elif pert_type in [
                "change_hard",
                "grow_soft",
                "charge_soft"
                            ]:
                # Bond data has already been changed, assume endpoints.
                file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                file.write("        initial_force  %.5f\n" % amber_bond.k())
                file.write("        initial_equil  %.5f\n" % amber_bond.r0())
                file.write("        final_force    %.5f\n" % amber_bond.k())
                file.write("        final_equil    %.5f\n" % amber_bond.r0())





            # End bond record.
            file.write("    endbond\n")

        # Now add records for the shared bonds.
        for idx0, idx1 in sorted(bonds_shared_idx.values(),
                                 key=lambda idx_pair: sort_bonds(bonds0, idx_pair[0])):
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
                    if pert_type in [
                        "discharge_soft",
                        "vanish_soft"
                                    ]:
                        # Bonds are not perturbed.
                        file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                        file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                        file.write("        initial_force  %.5f\n" % amber_bond0.k())
                        file.write("        initial_equil  %.5f\n" % amber_bond0.r0())
                        file.write("        final_force    %.5f\n" % amber_bond0.k())
                        file.write("        final_equil    %.5f\n" % amber_bond0.r0())
                    elif pert_type == "change_bonds" or pert_type == "standard":
                        # Bonds are perturbed.
                        file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                        file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                        file.write("        initial_force  %.5f\n" % amber_bond0.k())
                        file.write("        initial_equil  %.5f\n" % amber_bond0.r0())
                        file.write("        final_force    %.5f\n" % amber_bond1.k())
                        file.write("        final_equil    %.5f\n" % amber_bond1.r0())
                    elif pert_type in [
                        "change_hard",
                        "grow_soft",
                        "charge_soft"
                                    ]:
                        # Bonds are already perturbed.
                        file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                        file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                        file.write("        initial_force  %.5f\n" % amber_bond1.k())
                        file.write("        initial_equil  %.5f\n" % amber_bond1.r0())
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

        def sort_angles(angles, idx):
            # Get the angle potential.
            angle = angles[idx]

            # Get the AtomIdx for the atoms in the angle.
            idx0 = info.atomIdx(angle.atom0())
            idx1 = info.atomIdx(angle.atom1())
            idx2 = info.atomIdx(angle.atom2())

            return (mol.atom(idx1).name().value(),
                    mol.atom(idx0).name().value(),
                    mol.atom(idx2).name().value())

        # lambda = 0.
        for idx in sorted(angles0_unique_idx.values(),
                          key=lambda idx: sort_angles(angles0, idx)):
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

            if pert_type in ["standard", "change_bonds"]:
                # Angle data.
                file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                file.write("        atom2          %s\n" % mol.atom(idx2).name().value())
                file.write("        initial_force  %.5f\n" % amber_angle.k())
                file.write("        initial_equil  %.5f\n" % amber_angle.theta0())
                file.write("        final_force    %.5f\n" % 0.0)
                file.write("        final_equil    %.5f\n" % amber_angle.theta0())
            elif pert_type in [
                "discharge_soft",
                "vanish_soft",
                "change_bonds"]:
                # Angle data, unperturbed.
                file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                file.write("        atom2          %s\n" % mol.atom(idx2).name().value())
                file.write("        initial_force  %.5f\n" % amber_angle.k())
                file.write("        initial_equil  %.5f\n" % amber_angle.theta0())
                file.write("        final_force    %.5f\n" % amber_angle.k())
                file.write("        final_equil    %.5f\n" % amber_angle.theta0())                
            elif pert_type in [
                "change_hard",
                "grow_soft",
                "charge_soft"]:
                # Angle data, already perturbed.
                file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                file.write("        atom2          %s\n" % mol.atom(idx2).name().value())
                file.write("        initial_force  %.5f\n" % 0.0)
                file.write("        initial_equil  %.5f\n" % amber_angle.theta0())
                file.write("        final_force    %.5f\n" % 0.0)
                file.write("        final_equil    %.5f\n" % amber_angle.theta0()) 
            # End angle record.
            file.write("    endangle\n")

        # lambda = 1.
        for idx in sorted(angles1_unique_idx.values(),
                          key=lambda idx: sort_angles(angles1, idx)):
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

            if pert_type in ["standard", "change_bonds"]:
                # Angle data.
                file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                file.write("        atom2          %s\n" % mol.atom(idx2).name().value())
                file.write("        initial_force  %.5f\n" % 0.0)
                file.write("        initial_equil  %.5f\n" % amber_angle.theta0())
                file.write("        final_force    %.5f\n" % amber_angle.k())
                file.write("        final_equil    %.5f\n" % amber_angle.theta0())
            elif pert_type in [
                "discharge_soft",
                "vanish_soft",
                "change_bonds"]:
                # Angle data, unperturbed.
                file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                file.write("        atom2          %s\n" % mol.atom(idx2).name().value())
                file.write("        initial_force  %.5f\n" % 0.0)
                file.write("        initial_equil  %.5f\n" % amber_angle.theta0())
                file.write("        final_force    %.5f\n" % 0.0)
                file.write("        final_equil    %.5f\n" % amber_angle.theta0())                
            elif pert_type in [
                "change_hard",
                "grow_soft",
                "charge_soft"]:
                # Angle data, already perturbed.
                file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                file.write("        atom2          %s\n" % mol.atom(idx2).name().value())
                file.write("        initial_force  %.5f\n" % amber_angle.k())
                file.write("        initial_equil  %.5f\n" % amber_angle.theta0())
                file.write("        final_force    %.5f\n" % amber_angle.k())
                file.write("        final_equil    %.5f\n" % amber_angle.theta0()) 

            # End angle record.
            file.write("    endangle\n")

        # Now add records for the shared angles.
        for idx0, idx1 in sorted(angles_shared_idx.values(),
                                 key=lambda idx_pair: sort_angles(angles0, idx_pair[0])):
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
            if pert_type in ["standard", "change_bonds"]:
                # Angle data.
                file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                file.write("        atom2          %s\n" % mol.atom(idx2).name().value())
                file.write("        initial_force  %.5f\n" % amber_angle0.k())
                file.write("        initial_equil  %.5f\n" % amber_angle0.theta0())
                file.write("        final_force    %.5f\n" % amber_angle1.k())
                file.write("        final_equil    %.5f\n" % amber_angle1.theta0())
            elif pert_type in [
                "discharge_soft",
                "vanish_soft",
                "change_bonds"]:
                # Angle data, unperturbed.
                file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                file.write("        atom2          %s\n" % mol.atom(idx2).name().value())
                file.write("        initial_force  %.5f\n" % amber_angle0.k())
                file.write("        initial_equil  %.5f\n" % amber_angle0.theta0())
                file.write("        final_force    %.5f\n" % amber_angle0.k())
                file.write("        final_equil    %.5f\n" % amber_angle0.theta0())             
            elif pert_type in [
                "change_hard",
                "grow_soft",
                "charge_soft"]:
                # Angle data, already perturbed.
                file.write("        atom0          %s\n" % mol.atom(idx0).name().value())
                file.write("        atom1          %s\n" % mol.atom(idx1).name().value())
                file.write("        atom2          %s\n" % mol.atom(idx2).name().value())
                file.write("        initial_force  %.5f\n" % amber_angle1.k())
                file.write("        initial_equil  %.5f\n" % amber_angle1.theta0())
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

        def sort_dihedrals(dihedrals, idx):
            # Get the dihedral potential.
            dihedral = dihedrals[idx]

            # Get the AtomIdx for the atoms in the angle.
            idx0 = info.atomIdx(dihedral.atom0())
            idx1 = info.atomIdx(dihedral.atom1())
            idx2 = info.atomIdx(dihedral.atom2())
            idx3 = info.atomIdx(dihedral.atom3())

            return (mol.atom(idx1).name().value(),
                    mol.atom(idx2).name().value(),
                    mol.atom(idx0).name().value(),
                    mol.atom(idx3).name().value())

        # lambda = 0.
        for idx in sorted(dihedrals0_unique_idx.values(),
                          key=lambda idx: sort_dihedrals(dihedrals0, idx)):
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
            amber_dihedral_terms_sorted = sorted(
                amber_dihedral.terms(), key=lambda t: (t.k(), t.periodicity(), t.phase()))
            for term in amber_dihedral_terms_sorted:
                if pert_type in [
                    "discharge_soft",
                    "vanish_soft",
                    "change_bonds",
                    "standard"]:
                    file.write(" %5.4f %.1f %7.6f" % (term.k(), term.periodicity(), term.phase()))
                else: 
                    file.write(" %5.4f %.1f %7.6f" % (0.0, term.periodicity(), term.phase()))
            file.write("\n")
            file.write("        final form    ")
            for term in amber_dihedral_terms_sorted:
                if pert_type in [
                    "discharge_soft",
                    "vanish_soft"]:
                    file.write(" %5.4f %.1f %7.6f" % (term.k(), term.periodicity(), term.phase()))
                else:
                    file.write(" %5.4f %.1f %7.6f" % (0.0, term.periodicity(), term.phase()))
            file.write("\n")

            # End dihedral record.
            file.write("    enddihedral\n")

        # lambda = 1.
        for idx in sorted(dihedrals1_unique_idx.values(),
                          key=lambda idx: sort_dihedrals(dihedrals1, idx)):
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
            amber_dihedral_terms_sorted = sorted(
                amber_dihedral.terms(), key=lambda t: (t.k(), t.periodicity(), t.phase()))
            for term in amber_dihedral_terms_sorted:
                if pert_type in [
                    "discharge_soft",
                    "vanish_soft",
                    "change_bonds",
                    "standard"]:
                    file.write(" %5.4f %.1f %7.6f" % (0.0, term.periodicity(), term.phase()))
                else: 
                    file.write(" %5.4f %.1f %7.6f" % (term.k(), term.periodicity(), term.phase()))

            file.write("\n")
            file.write("        final_form    ")
            for term in amber_dihedral_terms_sorted:
                if pert_type in [
                    "discharge_soft",
                    "vanish_soft"]:
                    file.write(" %5.4f %.1f %7.6f" % (0.0, term.periodicity(), term.phase()))
                else:
                    file.write(" %5.4f %.1f %7.6f" % (term.k(), term.periodicity(), term.phase()))
            file.write("\n")

            # End dihedral record.
            file.write("    enddihedral\n")

        # Now add records for the shared dihedrals.
        for idx0, idx1 in sorted(dihedrals_shared_idx.values(),
                                 key=lambda idx_pair: sort_dihedrals(dihedrals0, idx_pair[0])):
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
                        for term in sorted(amber_dihedral0.terms(),
                                           key=lambda t: (t.k(), t.periodicity(), t.phase())):
                            if zero_k and has_dummy_initial:
                                k = 0.0
                            else:
                                k = term.k()
                            file.write(" %5.4f %.1f %7.6f" % (k, term.periodicity(), term.phase()))
                        file.write("\n")
                        file.write("        final_form    ")
                        for term in sorted(amber_dihedral1.terms(),
                                           key=lambda t: (t.k(), t.periodicity(), t.phase())):
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
            amber_dihedral_terms_sorted = sorted(
                amber_dihedral.terms(), key=lambda t: (t.k(), t.periodicity(), t.phase()))
            for term in amber_dihedral_terms_sorted:
                if pert_type in [
                    "discharge_soft",
                    "vanish_soft",
                    "change_bonds",
                    "standard"]:
                    file.write(" %5.4f %.1f %7.6f" % (term.k(), term.periodicity(), term.phase()))
                else:
                    file.write(" %5.4f %.1f %7.6f" % (0.0, term.periodicity(), term.phase()))
                
            file.write("\n")
            file.write("        final form    ")
            for term in amber_dihedral_terms_sorted:
                if pert_type in [
                    "discharge_soft",
                    "vanish_soft"]:
                    file.write(" %5.4f %.1f %7.6f" % (term.k(), term.periodicity(), term.phase()))
                else:
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
            for term in sorted(amber_dihedral0.terms(),
                               key=lambda t: (t.k(), t.periodicity(), t.phase())):
                if pert_type in [
                    "discharge_soft",
                    "vanish_soft",
                    "change_bonds",
                    "standard"]:
                    file.write(" %5.4f %.1f %7.6f" % (0.0, term.periodicity(), term.phase()))
                else:
                    file.write(" %5.4f %.1f %7.6f" % (term.k(), term.periodicity(), term.phase()))
            file.write("\n")
            file.write("        final_form    ")
            for term in sorted(amber_dihedral1.terms(),
                               key=lambda t: (t.k(), t.periodicity(), t.phase())):
                if pert_type in [
                    "discharge_soft",
                    "vanish_soft"]:
                    file.write(" %5.4f %.1f %7.6f" % (0.0, term.periodicity(), term.phase()))
                else:
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
                        for term in sorted(amber_dihedral0.terms(),
                                           key=lambda t: (t.k(), t.periodicity(), t.phase())):
                            if zero_k and has_dummy_initial:
                                k = 0.0
                            else:
                                k = term.k()
                            file.write(" %5.4f %.1f %7.6f" % (k, term.periodicity(), term.phase()))
                        file.write("\n")
                        file.write("        final_form    ")
                        for term in sorted(amber_dihedral1.terms(),
                                           key=lambda t: (t.k(), t.periodicity(), t.phase())):
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
