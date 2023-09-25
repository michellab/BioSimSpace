import itertools as _it
import os as _os
import shutil as _shutil
import tempfile

import numpy as _np
import parmed as _pmd
from sire.legacy import IO as _SireIO
from sire.legacy import Mol as _SireMol

from ._merge import _removeDummies
from ..IO import readMolecules as _readMolecules, saveMolecules as _saveMolecules
from .._SireWrappers import Molecule as _Molecule


def _squash(system, explicit_dummies=False):
    """Internal function which converts a merged BioSimSpace system into an AMBER-compatible format, where all perturbed
    molecules are represented sequentially, instead of in a mixed topology, like in GROMACS. In the current
    implementation, all perturbed molecules are moved at the end of the squashed system. For example, if we have an
    input system, containing regular molecules (M) and perturbed molecules (P):

    M0 - M1 - P0 - M2 - P1 - M3

    This function will return the following squashed system:

    M0 - M1 - M2 - M3 - P0_A - PO_B - P1_A - P1_B

    Where A and B denote the dummyless lambda=0 and lambda=1 states. In addition, we also
    return a mapping between the old unperturbed molecule indices and the new ones. This
    mapping can be used during coordinate update. Updating the coordinates of the perturbed
    molecules, however, has to be done manually through the Python layer.

    Parameters
    ----------

    system : BioSimSpace._SireWrappers.System
        The system.

    explicit_dummies : bool
        Whether to keep the dummy atoms explicit at the endstates or remove them.

    Returns
    -------

    system : BioSimSpace._SireWrappers.System
         The output squashed system.

    mapping : dict(Sire.Mol.MolIdx, Sire.Mol.MolIdx)
         The corresponding molecule-to-molecule mapping. Only the non-perturbable
         molecules are contained in this mapping as the perturbable ones do not
         have a one-to-one mapping and cannot be expressed as a dictionary.
    """
    # Create a copy of the original system.
    new_system = system.copy()

    # Get the perturbable molecules and their corresponding indices.
    pertmol_idxs = [
        i
        for i, molecule in enumerate(system.getMolecules())
        if molecule.isPerturbable()
    ]
    pert_mols = system.getPerturbableMolecules()

    # Remove the perturbable molecules from the system.
    new_system.removeMolecules(pert_mols)

    # Add them back at the end of the system. This is generally faster than keeping their order the same.
    new_indices = list(range(system.nMolecules()))
    for pertmol_idx, pert_mol in zip(pertmol_idxs, pert_mols):
        new_indices.remove(pertmol_idx)
        new_system += _squash_molecule(pert_mol, explicit_dummies=explicit_dummies)

    # Create the old molecule index to new molecule index mapping.
    mapping = {
        _SireMol.MolIdx(idx): _SireMol.MolIdx(i) for i, idx in enumerate(new_indices)
    }

    return new_system, mapping


def _squash_molecule(molecule, explicit_dummies=False):
    """This internal function converts a perturbed molecule to a system that is
    recognisable to the AMBER alchemical code. If the molecule contains a single
    residue, then the squashed system is just the two separate pure endstate
    molecules in order. If the molecule contains regular (R) and perturbable (P)
    resides of the form:

    R0 - R1 - P0 - R2 - P1 - R3

    Then a system containing a single molecule will be returned, which is generated
    by ParmEd's tiMerge as follows:

    R0 - R1 - P0_A - R2 - P1_A - R3 - P0_B - P1_B

    Where A and B denote the dummyless lambda=0 and lambda=1 states.

    Parameters
    ----------

    molecule : BioSimSpace._SireWrappers.Molecule
        The input molecule.

    explicit_dummies : bool
        Whether to keep the dummy atoms explicit at the endstates or remove them.

    Returns
    -------

    system : BioSimSpace._SireWrappers.System
         The output squashed system.
    """
    if not molecule.isPerturbable():
        return molecule

    if explicit_dummies:
        # We make sure we use the same coordinates at both endstates.
        c = molecule.copy()._sire_object.cursor()
        c["coordinates1"] = c["coordinates0"]
        molecule = _Molecule(c.commit())

    # Generate a "system" from the molecule at lambda = 0 and another copy at lambda = 1.
    if explicit_dummies:
        mol0 = molecule.copy()._toRegularMolecule(
            is_lambda1=False, convert_amber_dummies=True, generate_intrascale=True
        )
        mol1 = molecule.copy()._toRegularMolecule(
            is_lambda1=True, convert_amber_dummies=True, generate_intrascale=True
        )
    else:
        mol0 = _removeDummies(molecule, False)
        mol1 = _removeDummies(molecule, True)
    system = (mol0 + mol1).toSystem()

    # We only need to call tiMerge for multi-residue molecules
    if molecule.nResidues() == 1:
        return system

    # Perform the multi-residue squashing with ParmEd as it is much easier and faster.
    with tempfile.TemporaryDirectory() as tempdir:
        # Load in ParmEd.
        _saveMolecules(f"{tempdir}/temp", mol0 + mol1, "prm7,rst7")
        _shutil.move(f"{tempdir}/temp.prm7", f"{tempdir}/temp.parm7")
        parm = _pmd.load_file(f"{tempdir}/temp.parm7", xyz=f"{tempdir}/temp.rst7")

        # Determine the molecule masks.
        mol_mask0 = f"@1-{mol0.nAtoms()}"
        mol_mask1 = f"@{mol0.nAtoms() + 1}-{system.nAtoms()}"

        # Determine the residue masks.
        atom0_offset, atom1_offset = 0, mol0.nAtoms()
        res_atoms0, res_atoms1 = [], []
        for res0, res1, res01 in zip(
            mol0.getResidues(), mol1.getResidues(), molecule.getResidues()
        ):
            if _is_perturbed(res01) or molecule.nResidues() == 1:
                res_atoms0 += list(range(atom0_offset, atom0_offset + res0.nAtoms()))
                res_atoms1 += list(range(atom1_offset, atom1_offset + res1.nAtoms()))
            atom0_offset += res0.nAtoms()
            atom1_offset += res1.nAtoms()
        res_mask0 = _amber_mask_from_indices(res_atoms0)
        res_mask1 = _amber_mask_from_indices(res_atoms1)

        # Merge the residues.
        action = _pmd.tools.tiMerge(parm, mol_mask0, mol_mask1, res_mask0, res_mask1)
        action.output = open(_os.devnull, "w")  # Avoid some of the spam
        action.execute()

        # Reload into BioSimSpace.
        # TODO: prm7/rst7 doesn't work for some reason so we need to use gro/top
        parm.save(f"{tempdir}/squashed.gro", overwrite=True)
        parm.save(f"{tempdir}/squashed.top", overwrite=True)
        squashed_mol = _readMolecules(
            [f"{tempdir}/squashed.gro", f"{tempdir}/squashed.top"]
        )

    return squashed_mol


def _unsquash(system, squashed_system, mapping, **kwargs):
    """Internal function which converts an alchemical AMBER system where the perturbed molecules are
    defined sequentially and updates the coordinates and velocities of an input unsquashed system.
    Refer to the _squash() function documentation to see the structure of the squashed system
    relative to the unsquashed one.

    Parameters
    ----------

    system : BioSimSpace._SireWrappers.System
        The regular unsquashed system.

    squashed_system : BioSimSpace._SireWrappers.System
        The corresponding squashed system.

    mapping : dict(Sire.Mol.MolIdx, Sire.Mol.MolIdx)
        The molecule-molecule mapping generated by _squash().

    kwargs : dict
        A dictionary of optional keyword arguments to supply to _unsquash_molecule().

    Returns
    -------
    system : BioSimSpace._SireWrappers.System
         The output unsquashed system.
    """
    # Create a copy of the original new_system.
    new_system = system.copy()

    # Update the unperturbed molecule coordinates in the original new_system using the mapping.
    if mapping:
        new_system._sire_object, _ = _SireIO.updateCoordinatesAndVelocities(
            new_system._sire_object, squashed_system._sire_object, mapping
        )

    # From now on we handle all perturbed molecules.
    pertmol_idxs = [
        i
        for i, molecule in enumerate(new_system.getMolecules())
        if molecule.isPerturbable()
    ]

    # Get the molecule mapping and combine it with the lambda=0 molecule being prioritised
    molecule_mapping0 = _squashed_molecule_mapping(new_system, is_lambda1=False)
    molecule_mapping1 = _squashed_molecule_mapping(new_system, is_lambda1=True)
    molecule_mapping0_rev = {v: k for k, v in molecule_mapping0.items()}
    molecule_mapping1_rev = {v: k for k, v in molecule_mapping1.items()}
    molecule_mapping_rev = {**molecule_mapping1_rev, **molecule_mapping0_rev}
    molecule_mapping_rev = {
        k: v for k, v in molecule_mapping_rev.items() if v in pertmol_idxs
    }

    # Update the perturbed molecule coordinates based on the molecule mapping
    for merged_idx in set(molecule_mapping_rev.values()):
        pertmol = new_system[merged_idx]
        squashed_idx0 = molecule_mapping0[merged_idx]
        squashed_idx1 = molecule_mapping1[merged_idx]

        if squashed_idx0 == squashed_idx1:
            squashed_molecules = squashed_system[squashed_idx0].toSystem()
        else:
            squashed_molecules = (
                squashed_system[squashed_idx0] + squashed_system[squashed_idx1]
            ).toSystem()

        new_pertmol = _unsquash_molecule(pertmol, squashed_molecules, **kwargs)
        new_system.updateMolecule(merged_idx, new_pertmol)

    return new_system


def _unsquash_molecule(molecule, squashed_molecules, explicit_dummies=False):
    """This internal function loads the coordinates and velocities of squashed molecules
    as defined by the _squash_molecule() function into an unsquashed merged molecule.

    Parameters
    ----------

    molecule : BioSimSpace._SireWrappers.Molecule
        The unsquashed merged molecule whose coordinates and velocities are to be updated.

    squashed_molecules : BioSimSpace._SireWrappers.Molecules
        The corresponding squashed molecule(s) whose coordinates are to be used for updating.

    explicit_dummies : bool
        Whether to keep the dummy atoms explicit at the endstates or remove them.

    Returns
    -------

    molecule : BioSimSpace._SireWrappers.Molecule
         The output updated merged molecule.
    """
    # Get the atom mapping and combine it with the lambda=0 molecule being prioritised
    atom_mapping0 = _squashed_atom_mapping(
        molecule, is_lambda1=False, explicit_dummies=explicit_dummies
    )
    atom_mapping1 = _squashed_atom_mapping(
        molecule, is_lambda1=True, explicit_dummies=explicit_dummies
    )
    atom_mapping = {**atom_mapping1, **atom_mapping0}
    update_velocity = squashed_molecules[0]._sire_object.hasProperty("velocity")

    # Even though the two molecules should have the same coordinates, they might be PBC wrapped differently.
    # Here we take the first common core atom and translate the second molecule.
    if len(squashed_molecules) == 2:
        common_atoms = set(atom_mapping0.keys()) & set(atom_mapping1.keys())
        first_common_atom = list(sorted(common_atoms))[0]
        pertatom0 = squashed_molecules.getAtom(atom_mapping0[first_common_atom])
        pertatom1 = squashed_molecules.getAtom(atom_mapping1[first_common_atom])
        pertatom_coords0 = pertatom0._sire_object.property("coordinates")
        pertatom_coords1 = pertatom1._sire_object.property("coordinates")
        translation_vec = pertatom_coords1 - pertatom_coords0

    # Update the coordinates and velocities.
    siremol = molecule.copy()._sire_object.edit()
    for merged_atom_idx, squashed_atom_idx in atom_mapping.items():
        merged_atom = siremol.atom(_SireMol.AtomIdx(merged_atom_idx))
        squashed_atom = squashed_molecules.getAtom(squashed_atom_idx)

        # Update the coordinates.
        coordinates = squashed_atom._sire_object.property("coordinates")

        # Apply the translation if the atom is coming from the second molecule.
        if len(squashed_molecules) == 2 and squashed_atom_idx in atom_mapping1.values():
            coordinates -= translation_vec

        siremol = merged_atom.setProperty("coordinates0", coordinates).molecule()
        siremol = merged_atom.setProperty("coordinates1", coordinates).molecule()

        # Update the velocities.
        if update_velocity:
            velocities = squashed_atom._sire_object.property("velocity")
            siremol = merged_atom.setProperty("velocity0", velocities).molecule()
            siremol = merged_atom.setProperty("velocity1", velocities).molecule()

    return _Molecule(siremol.commit())


def _squashed_molecule_mapping(system, is_lambda1=False):
    """This internal function returns a dictionary whose keys correspond to the molecule
    index of the each molecule in the original merged system, and whose values
    contain the corresponding index of the same molecule at the specified endstate
    in the squashed system.

    Parameters
    ----------

    system : BioSimSpace._SireWrappers.System
        The input merged system.

    is_lambda1 : bool
        Whether to use the lambda=1 endstate.

    Returns
    -------

    mapping : dict(int, int)
        The corresponding molecule mapping.
    """
    # Get the perturbable molecules and their corresponding indices.
    pertmol_idxs = [i for i, molecule in enumerate(system) if molecule.isPerturbable()]

    # Add them back at the end of the system. This is generally faster than keeping their order the same.
    new_indices = list(range(system.nMolecules()))
    for pertmol_idx in pertmol_idxs:
        new_indices.remove(pertmol_idx)

        # Multi-residue molecules are squashed to one molecule with extra residues.
        if system[pertmol_idx].nResidues() > 1:
            new_indices.append(pertmol_idx)
        # Since we have two squashed molecules, we pick the first one at lambda=0 and the second one at lambda = 1.
        elif not is_lambda1:
            new_indices.extend([pertmol_idx, None])
        else:
            new_indices.extend([None, pertmol_idx])

    # Create the old molecule index to new molecule index mapping.
    mapping = {idx: i for i, idx in enumerate(new_indices) if idx is not None}

    return mapping


def _squashed_atom_mapping(system, is_lambda1=False, environment=True, **kwargs):
    """This internal function returns a dictionary whose keys correspond to the atom
    index of the each atom in the original merged system, and whose values
    contain the corresponding index of the same atom at the specified endstate
    in the squashed system.

    Parameters
    ----------

    system : BioSimSpace._SireWrappers.System
        The input merged system.

    is_lambda1 : bool
        Whether to use the lambda=1 endstate.

    environment : bool
        Whether to include all environment atoms (i.e. ones that are not perturbed).

    kwargs :
        Keyword arguments to pass to _squashed_atom_mapping_molecule().

    Returns
    -------

    mapping : dict(int, int)
        The corresponding atom mapping.
    """
    if isinstance(system, _Molecule):
        return _squashed_atom_mapping(
            system.toSystem(), is_lambda1=is_lambda1, environment=environment, **kwargs
        )

    # Both mappings start from 0 and we add all offsets at the end.
    atom_mapping = {}
    atom_idx, squashed_atom_idx, squashed_atom_idx_perturbed = 0, 0, 0
    squashed_offset = sum(x.nAtoms() for x in system if not x.isPerturbable())
    for molecule in system:
        if molecule.isPerturbable():
            residue_atom_mapping, n_squashed_atoms = _squashed_atom_mapping_molecule(
                molecule,
                offset_merged=atom_idx,
                offset_squashed=squashed_offset + squashed_atom_idx_perturbed,
                is_lambda1=is_lambda1,
                environment=environment,
                **kwargs,
            )
            atom_mapping.update(residue_atom_mapping)
            atom_idx += molecule.nAtoms()
            squashed_atom_idx_perturbed += n_squashed_atoms
        elif molecule.isDecoupled():
            residue_atom_mapping, n_squashed_atoms = _squashed_atom_mapping_molecule(
                molecule,
                offset_merged=atom_idx,
                offset_squashed=squashed_atom_idx,
                is_lambda1=is_lambda1,
                environment=environment,
                **kwargs,
            )
            atom_mapping.update(residue_atom_mapping)
            atom_idx += molecule.nAtoms()
            squashed_atom_idx += n_squashed_atoms
        else:
            atom_indices = _np.arange(atom_idx, atom_idx + molecule.nAtoms())
            squashed_atom_indices = _np.arange(
                squashed_atom_idx, squashed_atom_idx + molecule.nAtoms()
            )
            if environment:
                atom_mapping.update(dict(zip(atom_indices, squashed_atom_indices)))
            atom_idx += molecule.nAtoms()
            squashed_atom_idx += molecule.nAtoms()

    # Convert from NumPy integers to Python integers.
    return {int(k): int(v) for k, v in atom_mapping.items()}


def _squashed_atom_mapping_molecule(
    molecule,
    offset_merged=0,
    offset_squashed=0,
    is_lambda1=False,
    environment=True,
    common=True,
    dummies=True,
    explicit_dummies=False,
):
    """This internal function returns a dictionary whose keys correspond to the atom
    index of the each atom in the original merged molecule, and whose values
    contain the corresponding index of the same atom at the specified endstate
    in the squashed molecule at a particular offset.

    Parameters
    ----------

    molecule : BioSimSpace._SireWrappers.Molecule
        The input merged molecule.

    offset_merged : int
        The index at which to start the merged atom numbering.

    offset_squashed : int
        The index at which to start the squashed atom numbering.

    is_lambda1 : bool
        Whether to use the lambda=1 endstate.

    environment : bool
        Whether to include all environment atoms (i.e. ones that are not perturbed).

    common : bool
        Whether to include all common atoms (i.e. ones that are perturbed but are
        not dummies in the endstate of interest).

    dummies : bool
        Whether to include all dummy atoms (i.e. ones that are perturbed and are
        dummies in the endstate of interest).

    explicit_dummies : bool
        Whether to keep the dummy atoms explicit at the endstates or remove them.

    Returns
    -------

    mapping : dict(int, int)
        The corresponding atom mapping.

    n_atoms : int
        The number of squashed atoms that correspond to the squashed molecule.
    """
    if molecule.isDecoupled():
        # Check if the state 0 is coupled
        coupled_at_lambda0 = _check_decouple(molecule)
        if dummies is False:
            return {}, molecule.nAtoms()
        if coupled_at_lambda0 is (not is_lambda1):
            return {
                offset_merged + i: offset_squashed + i for i in range(molecule.nAtoms())
            }, molecule.nAtoms()
        else:
            return {}, molecule.nAtoms()
    elif not molecule.isPerturbable():
        if environment:
            return {
                offset_merged + i: offset_squashed + i for i in range(molecule.nAtoms())
            }, molecule.nAtoms()
        else:
            return {}, molecule.nAtoms()

    # Both mappings start from 0 and we add all offsets at the end.
    mapping, mapping_lambda1 = {}, {}
    atom_idx_merged, atom_idx_squashed, atom_idx_squashed_lambda1 = 0, 0, 0
    for residue in molecule.getResidues():
        if not (_is_perturbed(residue) or molecule.nResidues() == 1):
            # The residue is not perturbed.
            if common:
                mapping.update(
                    {
                        atom_idx_merged + i: atom_idx_squashed + i
                        for i in range(residue.nAtoms())
                    }
                )
            atom_idx_merged += residue.nAtoms()
            atom_idx_squashed += residue.nAtoms()
        else:
            # The residue is perturbed.

            # Determine the dummy and the non-dummy atoms.
            types0 = [
                atom._sire_object.property("ambertype0") for atom in residue.getAtoms()
            ]
            types1 = [
                atom._sire_object.property("ambertype1") for atom in residue.getAtoms()
            ]

            if explicit_dummies:
                # If both endstates are dummies then we treat them as common core atoms
                dummy0 = dummy1 = _np.asarray(
                    [("du" in x) or ("du" in y) for x, y in zip(types0, types1)]
                )
                common0 = common1 = ~dummy0
                in_mol0 = in_mol1 = _np.asarray([True] * residue.nAtoms())
            else:
                in_mol0 = _np.asarray(["du" not in x for x in types0])
                in_mol1 = _np.asarray(["du" not in x for x in types1])
                dummy0 = ~in_mol1
                dummy1 = ~in_mol0
                common0 = _np.logical_and(in_mol0, ~dummy0)
                common1 = _np.logical_and(in_mol1, ~dummy1)

            ndummy0 = residue.nAtoms() - sum(in_mol1)
            ndummy1 = residue.nAtoms() - sum(in_mol0)
            ncommon = residue.nAtoms() - ndummy0 - ndummy1
            natoms0 = ncommon + ndummy0
            natoms1 = ncommon + ndummy1

            # Determine the full mapping indices for the merged and squashed systems.
            if not is_lambda1:
                atom_indices = _np.arange(
                    atom_idx_merged, atom_idx_merged + residue.nAtoms()
                )[in_mol0]
                squashed_atom_indices = _np.arange(
                    atom_idx_squashed, atom_idx_squashed + natoms0
                )
                mapping_to_update = mapping
            else:
                atom_indices = _np.arange(
                    atom_idx_merged, atom_idx_merged + residue.nAtoms()
                )[in_mol1]
                squashed_atom_indices = _np.arange(
                    atom_idx_squashed_lambda1, atom_idx_squashed_lambda1 + natoms1
                )
                mapping_to_update = mapping_lambda1

            # Determine which atoms to return.
            in_mol_mask = in_mol1 if is_lambda1 else in_mol0
            common_mask = common1 if is_lambda1 else common0
            dummy_mask = dummy1 if is_lambda1 else dummy0
            update_mask = _np.asarray([False] * atom_indices.size)

            if common:
                update_mask = _np.logical_or(update_mask, common_mask[in_mol_mask])
            if dummies:
                update_mask = _np.logical_or(update_mask, dummy_mask[in_mol_mask])

            # Finally update the relevant mapping
            mapping_to_update.update(
                dict(zip(atom_indices[update_mask], squashed_atom_indices[update_mask]))
            )

            # Increment the offsets and continue.
            atom_idx_merged += residue.nAtoms()
            atom_idx_squashed += natoms0
            atom_idx_squashed_lambda1 += natoms1

    # Finally add the appropriate offsets
    if explicit_dummies:
        all_ndummy1 = 0
    else:
        all_ndummy1 = sum(
            "du" in x for x in molecule._sire_object.property("ambertype0").toVector()
        )

    offset_squashed_lambda1 = molecule.nAtoms() - all_ndummy1
    res = {
        **{offset_merged + k: offset_squashed + v for k, v in mapping.items()},
        **{
            offset_merged + k: offset_squashed + offset_squashed_lambda1 + v
            for k, v in mapping_lambda1.items()
        },
    }

    return res, atom_idx_squashed + atom_idx_squashed_lambda1


def _is_perturbed(residue):
    """This determines whether a merged residue is actually perturbed. Note that
    it is possible that this function returns false negatives.

    Parameters
    ----------

    residue : BioSimSpace._SireWrappers.Residue
        The input residue.

    Returns
    -------

    res : bool
        Whether the residue is perturbed.
    """
    # If the elements are different, then we are definitely perturbing.
    elem0 = [atom._sire_object.property("element0") for atom in residue.getAtoms()]
    elem1 = [atom._sire_object.property("element1") for atom in residue.getAtoms()]
    return elem0 != elem1


def _check_decouple(mol):
    """Check the decouple molecule and return whether the ligand is coupled in the state
    A.

    Parameters
    ----------

    mol : BioSimSpace._SireWrappers.Molecule
        The decoupled molecule.

    Returns
    -------

    bool
        Whether the molecule is coupled in state A (0).
    """
    charge_0 = mol._sire_object.property("decouple")["charge"][0].value()
    charge_1 = mol._sire_object.property("decouple")["charge"][1].value()
    LJ_0 = mol._sire_object.property("decouple")["LJ"][0].value()
    LJ_1 = mol._sire_object.property("decouple")["LJ"][1].value()
    if charge_0 != LJ_0 or charge_1 != LJ_1:
        raise ValueError(
            f"The LJ and charge needs to be changed at the same time: "
            f"charge_0:{charge_0};charge_1:{charge_1};LJ_0:{LJ_0};LJ_1:{LJ_1}."
        )
    if charge_0 == charge_1:
        raise ValueError(
            "The state A ({charge_0}) has to be different from state B ({charge_1})."
        )
    return bool(charge_0)


def _amber_mask_from_indices(atom_idxs):
    """Internal helper function to create an AMBER mask from a list of atom indices.

    Parameters
    ----------

    atom_idxs : [int]
        A list of atom indices.

    Returns
    -------

    mask : str
        The AMBER mask.
    """
    # AMBER has a restriction on the number of characters in the restraint
    # mask (not documented) so we can't just use comma-separated atom
    # indices. Instead we loop through the indices and use hyphens to
    # separate contiguous blocks of indices, e.g. 1-23,34-47,...

    if atom_idxs:
        # AMBER masks are 1-indexed, while BioSimSpace indices are 0-indexed.
        atom_idxs = [x + 1 for x in sorted(list(set(atom_idxs)))]
        if not all(isinstance(x, int) for x in atom_idxs):
            raise TypeError("'atom_idxs' must be a list of 'int' types.")
        groups = []
        initial_idx = atom_idxs[0]
        for prev_idx, curr_idx in _it.zip_longest(atom_idxs, atom_idxs[1:]):
            if curr_idx != prev_idx + 1 or curr_idx is None:
                if initial_idx == prev_idx:
                    groups += [str(initial_idx)]
                else:
                    groups += [f"{initial_idx}-{prev_idx}"]
                initial_idx = curr_idx
        mask = "@" + ",".join(groups)
    else:
        mask = ""

    return mask
