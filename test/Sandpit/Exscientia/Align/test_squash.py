import pickle
import pytest

import BioSimSpace.Sandpit.Exscientia as BSS


@pytest.fixture
def perturbed_system():
    # N_atoms are: 12, 15, 18, 21, 24, 27 and 30.
    mol_smiles = ["c1ccccc1", "c1ccccc1C", "c1ccccc1CC", "c1ccccc1CCC", "c1ccccc1CCCC", "c1ccccc1CCCCC",
                  "c1ccccc1CCCCCC"]
    mols = [BSS.Parameters.gaff(smi).getMolecule() for smi in mol_smiles]
    pert_mols = [mols[0], BSS.Align.merge(mols[1], mols[2]), mols[3], mols[4], BSS.Align.merge(mols[5], mols[6])]
    system = BSS._SireWrappers.System(pert_mols)
    return system


@pytest.fixture
def perturbed_tripeptide():
    return pickle.load(open("test/Sandpit/Exscientia/input/morphs/merged_tripeptide.pickle", "rb"))


def test_squash(perturbed_system):
    squashed_system, mapping = BSS.Align._squash._squash(perturbed_system)
    assert len(squashed_system) == 7
    n_atoms = [mol.nAtoms() for mol in squashed_system]
    assert squashed_system[-2].getResidues()[0].name() == "LIG"
    assert squashed_system[-1].getResidues()[0].name() == "LIG"
    # First we must have the unperturbed molecules, and then the perturbed ones.
    assert n_atoms == [12, 21, 24, 15, 18, 27, 30]
    python_mapping = {k.value(): v.value() for k, v in mapping.items()}
    assert python_mapping == {0: 0, 2: 1, 3: 2}


def test_squash_multires(perturbed_tripeptide):
    squashed_system, mapping = BSS.Align._squash._squash(perturbed_tripeptide)
    assert len(squashed_system) == 1
    assert len(squashed_system[0].getResidues()) == 4


@pytest.mark.parametrize("is_lambda1", [False, True])
def test_squashed_molecule_mapping(perturbed_system, is_lambda1):
    res = BSS.Align._squash._squashed_molecule_mapping(perturbed_system, is_lambda1=is_lambda1)
    if not is_lambda1:
        expected = {0: 0, 2: 1, 3: 2, 1: 3, 4: 5}
    else:
        expected = {0: 0, 2: 1, 3: 2, 1: 4, 4: 6}
    assert res == expected


@pytest.mark.parametrize("is_lambda1", [False, True])
def test_squashed_atom_mapping(perturbed_tripeptide, is_lambda1):
    res = BSS.Align._squash._squashed_atom_mapping(perturbed_tripeptide, is_lambda1=is_lambda1)
    if not is_lambda1:
        merged_indices = list(range(16)) + list(range(16, 30)) + list(range(43, 51))
        squashed_indices = list(range(16)) + list(range(16, 30)) + list(range(30, 38))
    else:
        merged_indices = list(range(16)) + list(range(16, 21)) + list(range(23, 26)) + \
                         list(range(30, 43)) + list(range(43, 51))
        squashed_indices = list(range(16)) + list(range(38, 59)) + list(range(30, 38))
    expected = dict(zip(merged_indices, squashed_indices))
    assert res == expected


def test_unsquash(perturbed_system):
    squashed_system, mapping = BSS.Align._squash._squash(perturbed_system)
    new_perturbed_system = BSS.Align._squash._unsquash(perturbed_system, squashed_system, mapping)
    assert [mol0.nAtoms() == mol1.nAtoms() for mol0, mol1 in zip(perturbed_system, new_perturbed_system)]
    assert [mol0.isPerturbable() == mol1.isPerturbable() for mol0, mol1 in zip(perturbed_system, new_perturbed_system)]


def test_unsquash_multires(perturbed_tripeptide):
    squashed_system, mapping = BSS.Align._squash._squash(perturbed_tripeptide)
    new_perturbed_system = BSS.Align._squash._unsquash(perturbed_tripeptide, squashed_system, mapping)
    assert [mol0.nAtoms() == mol1.nAtoms() for mol0, mol1 in zip(perturbed_tripeptide, new_perturbed_system)]
    assert [mol0.isPerturbable() == mol1.isPerturbable() for mol0, mol1 in zip(perturbed_tripeptide, new_perturbed_system)]
