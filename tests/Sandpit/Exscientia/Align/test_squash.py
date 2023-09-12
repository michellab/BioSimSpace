import pickle
import pytest

import BioSimSpace.Sandpit.Exscientia as BSS

from tests.Sandpit.Exscientia.conftest import has_amber


@pytest.fixture(scope="session")
def perturbed_system():
    # N_atoms are: 12, 15, 18, 21, 24, 27 and 30.
    mol_smiles = [
        "c1ccccc1",
        "c1ccccc1C",
        "c1ccccc1CC",
        "c1ccccc1CCC",
        "c1ccccc1CCCC",
        "c1ccccc1CCCCC",
        "c1ccccc1CCCCCC",
    ]
    mols = [BSS.Parameters.gaff(smi).getMolecule() for smi in mol_smiles]
    pert_mols = [
        mols[0],
        BSS.Align.merge(mols[1], mols[2]),
        mols[3],
        mols[4],
        BSS.Align.merge(mols[5], mols[6]),
    ]
    system = BSS._SireWrappers.System(pert_mols)
    return system


@pytest.fixture
def perturbed_tripeptide():
    return pickle.load(
        open("tests/Sandpit/Exscientia/input/merged_tripeptide.pickle", "rb")
    )


@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
@pytest.mark.parametrize(
    "explicit,expected_n_atoms",
    [
        (False, [12, 21, 24, 15, 18, 27, 30]),
        (True, [12, 21, 24, 18, 18, 30, 30]),
    ],
)
def test_squash(perturbed_system, explicit, expected_n_atoms):
    squashed_system, mapping = BSS.Align._squash._squash(
        perturbed_system, explicit_dummies=explicit
    )
    assert len(squashed_system) == 7
    n_atoms = [mol.nAtoms() for mol in squashed_system]
    assert squashed_system[-2].getResidues()[0].name() == "LIG"
    assert squashed_system[-1].getResidues()[0].name() == "LIG"
    # First we must have the unperturbed molecules, and then the perturbed ones.
    assert n_atoms == expected_n_atoms
    python_mapping = {k.value(): v.value() for k, v in mapping.items()}
    assert python_mapping == {0: 0, 2: 1, 3: 2}


@pytest.mark.parametrize("explicit", [False, True])
def test_squash_multires(perturbed_tripeptide, explicit):
    squashed_system, mapping = BSS.Align._squash._squash(
        perturbed_tripeptide, explicit_dummies=explicit
    )
    assert len(squashed_system) == 1
    assert len(squashed_system[0].getResidues()) == 4


@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
@pytest.mark.parametrize("is_lambda1", [False, True])
def test_squashed_molecule_mapping(perturbed_system, is_lambda1):
    res = BSS.Align._squash._squashed_molecule_mapping(
        perturbed_system, is_lambda1=is_lambda1
    )
    if not is_lambda1:
        expected = {0: 0, 2: 1, 3: 2, 1: 3, 4: 5}
    else:
        expected = {0: 0, 2: 1, 3: 2, 1: 4, 4: 6}
    assert res == expected


@pytest.mark.parametrize("is_lambda1", [False, True])
def test_squashed_atom_mapping_implicit(perturbed_tripeptide, is_lambda1):
    res = BSS.Align._squash._squashed_atom_mapping(
        perturbed_tripeptide, is_lambda1=is_lambda1, explicit_dummies=False
    )
    if not is_lambda1:
        merged_indices = list(range(16)) + list(range(16, 30)) + list(range(43, 51))
        squashed_indices = list(range(16)) + list(range(16, 30)) + list(range(30, 38))
    else:
        merged_indices = (
            list(range(16))
            + list(range(16, 21))
            + list(range(23, 26))
            + list(range(30, 43))
            + list(range(43, 51))
        )
        squashed_indices = list(range(16)) + list(range(38, 59)) + list(range(30, 38))
    expected = dict(zip(merged_indices, squashed_indices))
    assert res == expected


@pytest.mark.parametrize("is_lambda1", [False, True])
def test_squashed_atom_mapping_explicit(perturbed_tripeptide, is_lambda1):
    res = BSS.Align._squash._squashed_atom_mapping(
        perturbed_tripeptide, is_lambda1=is_lambda1, explicit_dummies=True
    )
    merged_indices = list(range(51))
    if not is_lambda1:
        squashed_indices = list(range(16)) + list(range(16, 43)) + list(range(43, 51))
    else:
        squashed_indices = list(range(16)) + list(range(51, 78)) + list(range(43, 51))
    expected = dict(zip(merged_indices, squashed_indices))
    assert res == expected


@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
@pytest.mark.parametrize("explicit", [False, True])
def test_unsquash(perturbed_system, explicit):
    squashed_system, mapping = BSS.Align._squash._squash(
        perturbed_system, explicit_dummies=explicit
    )
    new_perturbed_system = BSS.Align._squash._unsquash(
        perturbed_system, squashed_system, mapping, explicit_dummies=explicit
    )
    assert [
        mol0.nAtoms() == mol1.nAtoms()
        for mol0, mol1 in zip(perturbed_system, new_perturbed_system)
    ]
    assert [
        mol0.isPerturbable() == mol1.isPerturbable()
        for mol0, mol1 in zip(perturbed_system, new_perturbed_system)
    ]


@pytest.mark.parametrize("explicit", [False, True])
def test_unsquash_multires(perturbed_tripeptide, explicit):
    squashed_system, mapping = BSS.Align._squash._squash(
        perturbed_tripeptide, explicit_dummies=explicit
    )
    new_perturbed_system = BSS.Align._squash._unsquash(
        perturbed_tripeptide, squashed_system, mapping, explicit_dummies=explicit
    )
    assert [
        mol0.nAtoms() == mol1.nAtoms()
        for mol0, mol1 in zip(perturbed_tripeptide, new_perturbed_system)
    ]
    assert [
        mol0.isPerturbable() == mol1.isPerturbable()
        for mol0, mol1 in zip(perturbed_tripeptide, new_perturbed_system)
    ]


@pytest.fixture(scope="module")
def benzene():
    return BSS.Parameters.gaff("c1ccccc1").getMolecule()


@pytest.fixture(scope="module")
def decoupled_system(benzene):
    mol = BSS.Align.decouple(benzene, charge=(True, False), LJ=(True, False))
    return mol.copy()


@pytest.fixture(scope="module")
def redecoupled_system(benzene):
    return BSS.Align.decouple(benzene, charge=(False, True), LJ=(False, True))


@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
@pytest.mark.parametrize(
    "mol,is_lambda1,dummies,expect",
    [
        ("decoupled_system", False, True, True),
        ("decoupled_system", True, True, False),
        ("redecoupled_system", False, True, False),
        ("redecoupled_system", True, True, True),
        ("decoupled_system", False, False, False),
        ("decoupled_system", True, False, False),
        ("redecoupled_system", False, False, False),
        ("redecoupled_system", True, False, False),
    ],
)
def test_squashed_molecule_mapping_decouple_lambda(
    mol, is_lambda1, expect, dummies, request
):
    """Test if the result is generated for the correct side of the lambda.
    If mol is being decoupled, the mapping should be in lambda_0 but not in lambda_1.
    If mol is being recoupled, the mapping should be in lambda_1 but not in lambda_0.
    If dummies is False, don't return any mapping."""
    mol = request.getfixturevalue(mol)
    mapping = BSS.Align._squash._squashed_atom_mapping(
        mol, dummies=dummies, is_lambda1=is_lambda1
    )
    if expect:
        assert mapping[0] == 0
        assert mapping[11] == 11
    else:
        assert len(mapping) == 0


@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
@pytest.mark.parametrize(
    "first,second,third,order",
    [
        ("decouple", "merge", "normal", (0, 2, 1)),
        ("decouple", "normal", "merge", (0, 1, 2)),
        ("merge", "decouple", "normal", (2, 0, 1)),
        ("merge", "normal", "decouple", (2, 0, 1)),
        ("normal", "merge", "decouple", (0, 2, 1)),
        ("normal", "decouple", "merge", (0, 1, 2)),
    ],
)
def test_squashed_molecule_mapping_order(
    first, second, third, order, benzene, decoupled_system
):
    """Test if the molecules are in the correct order. For example, if I have a system
    that is in the order of perturbed molecule (0), decoupled molecule (1) and normal molecule (2),
    the order of the molecules in the unsquashed system would be 0, 1, 2.
    After squashing, the order would be 1, 2, 0. This test checks if the correct order has been generated.
    """
    mol_list = []
    for mol in [first, second, third]:
        if mol == "normal":
            mol_list.append(benzene)
        elif mol == "merge":
            mol_list.append(BSS.Align.merge(benzene, benzene))
        else:
            mol_list.append(decoupled_system)
    system = BSS._SireWrappers.System(mol_list)
    mapping = BSS.Align._squash._squashed_atom_mapping(system)
    for atom_idx, squashed_idx in enumerate(order):
        # The index is 0-based so need the offset
        assert mapping[(atom_idx + 1) * 12 - 1] == (squashed_idx + 1) * 12 - 1
