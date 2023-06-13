import pytest

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Align import make_ml

# Store the tutorial URL.
url = BSS.tutorialUrl()


@pytest.fixture
def mol():
    return BSS.IO.readMolecules(["tests/input/ala.top", "tests/input/ala.crd"])[0]


@pytest.fixture
def mol_ml(mol):
    return make_ml(mol)


@pytest.mark.parametrize("molecule,isML", [("mol", False), ("mol_ml", True)])
def test_isML(molecule, isML, request):
    assert request.getfixturevalue(molecule).isML() is isML


@pytest.mark.parametrize("molecule,n_ML", [("mol", 0), ("mol_ml", 1)])
def test_getMLMolecules(molecule, n_ML, request):
    assert len(request.getfixturevalue(molecule).toSystem().getMLMolecules()) == n_ML


@pytest.mark.parametrize("molecule,n_ML", [("mol", 0), ("mol_ml", 1)])
def test_nMLMolecules(molecule, n_ML, request):
    assert request.getfixturevalue(molecule).toSystem().nMLMolecules() == n_ML
