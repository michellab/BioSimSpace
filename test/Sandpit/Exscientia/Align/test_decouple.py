import pytest

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Align._decouple import decouple

@pytest.fixture()
def mol():
    # Benzene.
    mol = BSS.IO.readMolecules("test/input/amber/ala/*").getMolecule(0)
    return mol

def test_sanity(mol):
    assert mol.isDecoupled() is False

@pytest.mark.parametrize("charge0", [True, False])
@pytest.mark.parametrize("charge1", [True, False])
@pytest.mark.parametrize("LJ0", [True, False])
@pytest.mark.parametrize("LJ1", [True, False])
def test_set_property_map(mol, charge0, charge1, LJ0, LJ1):
    '''Test if the charge and LJ are set properly'''
    new = decouple(mol, charge=(charge0, charge1), LJ=(LJ0, LJ1))
    decouple_dict = new._sire_object.property("decouple")
    assert bool(decouple_dict["charge"][0]) is charge0
    assert bool(decouple_dict["charge"][1]) is charge1
    assert bool(decouple_dict["LJ"][0]) is LJ0
    assert bool(decouple_dict["LJ"][1]) is LJ1
    assert new.isDecoupled() is True
    assert decouple_dict['intramol'].value() is True

def test_ommit_property_map(mol):
    '''Test if the charge and LJ are set to True by default'''
    new = decouple(mol)
    decouple_dict = new._sire_object.property("decouple")
    assert bool(decouple_dict["charge"][0]) is True
    assert bool(decouple_dict["charge"][1]) is False
    assert bool(decouple_dict["LJ"][0]) is True
    assert bool(decouple_dict["LJ"][1]) is False
    assert new.isDecoupled() is True
    assert decouple_dict['intramol'].value() is True

def test_intramol(mol):
    '''Test the case of intramol=False'''
    new = decouple(mol, intramol=False)
    decouple_dict = new._sire_object.property("decouple")
    assert new.isDecoupled() is True
    assert decouple_dict['intramol'].value() is False

def test_getDecoupledMolecules(mol):
    '''Test the method of DecoupledMolecules'''
    new = decouple(mol)
    decoupled_mol_list = new.toSystem().getDecoupledMolecules()
    assert isinstance(decoupled_mol_list[0], BSS._SireWrappers.Molecule)
    assert new.toSystem().nDecoupledMolecules() == 1

def test_no_DecoupledMolecules(mol):
    '''Test the method of DecoupledMolecules when there is no decoupled mols.'''
    assert mol.toSystem().nDecoupledMolecules() == 0

def test_topology(mol, tmp_path):
    '''Test if the decoupled file could be written correctly.'''
    new = decouple(mol)
    BSS.IO.saveMolecules(str(tmp_path / 'topol'), new.toSystem(), 'grotop')
    assert (tmp_path / 'topol.top').is_file()
