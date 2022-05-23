import pytest

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Align._decouple import decouple


@pytest.fixture()
def mol():
    # Benzene.
    mol = BSS.Parameters.openff_unconstrained_2_0_0(
        "c1ccccc1").getMolecule()
    return mol

def test_sanity(mol):
    assert mol.isDecoupled() == False

@pytest.mark.parametrize("charge0", [True, False])
@pytest.mark.parametrize("charge1", [True, False])
@pytest.mark.parametrize("LJ0", [True, False])
@pytest.mark.parametrize("LJ1", [True, False])
def test_set_property_map(mol, charge0, charge1, LJ0, LJ1):
    '''Test if the charge and LJ are set properly'''
    property_map0 = {"charge": charge0, "LJ": LJ0}
    property_map1 = {"charge": charge1, "LJ": LJ1}
    new = decouple(mol, property_map0=property_map0,
                    property_map1=property_map1)
    assert new._sire_object.property('charge0').value() == charge0
    assert new._sire_object.property('charge1').value() == charge1
    assert new._sire_object.property('LJ0').value() == LJ0
    assert new._sire_object.property('LJ1').value() == LJ1
    assert new.isDecoupled() == True
    assert new._sire_object.property('annihilated').value() == False

def test_ommit_property_map(mol):
    '''Test if the charge and LJ are set to True by default'''
    new = decouple(mol)
    assert new._sire_object.property('charge0').value() == True
    assert new._sire_object.property('charge1').value() == False
    assert new._sire_object.property('LJ0').value() == True
    assert new._sire_object.property('LJ1').value() == False
    assert new.isDecoupled() == True
    assert new._sire_object.property('annihilated').value() == False

def test_intramol(mol):
    '''Test the case of intramol=False'''
    new = decouple(mol, intramol=False, property_map0={}, property_map1={})
    assert new.isDecoupled() == True
    assert new._sire_object.property('annihilated').value() == True

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
