import pytest
import itertools

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

def test_set_property_map(mol):
    '''Test if the charge and LJ are set properly'''
    for values in itertools.product(*[(True, False) for _ in range(4)]):
        property_map0 = {"charge": values[0], "LJ": values[1]}
        property_map1 = {"charge": values[2], "LJ": values[3]}
        new = decouple(mol, property_map0=property_map0,
                        property_map1=property_map1)
        assert new._sire_object.property('charge0').value() == values[0]
        assert new._sire_object.property('charge1').value() == values[2]
        assert new._sire_object.property('LJ0').value() == values[1]
        assert new._sire_object.property('LJ1').value() == values[3]
        assert new.isDecoupled() == True
        assert new._sire_object.property('is_decoupled').value() == True
        assert new._sire_object.property('intramol').value() == True

def test_ommit_property_map(mol):
    '''Test if the charge and LJ are set to True by default'''
    new = decouple(mol)
    assert new._sire_object.property('charge0').value() == True
    assert new._sire_object.property('charge1').value() == True
    assert new._sire_object.property('LJ0').value() == True
    assert new._sire_object.property('LJ1').value() == True
    assert new.isDecoupled() == True
    assert new._sire_object.property('is_decoupled').value() == True
    assert new._sire_object.property('intramol').value() == True

def test_intramol(mol):
    '''Test the case of intramol=False'''
    new = decouple(mol, intramol=False, property_map0={}, property_map1={})
    assert new.isDecoupled() == True
    assert new._sire_object.property('is_decoupled').value() == True
    assert new._sire_object.property('intramol').value() == False

def test_getDecoupledMolecules(mol):
    '''Test the case of intramol=False'''
    new = decouple(mol)
    decoupled_mol_list = new.toSystem().getDecoupledMolecules()
    assert isinstance(decoupled_mol_list[0], BSS._SireWrappers.Molecule)
    assert new.toSystem().nDecoupledMolecules() == 1

