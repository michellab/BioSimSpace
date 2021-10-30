import BioSimSpace as BSS

import os
import pytest

@pytest.fixture
def system():
    return BSS.IO.readMolecules("test/io/amber/ala/*")

def test_system(system):
    """Test streaming of a Sire system."""

    # Stream to file.
    BSS.Stream.save(system, "test")

    # Stream from file.
    s = BSS.Stream.load("test.s3")

    # Check that both systems contain the same number of molecules.
    assert system.nMolecules() == s.nMolecules()

    # Check that the molecules in both systems contain the same number
    # of residues and atoms.
    for m0, m1 in zip(system, s):
        assert m0.nResidues() == m1.nResidues()
        assert m0.nAtoms() == m1.nAtoms()

    # Remove the test file.
    os.remove("test.s3")

def test_molecule(system):
    """Test streaming of a Sire molecule."""

    # Extract the first molecule.
    m0 = system[0]

    # Stream to file.
    BSS.Stream.save(m0, "test")

    # Stream from file.
    m1 = BSS.Stream.load("test.s3")

    # Check that the molecules contain the same number of residues and atoms.
    assert m0.nResidues() == m1.nResidues()
    assert m0.nAtoms() == m1.nAtoms()

    # Remove the test file.
    os.remove("test.s3")

def test_molecules(system):
    """Test streaming of a Sire molecule group."""

    # Stream to file.
    BSS.Stream.save(system.getMolecules(), "test")

    # Stream from file.
    m = BSS.Stream.load("test.s3")

    # Check that the system contain the same number of molecules as the
    # molecule group.
    assert system.nMolecules() == len(m)

    # Check that the molecules in the system and molecule group contain the
    # same number of residues and atoms.
    for m0, m1 in zip(system, m):
        assert m0.nResidues() == m1.nResidues()
        assert m0.nAtoms() == m1.nAtoms()

    # Remove the test file.
    os.remove("test.s3")

def test_residue(system):
    """Test streaming of a Sire residue."""

    # Extract the first residue of the first molecule.
    r0 = system[0].getResidues()[0]

    # Stream to file.
    BSS.Stream.save(r0, "test")

    # Stream from file.
    r1 = BSS.Stream.load("test.s3")

    # Check that the residues contain the same number of atoms.
    assert r0.nAtoms() == r1.nAtoms()

    # Remove the test file.
    os.remove("test.s3")

def test_atom(system):
    """Test streaming of a Sire atom."""

    # Extract the first atom of the first molecule.
    a0 = system[0].getAtoms()[0]

    # Stream to file.
    BSS.Stream.save(a0, "test")

    # Stream from file.
    a1 = BSS.Stream.load("test.s3")

    # Check that the atom elements are the same.
    assert a0.element() == a1.element()

    # Remove the test file.
    os.remove("test.s3")

def test_select_result(system):
    """Test streaming of a Sire select result."""

    # Search for all carbon atoms.
    s0 = system.search("element C")

    # Stream to file.
    BSS.Stream.save(s0, "test")

    # Stream from file.
    s1 = BSS.Stream.load("test.s3")

    # Check that the number of search results is the same.
    assert s0.nResults() == s1.nResults()

    # Check that the results have the same element and index.
    for a0, a1 in zip(s0, s1):
        assert a0.index() == a1.index()
        assert a0.element() == a1.element()

    # Remove the test file.
    os.remove("test.s3")
