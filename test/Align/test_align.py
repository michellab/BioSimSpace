from Sire.MM import InternalFF, IntraCLJFF, IntraFF
from Sire.Mol import AtomIdx, PartialMolecule

import BioSimSpace as BSS

import pytest

# Parameterise the function with a set of valid atom pre-matches.
@pytest.mark.parametrize("prematch", [{3 : 1},
                                      {5 : 9},
                                      {4 : 5},
                                      {1 : 0}])
def test_prematch(prematch):
    # Load the ligands.
    s0 = BSS.IO.readMolecules(BSS.IO.glob("test/io/ligands/ligand01*"))
    s1 = BSS.IO.readMolecules(BSS.IO.glob("test/io/ligands/ligand02*"))

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Get the best mapping between the molecules that contains the prematch.
    mapping = BSS.Align.matchAtoms(m0, m1, timeout=BSS.Units.Time.second,
                                   prematch=prematch)

    # Check that the prematch key:value pair is in the mapping.
    for key, value in prematch.items():
        assert mapping[key] == value

# Parameterise the function with a set of invalid atom pre-matches.
@pytest.mark.parametrize("prematch", [{-1 :  1},
                                      {50 :  9},
                                      { 4 : 48},
                                      { 1 : -1}])
def test_invalid_prematch(prematch):
    # Load the ligands.
    s0 = BSS.IO.readMolecules(BSS.IO.glob("test/io/ligands/ligand01*"))
    s1 = BSS.IO.readMolecules(BSS.IO.glob("test/io/ligands/ligand02*"))

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Assert that the invalid prematch raises a ValueError.
    with pytest.raises(ValueError):
        mapping = BSS.Align.matchAtoms(m0, m1, timeout=BSS.Units.Time.second,
                                       prematch=prematch)

def test_merge():
    # Load the ligands.
    s0 = BSS.IO.readMolecules(BSS.IO.glob("test/io/ligands/ligand31*"))
    s1 = BSS.IO.readMolecules(BSS.IO.glob("test/io/ligands/ligand38*"))

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Get the best mapping between the molecules.
    mapping = BSS.Align.matchAtoms(m0, m1, timeout=BSS.Units.Time.second)

    # Align m0 to m1 based on the mapping.
    m0 = BSS.Align.rmsdAlign(m0, m1, mapping)

    # Create the merged molecule.
    m2 = BSS.Align.merge(m0, m1, mapping, allow_ring_breaking=True)

    # Store the number of atoms in m0.
    n0 = m0._sire_object.nAtoms()

    # Test that the intramolecular energies area the same.

    # IntraCLJFF:
    #  Old interface. Uses the "intrascale" matrix. Validate that this
    #  is consistent.
    # IntraFF:
    #  New interface. Uses atom "connectivity". Validate that the bonding
    #  is consistent.

    intraclj0 = IntraCLJFF("intraclj")
    intraclj0.add(m0._sire_object)

    intraff0 = IntraFF("intraclj")
    intraff0.add(m0._sire_object)

    intraclj1 = IntraCLJFF("intraclj")
    intraclj1.add(m1._sire_object)

    intraff1 = IntraFF("intraclj")
    intraff1.add(m1._sire_object)

    intraclj2 = IntraCLJFF("intraclj")
    intraff2 = IntraFF("intraclj")

    # Create maps between property names: { "prop" : "prop0" }, { "prop" : "prop1" }
    pmap0 = {}
    pmap1 = {}
    for prop in m2._sire_object.propertyKeys():
        if prop[-1] == "0":
            pmap0[prop[:-1]] = prop
        elif prop[-1] == "1":
            pmap1[prop[:-1]] = prop

    intraclj2.add(m2._sire_object, pmap0)
    intraff2.add(m2._sire_object, pmap0)

    assert intraclj0.energy().value() == pytest.approx(intraclj2.energy().value())
    assert intraff0.energy().value() == pytest.approx(intraff2.energy().value())

    intraclj2 = IntraCLJFF("intraclj")
    intraff2 = IntraFF("intraclj")

    intraclj2.add(m2._sire_object, pmap1)
    intraff2.add(m2._sire_object, pmap1)

    assert intraclj1.energy().value() == pytest.approx(intraclj2.energy().value())
    assert intraff1.energy().value() == pytest.approx(intraff2.energy().value())

    # Test that the internal energies are consistent. This will validate that
    # bond, angle, dihedral, and improper energies are correct.

    internalff0 = InternalFF("internal")
    internalff0.setStrict(True)
    internalff0.add(m0._sire_object)

    internalff1 = InternalFF("internal")
    internalff1.setStrict(True)
    internalff1.add(m1._sire_object)

    # First extract a partial molecule using the atoms from molecule0 in
    # the merged molecule.
    selection = m2._sire_object.selection()
    selection.deselectAll()
    for atom in m0._sire_object.atoms():
        selection.select(atom.index())
    partial_mol = PartialMolecule(m2._sire_object, selection)

    internalff2 = InternalFF("internal")
    internalff2.setStrict(True)
    internalff2.add(partial_mol, pmap0)

    assert internalff0.energy().value() == pytest.approx(internalff2.energy().value())

    # Now extract a partial molecule using the atoms from molecule1 in
    # the merged molecule.
    selection = m2._sire_object.selection()
    selection.deselectAll()
    for idx in mapping.keys():
        selection.select(AtomIdx(idx))
    for idx in range(n0, m2._sire_object.nAtoms()):
        selection.select(AtomIdx(idx))
    partial_mol = PartialMolecule(m2._sire_object, selection)

    internalff2 = InternalFF("internal")
    internalff2.setStrict(True)
    internalff2.add(partial_mol, pmap1)

    assert internalff1.energy().value() == pytest.approx(internalff2.energy().value())

def test_ring_breaking_three_membered():
    # Load the ligands.
    s0 = BSS.IO.readMolecules(BSS.IO.glob("test/io/ligands/CAT-13a*"))
    s1 = BSS.IO.readMolecules(BSS.IO.glob("test/io/ligands/CAT-17g*"))

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Generate the mapping.
    mapping = BSS.Align.matchAtoms(m0, m1)

    # Align m0 to m1 based on the mapping.
    m0 = BSS.Align.rmsdAlign(m0, m1, mapping)

    # Try to merge the molecule without allowing ring breaking.
    with pytest.raises(BSS._Exceptions.IncompatibleError):
        m2 = BSS.Align.merge(m0, m1, mapping)

    # Now check that we can merge if we allow ring breaking.
    m2 = BSS.Align.merge(m0, m1, mapping, allow_ring_breaking=True)

def test_ring_breaking_five_membered():
    # Load the ligands.
    s0 = BSS.IO.readMolecules(BSS.IO.glob("test/io/ligands/ligand31*"))
    s1 = BSS.IO.readMolecules(BSS.IO.glob("test/io/ligands/ligand04*"))

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Load the pre-defined mapping.
    mapping = BSS.Align.matchAtoms(m0, m1)

    # Align m0 to m1 based on the mapping.
    m0 = BSS.Align.rmsdAlign(m0, m1, mapping)

    # Try to merge the molecule without allowing ring breaking.
    with pytest.raises(BSS._Exceptions.IncompatibleError):
        m2 = BSS.Align.merge(m0, m1, mapping)

    # Now check that we can merge if we allow ring breaking.
    m2 = BSS.Align.merge(m0, m1, mapping, allow_ring_breaking=True)

def test_ring_breaking_six_membered():
    # Load the ligands.
    s0 = BSS.IO.readMolecules(BSS.IO.glob("test/io/ligands/ligand31*"))
    s1 = BSS.IO.readMolecules(BSS.IO.glob("test/io/ligands/ligand38*"))

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Load the pre-defined mapping.
    mapping = BSS.Align.matchAtoms(m0, m1)

    # Align m0 to m1 based on the mapping.
    m0 = BSS.Align.rmsdAlign(m0, m1, mapping)

    # Try to merge the molecule without allowing ring breaking.
    with pytest.raises(BSS._Exceptions.IncompatibleError):
        m2 = BSS.Align.merge(m0, m1, mapping)

    # Now check that we can merge if we allow ring breaking.
    m2 = BSS.Align.merge(m0, m1, mapping, allow_ring_breaking=True)

@pytest.mark.parametrize("ligands", [["CAT-13c", "CAT-17i"],
                                     ["CAT-13e", "CAT-17g"]])
def test_ring_size_change(ligands):
    # Load the ligands.
    s0 = BSS.IO.readMolecules(BSS.IO.glob("test/io/ligands/%s.*" % ligands[0]))
    s1 = BSS.IO.readMolecules(BSS.IO.glob("test/io/ligands/%s.*" % ligands[1]))

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Generate the mapping.
    mapping = BSS.Align.matchAtoms(m0, m1)

    # Align m0 to m1 based on the mapping.
    m0 = BSS.Align.rmsdAlign(m0, m1, mapping)

    # Try to merge the molecule without allowing ring breaking.
    with pytest.raises(BSS._Exceptions.IncompatibleError):
        m2 = BSS.Align.merge(m0, m1, mapping)

    # Now check that we can merge if we allow ring breaking.
    m2 = BSS.Align.merge(m0, m1, mapping, allow_ring_breaking=True, allow_ring_size_change=True)

# Parameterise the function with a valid mapping.
@pytest.mark.parametrize("mapping", [{
    2 : 21,
    4 : 23,
    6 : 25,
    8 : 27,
    10 : 18,
    1 : 19,
    0 : 20,
    11 : 16,
    12 : 17,
    13 : 14,
    15 : 13,
    18 : 11,
    20 : 9,
    22 : 8,
    23 : 5,
    16 : 6,
    24 : 3,
    26 : 1,
    27 : 0,
    9 : 28,
    5 : 24,
    3 : 22,
    7 : 26,
    14 : 15,
    19 : 12,
    21 : 10,
    17 : 7,
    25 : 4}])
def test_grow_whole_ring(mapping):
    # Load the ligands.
    s0 = BSS.IO.readMolecules(BSS.IO.glob("test/io/ligands/grow1*"))
    s1 = BSS.IO.readMolecules(BSS.IO.glob("test/io/ligands/grow2*"))

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Align m0 to m1 based on the mapping.
    m0 = BSS.Align.rmsdAlign(m0, m1, mapping)

    # Check that we can merge without allowing ring breaking.
    m2 = BSS.Align.merge(m0, m1, mapping)
