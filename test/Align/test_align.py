from Sire.MM import InternalFF, IntraCLJFF, IntraFF
from Sire.Mol import AtomIdx, PartialMolecule

import BioSimSpace as BSS

import pytest

# Parameterise the function with a set of valid atom prematches.
@pytest.mark.parametrize("prematch", [{AtomIdx(3) : AtomIdx(1)},
                                      {AtomIdx(5) : AtomIdx(9)},
                                      {AtomIdx(4) : AtomIdx(5)},
                                      {AtomIdx(1) : AtomIdx(0)}])
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


# Parameterise the function with a set of invalid atom prematches.
@pytest.mark.parametrize("prematch", [{AtomIdx(-1) : AtomIdx(1)},
                                      {AtomIdx(50) : AtomIdx(9)},
                                      {AtomIdx(4) : AtomIdx(48)},
                                      {AtomIdx(1) : AtomIdx(-1)}])
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
    s0 = BSS.IO.readMolecules(BSS.IO.glob("test/io/ligands/ligand01*"))
    s1 = BSS.IO.readMolecules(BSS.IO.glob("test/io/ligands/ligand02*"))

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Get the best mapping between the molecules.
    mapping = BSS.Align.matchAtoms(m0, m1, timeout=BSS.Units.Time.second)

    # Align m0 to m1 based on the mapping.
    m0 = BSS.Align.rmsdAlign(m0, m1, mapping)

    # Create the merged molecule.
    m2 = BSS.Align.merge(m0, m1, mapping)

    # Store the number of atoms in m0.
    n0 = m0._sire_molecule.nAtoms()

    # Test that the intramolecular energies area the same.

    # IntraCLJFF:
    #  Old interface. Uses the "intrascale" matrix. Validate that this
    #  is consistent.
    # IntraFF:
    #  New interface. Uses atom "connectivity". Validate that the bonding
    #  is consistent.

    intraclj0 = IntraCLJFF("intraclj")
    intraclj0.add(m0._sire_molecule)

    intraff0 = IntraFF("intraclj")
    intraff0.add(m0._sire_molecule)

    intraclj1 = IntraCLJFF("intraclj")
    intraclj1.add(m1._sire_molecule)

    intraff1 = IntraFF("intraclj")
    intraff1.add(m1._sire_molecule)

    intraclj2 = IntraCLJFF("intraclj")
    intraff2 = IntraFF("intraclj")

    # Create maps between property names: { "prop" : "prop0" }, { "prop" : "prop1" }
    pmap0 = {}
    pmap1 = {}
    for prop in m2._sire_molecule.propertyKeys():
        if prop[-1] == "0":
            pmap0[prop[:-1]] = prop
        elif prop[-1] == "1":
            pmap1[prop[:-1]] = prop

    intraclj2.add(m2._sire_molecule, pmap0)
    intraff2.add(m2._sire_molecule, pmap0)

    assert intraclj0.energy().value() == pytest.approx(intraclj2.energy().value())
    assert intraff0.energy().value() == pytest.approx(intraff2.energy().value())

    intraclj2 = IntraCLJFF("intraclj")
    intraff2 = IntraFF("intraclj")

    intraclj2.add(m2._sire_molecule, pmap1)
    intraff2.add(m2._sire_molecule, pmap1)

    assert intraclj1.energy().value() == pytest.approx(intraclj2.energy().value())
    assert intraff1.energy().value() == pytest.approx(intraff2.energy().value())

    # Test that the internal energies are consistent. This will validate that
    # bond, angle, dihedral, and improper energies are correct.

    internalff0 = InternalFF("internal")
    internalff0.setStrict(True)
    internalff0.add(m0._sire_molecule)

    internalff1 = InternalFF("internal")
    internalff1.setStrict(True)
    internalff1.add(m1._sire_molecule)

    # First extract a partial molecule using the atoms from molecule0 in
    # the merged molecule.
    selection = m2._sire_molecule.selection()
    selection.deselectAll()
    for atom in m0._sire_molecule.atoms():
        selection.select(atom.index())
    partial_mol = PartialMolecule(m2._sire_molecule, selection)

    internalff2 = InternalFF("internal")
    internalff2.setStrict(True)
    internalff2.add(partial_mol, pmap0)

    assert internalff0.energy().value() == pytest.approx(internalff2.energy().value())

    # Now extract a partial molecule using the atoms from molecule1 in
    # the merged molecule.
    selection = m2._sire_molecule.selection()
    selection.deselectAll()
    for idx in mapping.keys():
        selection.select(idx)
    for idx in range(n0, m2._sire_molecule.nAtoms()):
        selection.select(AtomIdx(idx))
    partial_mol = PartialMolecule(m2._sire_molecule, selection)

    internalff2 = InternalFF("internal")
    internalff2.setStrict(True)
    internalff2.add(partial_mol, pmap1)

    assert internalff1.energy().value() == pytest.approx(internalff2.energy().value())
