import BioSimSpace as BSS

from sire.legacy import Mol as SireMol

import pytest


@pytest.fixture(scope="session")
def system():
    return BSS.IO.readMolecules(["test/input/ala.crd", "test/input/ala.top"])


def test_system(system):
    """
    Check that system conversions work as expected.
    """

    # Store the number of molecules, residues, and atoms.
    num_mol = system.nMolecules()
    num_res = system.nResidues()
    num_atm = system.nAtoms()

    # Convert to Sire format.
    sire_system = BSS.Convert.toSire(system)

    assert sire_system.nMolecules() == num_mol
    assert sire_system.nResidues() == num_res
    assert sire_system.nAtoms() == num_atm

    # Convert to RDKit format. The result will be a container of
    # RDKit molecules.
    rdmols = BSS.Convert.toRDKit(sire_system)

    assert len(rdmols) == num_mol

    # Work out atom total. (No concept of a residue.)
    nats = 0
    for rdmol in rdmols:
        nats += rdmol.GetNumAtoms()

    assert nats == num_atm

    # Convert back to BioSimSpace.
    bss_system = BSS.Convert.toBioSimSpace(rdmols).toSystem()

    assert bss_system.nMolecules() == num_mol
    assert bss_system.nAtoms() == num_atm


def test_residue(system):
    """
    Check that we can convert single residues.
    """

    # Convert to Sire format.
    sire_residue = BSS.Convert.toSire(system[0].getResidues()[0])

    assert isinstance(sire_residue, SireMol.Residue)

    # Now convert back.
    bss_residue = BSS.Convert.toBioSimSpace(sire_residue)

    assert isinstance(bss_residue, BSS._SireWrappers.Residue)


def test_atom(system):
    """
    Check that we can convert single atoms.
    """

    # Convert to Sire format.
    sire_atom = BSS.Convert.toSire(system[0].getAtoms()[0])

    assert isinstance(sire_atom, SireMol.Atom)

    # Now convert back.
    bss_atom = BSS.Convert.toBioSimSpace(sire_atom)

    assert isinstance(bss_atom, BSS._SireWrappers.Atom)
