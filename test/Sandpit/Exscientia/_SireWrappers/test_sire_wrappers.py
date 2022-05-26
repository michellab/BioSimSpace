import BioSimSpace.Sandpit.Exscientia as BSS

import pytest

@pytest.fixture
def system(scope="session"):
    """Re-use the same molecuar system for each test."""
    return BSS.IO.readMolecules("test/input/amber/ala/*")

# Parameterise the function with a set of molecule indices.
@pytest.mark.parametrize("index", [0, -1])
def test_molecule_equivalence(system, index):
    # Make sure that we get the same molecule, however we extract it from
    # the system.
    assert system[index] == system.getMolecules()[index] == system.getMolecule(index)

def test_iterators(system):
    # Iterate over all molecules in the system, either directly,
    # or after calling getMolecules() and make sure they are equivalent.

    # Get the molecules from the system.
    molecules = system.getMolecules()

    # First iterate over the system object directly.
    for idx, mol in enumerate(system):
        assert mol == molecules[idx] == system.getMolecule(idx)

    # Now iterate over the molecules instead.
    for idx, mol in enumerate(molecules):
        assert mol == system[idx] == system.getMolecule(idx)

    # Now search for all molecules in the system.
    search_result = system.search("all")

    # Iterate over the search result and make sure we match all molecules.
    for idx, mol in enumerate(search_result):
        # Single residue molecules will have been converted to a Residue object.
        if type(mol) is BSS._SireWrappers.Residue:
            assert mol.toMolecule() == system[idx] == molecules[idx] == system.getMolecule(idx)
        else:
            assert mol == system[idx] == molecules[idx] == system.getMolecule(idx)

def test_atom_reindexing(system):
    # Search for all oxygen atoms in water molecules water molecules within
    # the system.
    results = system.search("waters and element oxygen")

    # There are 22 atoms in the alanine-dipeptide, then three in each water
    # molecule. The oxygen atoms come first in each water molecule, so
    # absolute indices should start at 22 and increment by 3.

    # Starting index.
    index = 22

    for atom in results:
        # Ensure the absolute index matches.
        assert system.getIndex(atom) == index

        # Increment the index. (3-point water.)
        index += 3

def test_residue_reindexing(system):
    # Search for all waters by residue name.
    results = system.search("resname WAT")

    # There are 3 residues in the alanine-dipeptide, then one in each water
    # molecule. This means that residue indexing should start at 3 and
    # increment by 1.

    # Starting index.
    index = 3

    for residue in results:
        # Ensure the absolute index matches.
        assert system.getIndex(residue) == index

        index += 1

def test_molecule_reindexing(system):
    # Search for all waters by residue name.
    results = system.search("resname WAT")

    # There are 631 molecules in the system: an alanine-dipeptide, followed by
    # 630 water molecules. This means that molecule indexing should start at 1
    # and increment by 1.

    # Starting index.
    index = 1

    # Note that the waters will be returned as residue objects, since this
    # is the minimal representation, i.e. a water contains a single residue.
    # As such, we convert each result to a molecule.
    for residue in results:
        # Ensure the absolute index matches.
        assert system.getIndex(residue.toMolecule()) == index

        index += 1

def test_contains(system):
    # Extract the first molecule.
    m = system[0]

    # Make sure the molecule is in the system.
    assert m in system

    # Make sure a copy of the molecule isn't in the system.
    assert m.copy() not in system

    # Extract the first residue of the molecule.
    r = m.getResidues()[0]

    # Extract the first atom of the residue.
    a = r[0]

    # Make sure the residue and atom are in the system.
    assert r in system
    assert a in system

    # Make sure the residue and atom are in the molecule.
    assert r in m
    assert a in m

    # Make sure that the atom is in the residue.
    assert a in r

    # Get an atom from a different residue.
    a = system[1].getResidues()[0][0]

    # Make sure the atom isn't in the residue.
    assert a not in r
