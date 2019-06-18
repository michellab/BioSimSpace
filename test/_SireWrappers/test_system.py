import BioSimSpace as BSS

# Glob the input files.
files = BSS.IO.glob("test/io/amber/ala/*")

# Load the molecular system.
# This contains an alanin-dipeptide and 630 water molecules.
system = BSS.IO.readMolecules(files)

def test_atom_reindexing():
    # Search for all oxygen atoms in water molecules water molecules within
    # the system.
    result = system.search("waters and element oxygen")

    # There are 22 atoms in the alanine-dipeptide, then three in each water
    # molecule. The oxygen atoms come first in each water molecule, so
    # absolute indices should start at 22 and increment by 3.

    # Starting index.
    index = 22

    for atom in result.getResults():
        # Ensure the absolute index matches.
        assert system.getIndex(atom) == index

        # Increment the index. (3-point water.)
        index += 3

def test_residue_reindexing():
    # Search for all waters by residue name.
    result = system.search("resname WAT")

    # There are 3 residues in the alanine-dipeptide, then one in each water
    # molecule. This means that residue indexing should start at 3 and
    # increment by 1.

    # Starting index.
    index = 3

    for residue in result.getResults():
        # Ensure the absolute index matches.
        assert system.getIndex(residue) == index

        index += 1

def test_molecule_reindexing():
    # Search for all waters by residue name.
    result = system.search("resname WAT")

    # There are 631 molecules in the system: an alanine-dipeptide, followed by
    # 630 water molecules. This means that molecule indexing should start at 1
    # and increment by 1.

    # Starting index.
    index = 1

    # Note that the waters will be returned as residue objects, since this
    # is the minimal representation, i.e. a water contains a single residue.
    # As such, we convert each result to a molecule.
    for residue in result.getResults():
        # Ensure the absolute index matches.
        assert system.getIndex(residue.toMolecule()) == index

        index += 1
