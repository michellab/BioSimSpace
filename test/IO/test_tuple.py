import BioSimSpace as BSS

def test_tuple():
    """Check that we can read from a tuple of files."""
    BSS.IO.readMolecules(("test/input/ala.crd", "test/input/ala.top"))
