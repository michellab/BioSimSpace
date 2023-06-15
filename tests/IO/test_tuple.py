import BioSimSpace as BSS


def test_tuple():
    """Check that we can read from a tuple of files."""
    BSS.IO.readMolecules(("tests/input/ala.crd", "tests/input/ala.top"))
