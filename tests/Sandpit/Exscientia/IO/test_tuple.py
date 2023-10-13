import BioSimSpace.Sandpit.Exscientia as BSS
from tests.conftest import root_fp


def test_tuple():
    """Check that we can read from a tuple of files."""
    BSS.IO.readMolecules((f"{root_fp}/input/ala.crd", f"{root_fp}/input/ala.top"))
