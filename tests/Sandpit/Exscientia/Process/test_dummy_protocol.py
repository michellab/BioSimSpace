import pytest

from tests.Sandpit.Exscientia.conftest import has_amber, has_gromacs

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Process._process import Process


@pytest.fixture(scope="session")
def system():
    """Re-use the same molecuar system for each test."""
    return BSS.IO.readMolecules(["tests/input/ala.top", "tests/input/ala.crd"])


@pytest.mark.skipif(
    has_amber is False or has_gromacs is False,
    reason="Requires AMBER and GROMACS to be installed.",
)
@pytest.mark.parametrize("process", ["Gromacs", "Amber"])
def test_process(system, process):
    protocol = BSS.Protocol.Dummy()
    process_cls = getattr(BSS.Process, process)
    process = process_cls(system, protocol)
    assert isinstance(process, Process)
