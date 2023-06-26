import os
import pytest

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Process._process import Process

# Make sure AMBER is installed.
if BSS._amber_home is not None:
    exe = "%s/bin/sander" % BSS._amber_home
    if os.path.isfile(exe):
        has_amber = True
    else:
        has_amber = False
else:
    has_amber = False

# Check whether GROMACS is installed.
has_gromacs = BSS._gmx_exe is not None


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
