import BioSimSpace as BSS

from Sire.Base import findExe

import os
import pytest

# Make sure AMBER is installed.
if BSS._amber_home is not None:
    exe = "%s/bin/sander" % BSS._amber_home
    if os.path.isfile(exe):
        has_amber = True
    else:
        has_amber = False
else:
    has_amber = False

# Make sure GROMACS is installed.
if BSS._gmx_exe is None:
    has_gromacs = False
else:
    has_gromacs = True

# Make sure NAMD is installed.
try:
    findExe("namd2")
    has_namd = True
except:
    has_namd = False

@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
def test_amber():
    """Test a short AMBER minimisation protocol with the MD driver."""

    # Create a short minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=100)

    # Glob the input files.
    files = BSS.IO.glob("test/io/amber/ala/*")

    # Load the molecular system.
    system = BSS.IO.readMolecules(files)

    # Initialise the AMBER process.
    process = BSS.MD.run(system, protocol, name="test")

    # Wait for the process to end.
    process.wait()

    # Check that the process finishes without error.
    assert not process.isError()

@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_gromacs():
    """Test a short GROMACS minimisation protocol with the MD driver."""

    # Create a short minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=100)

    # Glob the input files.
    files = BSS.IO.glob("test/io/gromacs/kigaki/*")

    # Load the molecular system.
    system = BSS.IO.readMolecules(files)

    # Initialise the GROMACS process.
    process = BSS.MD.run(system, protocol, name="test")

    # Wait for the process to end.
    process.wait()

    # Check that the process finishes without error.
    assert not process.isError()

@pytest.mark.skipif(has_namd is False, reason="Requires NAMD to be installed.")
def test_namd():
    """Test a short NAMD minimisation protocol with the MD driver."""

    # Create a short minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=100)

    # Glob the input files.
    files = BSS.IO.glob("test/io/namd/alanin/*")

    # Load the molecular system.
    system = BSS.IO.readMolecules(files)

    # Initialise the NAMD process.
    process = BSS.MD.run(system, protocol, name="test")

    # Wait for the process to end.
    process.wait()

    # Check that the process finishes without error.
    assert not process.isError()
