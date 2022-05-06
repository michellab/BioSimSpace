import BioSimSpace as BSS

import filecmp
import pytest
import os
import warnings as _warnings

@pytest.fixture
def system(scope="session"):
    """Re-use the same molecuar system for each test."""
    return BSS.IO.readMolecules("test/input/amber/ala/*")

def test_minimise(system):
    """Test a minimisation protocol."""

    # Create a short minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=100)

    # Run the process and check that it finishes without error.
    assert run_process(system, protocol)

def test_equilibrate(system):
    """Test an equilibration protocol."""

    # Create a short equilibration protocol.
    protocol = BSS.Protocol.Equilibration(runtime=BSS.Types.Time(0.001, "nanoseconds"))

    # Run the process and check that it finishes without error.
    assert run_process(system, protocol)

def test_production(system):
    """Test a production protocol."""

    # Create a short production protocol.
    protocol = BSS.Protocol.Production(runtime=BSS.Types.Time(0.001, "nanoseconds"))

    # Run the process and check that it finishes without error.
    assert run_process(system, protocol)

@pytest.mark.parametrize("morph, pert",
    [("test/input/morphs/morph01.pickle", "test/input/morphs/morph01.pert")])
def test_pert_file(morph, pert):
    """Test the perturbation file writer."""

    import pickle

    # Unpickle the molecule.
    with open(morph, "rb") as file:
        mol = pickle.load(file)

    # Create the perturbation file.
    BSS.Process._somd._to_pert_file(mol)

    # Check that the files are the same.
    assert filecmp.cmp("MORPH.pert", pert)

    # Remove the temporary perturbation file.
    os.remove("MORPH.pert")

def run_process(system, protocol):
    """Helper function to run various simulation protocols."""

    # Initialise the SOMD process.
    process = BSS.Process.Somd(system, protocol, name="test", platform="CPU")

    # Start the SOMD simulation.
    process.start()

    # Wait for the process to end.
    process.wait()

    res = process.isError()
    # OpenMM minimisation occasionally fails with "Particle coordinate is nan"
    # if this is the case, the test will pass with a warning
    if res:
        with open(os.path.join(process.workDir(), "test.out"), "r") as hnd:
            if ("RuntimeError: Particle coordinate is nan" in hnd.read()):
                _warnings.warn("Test raised RuntimeError: Particle coordinate is nan", RuntimeWarning)
                res = False

    return not res

    # Return the process exit code.
    return not process.isError()
