import BioSimSpace as BSS

import pytest

# Make sure GROMSCS is installed.
has_gromacs = BSS._gmx_exe is not None

@pytest.fixture
def system(scope="session"):
    """Re-use the same molecuar system for each test."""
    return BSS.IO.readMolecules("test/input/amber/ala/*")

@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_minimise(system):
    """Test a minimisation protocol."""

    # Create a short minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=100)

    # Run the process and check that it finishes without error.
    assert run_process(system, protocol)

@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_equilibrate(system):
    """Test an equilibration protocol."""

    # Create a short equilibration protocol.
    protocol = BSS.Protocol.Equilibration(runtime=BSS.Types.Time(0.001, "nanoseconds"))

    # Run the process and check that it finishes without error.
    assert run_process(system, protocol)

@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_heat(system):
    """Test a heating protocol."""

    # Create a short heating protocol.
    protocol = BSS.Protocol.Equilibration(runtime=BSS.Types.Time(0.001, "nanoseconds"),
                                          temperature_start=BSS.Types.Temperature(0, "kelvin"),
                                          temperature_end=BSS.Types.Temperature(300, "kelvin"))

    # Run the process and check that it finishes without error.
    assert run_process(system, protocol)

@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_cool(system):
    """Test a cooling protocol."""

    # Create a short heating protocol.
    protocol = BSS.Protocol.Equilibration(runtime=BSS.Types.Time(0.001, "nanoseconds"),
                                          temperature_start=BSS.Types.Temperature(300, "kelvin"),
                                          temperature_end=BSS.Types.Temperature(0, "kelvin"))

    # Run the process and check that it finishes without error.
    assert run_process(system, protocol)

@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_production(system):
    """Test a production protocol."""

    # Create a short production protocol.
    protocol = BSS.Protocol.Production(runtime=BSS.Types.Time(0.001, "nanoseconds"))

    # Run the process and check that it finishes without error.
    assert run_process(system, protocol)

def run_process(system, protocol):
    """Helper function to run various simulation protocols."""

    # Initialise the GROMACS process.
    process = BSS.Process.Gromacs(system, protocol, name="test")

    # Only run on a single MPI rank.
    process.setArg("-ntmpi", 1)

    # Start the GROMACS simulation.
    process.start()

    # Wait for the process to end.
    process.wait()

    # Return the process exit code.
    return not process.isError()
