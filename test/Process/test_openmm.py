import BioSimSpace as BSS

import pytest

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

def test_heat(system):
    """Test a heating protocol."""

    # Create a short heating protocol.
    protocol = BSS.Protocol.Equilibration(runtime=BSS.Types.Time(0.001, "nanoseconds"),
                                          temperature_start=BSS.Types.Temperature(0, "kelvin"),
                                          temperature_end=BSS.Types.Temperature(300, "kelvin"))

    # Run the process and check that it finishes without error.
    assert run_process(system, protocol)

def test_cool(system):
    """Test a cooling protocol."""

    # Create a short heating protocol.
    protocol = BSS.Protocol.Equilibration(runtime=BSS.Types.Time(0.001, "nanoseconds"),
                                          temperature_start=BSS.Types.Temperature(300, "kelvin"),
                                          temperature_end=BSS.Types.Temperature(0, "kelvin"))

    # Run the process and check that it finishes without error.
    assert run_process(system, protocol)

def test_production(system):
    """Test a production protocol."""

    # Create a short production protocol.
    protocol = BSS.Protocol.Production(runtime=BSS.Types.Time(0.001, "nanoseconds"))

    # Run the process and check that it finishes without error.
    assert run_process(system, protocol)

def run_process(system, protocol):
    """Helper function to run various simulation protocols."""

    # Initialise the OpenMM  process.
    process = BSS.Process.OpenMM(system, protocol, name="test")

    # Start the OpenMM simulation.
    process.start()

    # Wait for the process to end.
    process.wait()

    # Return the process exit code.
    return not process.isError()
