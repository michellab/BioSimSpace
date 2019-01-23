import BioSimSpace as BSS

import pytest

# Make sure GROMSCS is installed.
has_gromacs = BSS._gmx_exe is not None

@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_minimise():
    """Test a minimisation protocol."""

    # Create a short minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=100)

    # Run the process and check that it finishes without error.
    assert run_process(protocol)

@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_equilibrate():
    """Test an equilibration protocol."""

    # Create a short equilibration protocol.
    protocol = BSS.Protocol.Equilibration(runtime=BSS.Types.Time(0.001, "nanoseconds"))

    # Run the process and check that it finishes without error.
    assert run_process(protocol)

@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_heat():
    """Test a heating protocol."""

    # Create a short heating protocol.
    protocol = BSS.Protocol.Equilibration(runtime=BSS.Types.Time(0.001, "nanoseconds"),
                                          temperature_start=BSS.Types.Temperature(0, "kelvin"),
                                          temperature_end=BSS.Types.Temperature(300, "kelvin"))

    # Run the process and check that it finishes without error.
    assert run_process(protocol)

@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_cool():
    """Test a cooling protocol."""

    # Create a short heating protocol.
    protocol = BSS.Protocol.Equilibration(runtime=BSS.Types.Time(0.001, "nanoseconds"),
                                          temperature_start=BSS.Types.Temperature(300, "kelvin"),
                                          temperature_end=BSS.Types.Temperature(0, "kelvin"))

    # Run the process and check that it finishes without error.
    assert run_process(protocol)

@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_production():
    """Test a production protocol."""

    # Create a short production protocol.
    protocol = BSS.Protocol.Production(runtime=BSS.Types.Time(0.001, "nanoseconds"))

    # Run the process and check that it finishes without error.
    assert run_process(protocol)

def create_process(protocol):
    """Create an Amber process for a given prototol."""

    # Glob the input files.
    files = BSS.IO.glob("test/io/amber/ala/*")

    # Load the molecular system.
    system = BSS.IO.readMolecules(files)

    # Initialise the GROMACS process.
    return BSS.Process.Gromacs(system, protocol, name="test")

def run_process(protocol):
    """Helper function to run various simulation protocols."""

    # Initialise the GROMACS process.
    process = create_process(protocol)

    # Start the GROMACS simulation.
    process.start()

    # Wait for the process to end.
    process.wait()

    # Return the process exit code.
    return not process.isError()
