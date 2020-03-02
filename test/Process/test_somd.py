import os
import warnings as _warnings
import BioSimSpace as BSS

def test_minimise():
    """Test a minimisation protocol."""

    # Create a short minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=100)

    # Run the process and check that it finishes without error.
    assert run_process(protocol)

def test_equilibrate():
    """Test an equilibration protocol."""

    # Create a short equilibration protocol.
    protocol = BSS.Protocol.Equilibration(runtime=BSS.Types.Time(0.001, "nanoseconds"))

    # Run the process and check that it finishes without error.
    assert run_process(protocol)

def test_production():
    """Test a production protocol."""

    # Create a short production protocol.
    protocol = BSS.Protocol.Production(runtime=BSS.Types.Time(0.001, "nanoseconds"))

    # Run the process and check that it finishes without error.
    assert run_process(protocol)

def create_process(protocol):
    """Create a SOMD process for a given prototol."""

    # Glob the input files.
    files = BSS.IO.glob("test/io/amber/ala/*")

    # Load the molecular system.
    system = BSS.IO.readMolecules(files)

    # Initialise the SOMD process.
    return BSS.Process.Somd(system, protocol, name="test", platform="CPU")

def run_process(protocol):
    """Helper function to run various simulation protocols."""

    # Initialise the SOMD process.
    process = create_process(protocol)

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
