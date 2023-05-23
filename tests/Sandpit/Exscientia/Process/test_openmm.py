import BioSimSpace.Sandpit.Exscientia as BSS

from BioSimSpace.Sandpit.Exscientia._Utils import _try_import, _have_imported

import pytest

# Make sure GROMACS is installed.
has_gromacs = BSS._gmx_exe is not None

# Make sure openff is installed.
_openff = _try_import("openff")
has_openff = _have_imported(_openff)

# Store the tutorial URL.
url = BSS.tutorialUrl()


@pytest.fixture(scope="session")
def system():
    """Re-use the same molecuar system for each test."""
    return BSS.IO.readMolecules(["tests/input/ala.top", "tests/input/ala.crd"])


def test_minimise(system):
    """Test a minimisation protocol."""

    # Create a short minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=100)

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


@pytest.mark.parametrize("restraint", ["backbone", "heavy", "all", "none"])
def test_equilibrate(system, restraint):
    """Test an equilibration protocol."""

    # Create a short equilibration protocol.
    protocol = BSS.Protocol.Equilibration(
        runtime=BSS.Types.Time(0.001, "nanoseconds"), restraint=restraint
    )

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


def test_heat(system):
    """Test a heating protocol."""

    # Create a short heating protocol.
    protocol = BSS.Protocol.Equilibration(
        runtime=BSS.Types.Time(0.001, "nanoseconds"),
        temperature_start=BSS.Types.Temperature(0, "kelvin"),
        temperature_end=BSS.Types.Temperature(300, "kelvin"),
    )

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


def test_cool(system):
    """Test a cooling protocol."""

    # Create a short heating protocol.
    protocol = BSS.Protocol.Equilibration(
        runtime=BSS.Types.Time(0.001, "nanoseconds"),
        temperature_start=BSS.Types.Temperature(300, "kelvin"),
        temperature_end=BSS.Types.Temperature(0, "kelvin"),
    )

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


def test_production(system):
    """Test a production protocol."""

    # Create a short production protocol.
    protocol = BSS.Protocol.Production(runtime=BSS.Types.Time(0.001, "nanoseconds"))

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


@pytest.mark.skipif(
    has_gromacs is False or has_openff is False,
    reason="Requires GROMACS and OpenFF to be installed",
)
def test_rhombic_dodecahedron():
    """Test that OpenMM can load and run rhombic dodecahedral triclinic spaces."""

    # Create a methane molecule.
    mol = BSS.Parameters.openff_unconstrained_2_0_0("C").getMolecule()

    # Generate box dimensions and angles for a hexagonal rhombic dodecahedron.
    box, angles = BSS.Box.rhombicDodecahedronHexagon(5 * BSS.Units.Length.nanometer)

    # Create a solvated system.
    solvated = BSS.Solvent.tip3p(mol, box=box, angles=angles)

    # Create a short minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=100)

    # Run the process, check that it finished without error, and returns a system.
    run_process(solvated, protocol)


def run_process(system, protocol):
    """Helper function to run various simulation protocols."""

    # Need to adjust protocol defaults to get trajectory frames
    # in short test runs.
    if not isinstance(protocol, BSS.Protocol.Minimisation):
        protocol.setReportInterval(100)
        protocol.setRestartInterval(100)

    # Initialise the OpenMM  process.
    process = BSS.Process.OpenMM(system, protocol, name="test")

    # Start the OpenMM simulation.
    process.start()

    # Wait for the process to end.
    process.wait()

    # Make sure the process didn't error.
    assert not process.isError()

    # Make sure that we get a molecular system back.
    assert process.getSystem() is not None
