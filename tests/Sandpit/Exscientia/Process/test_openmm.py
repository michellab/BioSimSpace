import math
import re

import BioSimSpace.Sandpit.Exscientia as BSS
import pytest

# Store the tutorial URL.
url = BSS.tutorialUrl()

@pytest.fixture(scope="session")
def system():
    """Re-use the same molecuar system for each test."""
    return BSS.IO.readMolecules(["tests/input/ala.top", "tests/input/ala.crd"])


@pytest.mark.parametrize(
    "restraint,dry_run",
    [("backbone", False), ("heavy", True), ("all", True), ("none", True)]
)
def test_minimise(system, restraint, dry_run):
    """Test a minimisation protocol."""

    # Create a short minimisation protocol.
    protocol = BSS.Protocol.Minimisation(
        steps=100, restraint=restraint, force_constant=1000
        )

    # Run the process, check that it finished without error, and returns a system.
    new_system, process = run_process(
        system, protocol, dry_run=dry_run
        )

    if restraint == "none":
        return

    # Check if restrained atoms stayed in place
    restraint_atoms = system.getRestraintAtoms(restraint)

    expected_lines = [
        f"restrained_atoms = {restraint_atoms}",
        "restraint.addBond(i, j, 0*nanometers, 418400.0 * kilojoules_per_mole / nanometer**2)"
    ]

    check_lines_in_file(expected_lines, f"{process._work_dir}/{process._name}_script.py")

    if new_system is None:
        return

    did_not_move = [
        get_dist_atoms(system.getAtom(x), new_system.getAtom(x)) < 0.05
        for x in restraint_atoms
    ]

    assert all(did_not_move)


@pytest.mark.parametrize(
    "restraint,dry_run",
    [("backbone", False), ("heavy", True), ("all", True), ("none", True)]
)
def test_equilibrate(system, restraint, dry_run):
    """Test an equilibration protocol."""

    # Create a short equilibration protocol.
    protocol = BSS.Protocol.Equilibration(
        runtime=BSS.Types.Time(0.001, "nanoseconds"),
        restraint=restraint,
        force_constant=1000
    )

    # Run the process, check that it finished without error, and returns a system.
    new_system, process = run_process(
        system, protocol, dry_run=dry_run
        )

    if restraint == "none":
        return

    # Check if restrained atoms stayed in place
    restraint_atoms = system.getRestraintAtoms(restraint)

    expected_lines = [
        f"restrained_atoms = {restraint_atoms}",
        "restraint.addBond(i, j, 0*nanometers, 418400.0 * kilojoules_per_mole / nanometer**2)"
    ]

    check_lines_in_file(expected_lines, f"{process._work_dir}/{process._name}_script.py")

    if new_system is None:
        return

    did_not_move = [
        get_dist_atoms(system.getAtom(x), new_system.getAtom(x)) < 0.5
        for x in restraint_atoms
    ]

    assert all(did_not_move)


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


@pytest.mark.parametrize(
    "restraint,dry_run",
    [("backbone", False), ("heavy", True), ("all", True), ("none", True)]
)
def test_production(system, restraint, dry_run):
    """Test a production protocol."""

    # Create a short production protocol.
    protocol = BSS.Protocol.Production(
        runtime=BSS.Types.Time(0.001, "nanoseconds"),
        restraint=restraint,
        force_constant=1000
    )

    # Run the process, check that it finished without error, and returns a system.
    new_system, process = run_process(
        system, protocol, dry_run=dry_run
        )

    if restraint == "none":
        return

    # Check if restrained atoms stayed in place
    restraint_atoms = system.getRestraintAtoms(restraint)

    expected_lines = [
        f"restrained_atoms = {restraint_atoms}",
        "restraint.addBond(i, j, 0*nanometers, 418400.0 * kilojoules_per_mole / nanometer**2)"
    ]

    check_lines_in_file(expected_lines, f"{process._work_dir}/{process._name}_script.py")

    if new_system is None:
        return

    did_not_move = [
        get_dist_atoms(system.getAtom(x), new_system.getAtom(x)) < 1
        for x in restraint_atoms
    ]

    assert all(did_not_move)


def run_process(system, protocol, dry_run=False):
    """Helper function to run various simulation protocols.

    Args:
        system: A BSS test system instance
        protocol: A BSS protocol instance
        file_regression: pytest file-regression fixture to check generated
            OpenMM script
        request: Pytest request fixture to determine name of generated
            OpenMM script
    """

    # Need to adjust protocol defaults to get trajectory frames
    # in short test runs.
    if not isinstance(protocol, BSS.Protocol.Minimisation):
        protocol.setReportInterval(100)
        protocol.setRestartInterval(100)

    # Initialise the OpenMM  process.
    process = BSS.Process.OpenMM(system, protocol, name="test")

    if dry_run:
        return None, process

    # Start the OpenMM simulation.
    process.start()

    # Wait for the process to end.
    process.wait()

    # Make sure the process didn't error.
    assert not process.isError()

    # Make sure that we get a molecular system back.
    new_system = process.getSystem()
    assert new_system is not None

    return new_system, process


def get_dist_atoms(a, b):
    """Return the distance between two `BioSimSpace.Atom` instances"""
    v = a.coordinates() - b.coordinates()
    return math.sqrt(
        sum([
            v.x().angstroms().value()**2,
            v.y().angstroms().value()**2,
            v.z().angstroms().value()**2
            ])
        )


def check_lines_in_file(expected_lines, filename):
    for line in expected_lines:
        with open(filename) as fp:
            for found_line in fp:
                found_line = found_line.strip()
                if found_line == line:
                    break
            else:
                raise AssertionError(
                    f"Expected line not found in script: {line}")
