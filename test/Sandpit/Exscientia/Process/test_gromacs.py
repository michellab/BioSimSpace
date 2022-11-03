import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Align import decouple
from BioSimSpace.Sandpit.Exscientia.Units.Length import angstrom
from BioSimSpace.Sandpit.Exscientia.Units.Angle import radian
from BioSimSpace.Sandpit.Exscientia.Units.Energy import kcal_per_mol
from BioSimSpace.Sandpit.Exscientia.FreeEnergy import Restraint
from BioSimSpace.Sandpit.Exscientia.Units.Temperature import kelvin

import pytest

# Make sure GROMSCS is installed.
has_gromacs = BSS._gmx_exe is not None

@pytest.fixture
def system(scope="session"):
    """Re-use the same molecuar system for each test."""
    return BSS.IO.readMolecules("test/Sandpit/Exscientia/input/amber/ala/*")

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
    assert run_process(system, protocol, extra_options={'verlet-buffer-tolerance': '2e-07'})

@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_production(system):
    """Test a production protocol."""

    # Create a short production protocol.
    protocol = BSS.Protocol.Production(runtime=BSS.Types.Time(0.001, "nanoseconds"))

    # Run the process and check that it finishes without error.
    assert run_process(system, protocol)

@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_vacuum_water(system):
    """Regression test for ensuring the water topology is swapped for vacuum simulations."""

    # Create a short production protocol.
    protocol = BSS.Protocol.Minimisation(steps=100)

    # Create a new system using the first two molecules of 'system'.
    # This will be an alanine-dipeptide and water in vacuum.
    new_system = (system[0] + system[1]).toSystem()

    # Run the process and check that it finishes without error.
    assert run_process(new_system, protocol)

@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
@pytest.mark.parametrize("restraint", ["backbone", "heavy"])
def test_restraints(restraint):
    """Regression test for correct injection of restraint file into GROMACS topology."""

    # Load the perturbable system.
    system = BSS.IO.readPerturbableSystem(
        "test/Sandpit/Exscientia/input/morphs/complex_vac0.prm7",
        "test/Sandpit/Exscientia/input/morphs/complex_vac0.rst7",
        "test/Sandpit/Exscientia/input/morphs/complex_vac1.prm7",
        "test/Sandpit/Exscientia/input/morphs/complex_vac1.rst7"
    )

    # Create an equilibration protocol with backbone restraints.
    protocol = BSS.Protocol.Equilibration(restraint=restraint)

    # Create the simulation process.
    process = BSS.Process.Gromacs(system, protocol)


def test_write_restraint(system, tmp_path):
    """Test if the restraint has been written in a way that could be processed
    correctly."""
    ligand = ligand = BSS.IO.readMolecules(BSS.IO.glob("test/input/ligands/ligand01*")).getMolecule(0)
    decoupled_ligand = decouple(ligand)
    l1 = decoupled_ligand.getAtoms()[0]
    l2 = decoupled_ligand.getAtoms()[1]
    l3 = decoupled_ligand.getAtoms()[2]
    ligand_2 = BSS.IO.readMolecules(BSS.IO.glob("test/input/ligands/ligand04*")).getMolecule(0)
    r1 = ligand_2.getAtoms()[0]
    r2 = ligand_2.getAtoms()[1]
    r3 = ligand_2.getAtoms()[2]
    system = (decoupled_ligand+ligand_2).toSystem()

    restraint_dict = {
        "anchor_points":{"r1":r1, "r2":r2, "r3":r3,
                         "l1":l1, "l2":l2, "l3":l3},
        "equilibrium_values":{"r0":7.84 * angstrom,
                              "thetaA0":0.81 * radian,
                              "thetaB0":1.74 * radian,
                              "phiA0":2.59 * radian,
                              "phiB0":-1.20 * radian,
                              "phiC0":2.63 * radian},
        "force_constants":{"kr":10 * kcal_per_mol / angstrom ** 2,
                           "kthetaA":10 * kcal_per_mol / (radian * radian),
                           "kthetaB":10 * kcal_per_mol / (radian * radian),
                           "kphiA":10 * kcal_per_mol / (radian * radian),
                           "kphiB":10 * kcal_per_mol / (radian * radian),
                           "kphiC":10 * kcal_per_mol / (radian * radian)}}
    restraint = Restraint(system, restraint_dict, 300 * kelvin, rest_type='Boresch')

    # Create a short production protocol.
    protocol = BSS.Protocol.Production(runtime=BSS.Types.Time(0.0001, "nanoseconds"))

    # Run the process and check that it finishes without error.
    assert run_process(system, protocol, restraint=restraint,
                       work_dir=str(tmp_path))
    with open(tmp_path / 'test.top', 'r') as f:
        assert 'intermolecular_interactions' in f.read()


def run_process(system, protocol, **kwargs):
    """Helper function to run various simulation protocols."""

    # Initialise the GROMACS process.
    process = BSS.Process.Gromacs(system, protocol, name="test", **kwargs)

    # Only run on a single MPI rank.
    process.setArg("-ntmpi", 1)

    # Start the GROMACS simulation.
    process.start()

    # Wait for the process to end.
    system = process.getSystem(block=True)
    assert system is not None

    # Return the process exit code.
    return not process.isError()
