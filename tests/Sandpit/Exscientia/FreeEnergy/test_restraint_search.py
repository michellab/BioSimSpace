import pytest

import numpy as np

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Align import decouple
from BioSimSpace.Sandpit.Exscientia._Exceptions import AnalysisError
from BioSimSpace.Sandpit.Exscientia.FreeEnergy import (
    RestraintSearch,
    Restraint,
)
from BioSimSpace.Sandpit.Exscientia.Trajectory import Trajectory
from BioSimSpace.Sandpit.Exscientia.Units.Angle import degree
from BioSimSpace.Sandpit.Exscientia.Units.Energy import kcal_per_mol
from BioSimSpace.Sandpit.Exscientia.Units.Length import angstrom
from BioSimSpace.Sandpit.Exscientia.Units.Length import nanometer

from tests.Sandpit.Exscientia.conftest import (
    url,
    has_gromacs,
    has_mdanalysis,
    has_mdrestraints_generator,
)


@pytest.fixture(scope="module")
def setup_system():
    "Setup the system for the tests"
    ligand = BSS.IO.readMolecules(
        [
            f"{url}/ligand04.prm7.bz2",
            f"{url}/ligand04.rst7.bz2",
        ]
    ).getMolecule(0)
    decouple_system = decouple(ligand).toSystem()
    protocol = BSS.Protocol.Production(runtime=BSS.Types.Time(8, "FEMTOSECOND"))
    return ligand, decouple_system, protocol


# Make sure GROMSCS is installed.
has_gromacs = BSS._gmx_exe is not None


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_run_Gromacs(setup_system):
    "Test if the normal run works on Gromacs"
    ligand, decouple_system, protocol = setup_system
    restraint_search = RestraintSearch(
        decouple_system, protocol=protocol, engine="GROMACS"
    )
    restraint_search.start()
    restraint_search.wait()
    assert not restraint_search._process.isError()


def test_run_Somd(setup_system):
    "Test if the normal run works with SOMD"
    ligand, decouple_system, protocol = setup_system
    restraint_search = RestraintSearch(
        decouple_system, protocol=protocol, engine="SOMD"
    )
    restraint_search.start()
    restraint_search.wait()
    assert not restraint_search._process.isError()


@pytest.mark.skipif(
    (
        has_mdrestraints_generator is False
        or has_gromacs is False
        or has_mdanalysis is False
    ),
    reason="Requires MDRestraintsGenerator, MDAnalysis and Gromacs to be installed.",
)
class TestMDRestraintsGenerator_analysis:
    @staticmethod
    @pytest.fixture(scope="class")
    def restraint_search(tmp_path_factory):
        outdir = tmp_path_factory.mktemp("out")
        system = BSS.IO.readMolecules(
            [
                f"{url}/crd.gro.bz2",
                f"{url}/complex.top.bz2",
            ]
        )
        ligand = system.getMolecule(1)
        decoupled_ligand = decouple(ligand)
        protein = system.getMolecule(0)
        new_system = (protein + decoupled_ligand).toSystem()

        protocol = BSS.Protocol.Production()
        restraint_search = BSS.FreeEnergy.RestraintSearch(
            new_system,
            protocol=protocol,
            engine="GROMACS",
            work_dir=str(outdir),
        )
        traj, top = BSS.IO.expand(url, ["traj.xtc", "complex.tpr"], ".bz2")
        restraint_search._process.getTrajectory = lambda: Trajectory(
            trajectory=traj, topology=top
        )
        restraint = restraint_search.analyse(
            method="MDRestraintsGenerator",
            restraint_type="Boresch",
            block=False,
        )
        return restraint, outdir

    def test_sanity(self, restraint_search):
        restraint, outdir = restraint_search
        assert isinstance(restraint, Restraint)

    def test_plots(self, restraint_search):
        """Test if all the plots has been generated correctly."""
        restraint, outdir = restraint_search
        assert (outdir / "bond_1.png").is_file()
        assert (outdir / "angle_1.png").is_file()
        assert (outdir / "angle_2.png").is_file()
        assert (outdir / "dihedral_1.png").is_file()
        assert (outdir / "dihedral_2.png").is_file()
        assert (outdir / "dihedral_3.png").is_file()

    def test_dG_off(self, restraint_search):
        """Test if the restraint generated is valid."""
        restraint, outdir = restraint_search
        assert (outdir / "dG_off.dat").is_file()
        dG = np.loadtxt(outdir / "dG_off.dat")
        assert isinstance(dG, np.ndarray)

    def test_top(self, restraint_search):
        """Test if the restraint generated has the same energy."""
        restraint, outdir = restraint_search
        assert (outdir / "BoreschRestraint.top").is_file()
        with open(outdir / "BoreschRestraint.top", "r") as f:
            assert "intermolecular_interactions" in f.read()


@pytest.mark.skipif(
    (has_gromacs is False or has_mdanalysis is False),
    reason="Requires MDAnalysis and Gromacs to be installed.",
)
class TestBSS_analysis:
    """Test selection of restraints using the inbuilt BSS method."""

    @staticmethod
    @pytest.fixture(scope="class")
    def _restraint_search(tmp_path_factory):
        outdir = tmp_path_factory.mktemp("out")
        system = BSS.IO.readMolecules(
            [
                f"{url}/crd.gro.bz2",
                f"{url}/complex.top.bz2",
            ]
        )
        ligand = system.getMolecule(1)
        decoupled_ligand = decouple(ligand)
        protein = system.getMolecule(0)
        new_system = (protein + decoupled_ligand).toSystem()

        protocol = BSS.Protocol.Production()
        restraint_search = BSS.FreeEnergy.RestraintSearch(
            new_system, protocol=protocol, engine="GROMACS", work_dir=str(outdir)
        )
        traj, top = BSS.IO.expand(url, ["traj.xtc", "complex.tpr"], ".bz2")
        restraint_search._process.getTrajectory = lambda: Trajectory(
            trajectory=traj, topology=top
        )
        return restraint_search, outdir

    @staticmethod
    @pytest.fixture(scope="class")
    def restraint(_restraint_search):
        restraint_search, outdir = _restraint_search
        restraint = restraint_search.analyse(
            method="BSS", restraint_type="Boresch", block=False
        )
        return restraint, outdir

    @staticmethod
    @pytest.fixture(scope="class")
    def restraint_k20(_restraint_search):
        restraint_search, outdir = _restraint_search
        restraint = restraint_search.analyse(
            method="BSS",
            restraint_type="Boresch",
            block=False,
            force_constant=20 * kcal_per_mol / angstrom**2,
        )
        return restraint, outdir

    def test_sanity(self, restraint):
        restraint, _ = restraint
        assert isinstance(restraint, Restraint)

    def test_plots(self, restraint):
        """Test if all the plots have been generated correctly"""
        restraint, outdir = restraint
        assert (outdir / "restraint_idx0_dof_time.png").is_file()
        assert (outdir / "restraint_idx0_dof_hist.png").is_file()

    def test_dG_off(self, restraint):
        """Test if the restraint generated has the same energy"""
        restraint, _ = restraint
        assert np.isclose(-9.8955, restraint.correction.value(), atol=0.01)

    def test_bond(self, restraint):
        restraint, _ = restraint
        equilibrium_values_r0 = (
            restraint._restraint_dict["equilibrium_values"]["r0"] / nanometer
        )
        assert np.isclose(0.6057, equilibrium_values_r0, atol=0.001)

    def test_angles(self, restraint):
        restraint, _ = restraint
        equilibrium_values_thetaA0 = (
            restraint._restraint_dict["equilibrium_values"]["thetaA0"] / degree
        )
        assert np.isclose(140.6085, equilibrium_values_thetaA0, atol=0.001)
        equilibrium_values_thetaB0 = (
            restraint._restraint_dict["equilibrium_values"]["thetaB0"] / degree
        )
        assert np.isclose(56.4496, equilibrium_values_thetaB0, atol=0.001)

    def test_dihedrals(self, restraint):
        restraint, _ = restraint
        equilibrium_values_phiA0 = (
            restraint._restraint_dict["equilibrium_values"]["phiA0"] / degree
        )
        assert np.isclose(21.6173, equilibrium_values_phiA0, atol=0.001)
        equilibrium_values_phiB0 = (
            restraint._restraint_dict["equilibrium_values"]["phiB0"] / degree
        )
        assert np.isclose(-19.4394, equilibrium_values_phiB0, atol=0.001)
        equilibrium_values_phiC0 = (
            restraint._restraint_dict["equilibrium_values"]["phiC0"] / degree
        )
        assert np.isclose(71.3148, equilibrium_values_phiC0, atol=0.001)

    def test_index(self, restraint):
        restraint, _ = restraint
        idxs = {
            k: restraint._restraint_dict["anchor_points"][k].index()
            for k in restraint._restraint_dict["anchor_points"]
        }
        assert idxs == {"r1": 1560, "r2": 1558, "r3": 1562, "l1": 10, "l2": 9, "l3": 11}

    def test_force_constant(self, restraint_k20):
        restraint, _ = restraint_k20
        for force_constant in restraint._restraint_dict["force_constants"].values():
            assert np.isclose(20, force_constant.value(), atol=0.01)

    def test_analysis_failure(self, _restraint_search):
        restraint_search, _ = _restraint_search
        with pytest.raises(AnalysisError):
            restraint = restraint_search.analyse(
                method="BSS",
                restraint_type="Boresch",
                block=False,
                cutoff=0.1 * angstrom,
            )
