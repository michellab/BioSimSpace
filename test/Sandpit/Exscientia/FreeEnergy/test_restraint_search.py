import pytest

import numpy as np

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Align import decouple
from BioSimSpace.Sandpit.Exscientia.FreeEnergy import (
    RestraintSearch,
    Restraint,
)
from BioSimSpace.Sandpit.Exscientia.Trajectory import Trajectory
from BioSimSpace.Sandpit.Exscientia.Units.Length import nanometer
from BioSimSpace.Sandpit.Exscientia.Units.Angle import degree

from BioSimSpace.Sandpit.Exscientia._Utils import _try_import, _have_imported

mda = _try_import("MDAnalysis")

has_mdanalysis = _have_imported(mda)


_MDRestraintsGenerator = _try_import(
    "MDRestraintsGenerator",
    install_command="pip install MDRestraintsGenerator",
)

is_MDRestraintsGenerator = _have_imported(_MDRestraintsGenerator)

# Make sure GROMACS is installed.
has_gromacs = BSS._gmx_exe is not None

# Store the tutorial URL.
url = BSS.tutorialUrl()


@pytest.mark.skipif(
    (has_gromacs is False or has_mdanalysis is False),
    reason="Requires GROMACS and MDAnalysis to be installed.",
)
def test_run_Gromacs():
    """Test if the normal run works on Gromacs."""
    ligand = BSS.IO.readMolecules(
        [
            f"{url}/ligand04.prm7.bz2",
            f"{url}/ligand04.rst7.bz2",
        ]
    ).getMolecule(0)
    decouple_ligand = decouple(ligand)
    protocol = BSS.Protocol.Production(runtime=BSS.Types.Time(8, "FEMTOSECOND"))
    restraint_search = RestraintSearch(
        decouple_ligand.toSystem(), protocol=protocol, engine="GROMACS"
    )
    restraint_search.start()
    restraint_search.wait()
    assert not restraint_search._process.isError()


@pytest.mark.skipif(
    (
        is_MDRestraintsGenerator is False
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
