import pytest

import MDAnalysis as mda
import numpy as np

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Align import decouple
from BioSimSpace.Sandpit.Exscientia.FreeEnergy import RestraintSearch, Restraint
from BioSimSpace.Sandpit.Exscientia.Trajectory import Trajectory
from BioSimSpace.Sandpit.Exscientia.Units.Length import nanometer
from BioSimSpace.Sandpit.Exscientia.Units.Angle import degree

from BioSimSpace.Sandpit.Exscientia._Utils import _try_import, _have_imported

_MDRestraintsGenerator = _try_import(
    "MDRestraintsGenerator", install_command="pip install MDRestraintsGenerator"
)

is_MDRestraintsGenerator = _have_imported(_MDRestraintsGenerator)

# Make sure GROMSCS is installed.
has_gromacs = BSS._gmx_exe is not None


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_run_Gromacs():
    """Test if the normal run works on Gromacs."""
    ligand = BSS.IO.readMolecules(
        BSS.IO.glob("test/input/ligands/ligand04*")
    ).getMolecule(0)
    decouple_ligand = decouple(ligand)
    protocol = BSS.Protocol.Production(runtime=BSS.Types.Time(8, "FEMTOSECOND"))
    restraint_search = RestraintSearch(
        decouple_ligand.toSystem(), protocol=protocol, engine="GROMACS"
    )
    restraint_search.start()
    restraint_search.wait()
    assert not restraint_search._process.isError()


class Trajectory(Trajectory):
    def __init__(self):
        pass

    def getTrajectory(self, format="mdanalysis"):
        return mda.Universe(
            "test/Sandpit/Exscientia/input/protein_ligand/complex.tpr",
            "test/Sandpit/Exscientia/input/protein_ligand/traj.xtc",
        )


@pytest.mark.skipif(
    (is_MDRestraintsGenerator is False or has_gromacs is False),
    reason="Requires MDRestraintsGenerator and Gromacs to be installed.",
)
class TestMDRestraintsGenerator_analysis:
    @staticmethod
    @pytest.fixture(scope="class")
    def restraint_search(tmp_path_factory):
        outdir = tmp_path_factory.mktemp("out")
        system = BSS.IO.readMolecules(
            [
                "test/Sandpit/Exscientia/input/protein_ligand/crd.gro",
                "test/Sandpit/Exscientia/input/protein_ligand/complex.top",
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
        restraint_search._process.getTrajectory = lambda: Trajectory()
        restraint = restraint_search.analyse(
            method="MDRestraintsGenerator", restraint_type="Boresch", block=False
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
