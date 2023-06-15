import pathlib

import pytest
import pandas as pd
import numpy as np
import bz2

try:
    from alchemlyb.parsing.gmx import extract_u_nk

    is_alchemlyb = True
except ModuleNotFoundError:
    is_alchemlyb = False

try:
    from alchemtest.gmx import load_ABFE
    from alchemtest.amber import load_bace_example

    is_alchemtest = True
except ModuleNotFoundError:
    is_alchemtest = False

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Protocol import FreeEnergyEquilibration
from BioSimSpace.Sandpit.Exscientia.Align._decouple import decouple
from BioSimSpace.Sandpit.Exscientia import Types as _Types

# Make sure GROMSCS is installed.
has_gromacs = BSS._gmx_exe is not None

# Store the tutorial URL.
url = BSS.tutorialUrl()


@pytest.mark.skipif(
    is_alchemtest is False, reason="Requires alchemtest and alchemlyb to be installed."
)
class TestRelativeAnalysis:
    @staticmethod
    @pytest.fixture(scope="class")
    def gmx_data(tmp_path_factory):
        outdir = tmp_path_factory.mktemp("gromacs")
        for leg in ["complex", "ligand"]:
            for index, filename in enumerate(load_ABFE().data[leg]):
                dir = outdir / leg / f"lambda_{index}"
                dir.mkdir(parents=True)
                eng = dir / "gromacs.xvg"
                eng.symlink_to(filename)
        return outdir

    @staticmethod
    @pytest.fixture(scope="class")
    def amber_data(tmp_path_factory):
        outdir = tmp_path_factory.mktemp("amber")
        for leg, stage in [("complex", "decharge"), ("solvated", "vdw")]:
            for index, filename in enumerate(load_bace_example().data[leg][stage]):
                dir = outdir / f"{leg}_{stage}" / f"lambda_{index}"
                dir.mkdir(parents=True)
                with open(dir / "amber.out", "w") as f:
                    with bz2.open(filename, "rt") as bz_file:
                        f.write(bz_file.read())
        return outdir

    @staticmethod
    @pytest.fixture(scope="class")
    def gmx_complex(gmx_data):
        outdir = gmx_data
        complex, _ = BSS.FreeEnergy.Relative.analyse(
            work_dir=str(outdir / "complex"),
            temperature=310 * BSS.Units.Temperature.kelvin,
        )
        return complex

    @staticmethod
    @pytest.fixture(scope="class")
    def gmx_ligand(gmx_data):
        outdir = gmx_data
        ligand, _ = BSS.FreeEnergy.Relative.analyse(
            work_dir=str(outdir / "ligand"),
            temperature=310 * BSS.Units.Temperature.kelvin,
        )
        return ligand

    @staticmethod
    @pytest.fixture(scope="class")
    def amber_complex_decharge(amber_data):
        outdir = amber_data
        complex, _ = BSS.FreeEnergy.Relative.analyse(
            work_dir=str(outdir / "complex_decharge"),
            temperature=298 * BSS.Units.Temperature.kelvin,
        )
        return complex

    @staticmethod
    @pytest.fixture(scope="class")
    def amber_solvated_vdw(amber_data):
        outdir = amber_data
        complex, _ = BSS.FreeEnergy.Relative.analyse(
            work_dir=str(outdir / "solvated_vdw"),
            temperature=298 * BSS.Units.Temperature.kelvin,
        )
        return complex

    @pytest.mark.parametrize(
        "fixture,length,energy",
        [
            ("gmx_ligand", 20, 7.654472744451637),
            ("gmx_complex", 30, 21.819752),
            ("amber_complex_decharge", 5, -5.25352),
            ("amber_solvated_vdw", 12, 2.261816),
        ],
    )
    def test_pmf(self, fixture, length, energy, request):
        pmf = request.getfixturevalue(fixture)
        assert len(pmf) == length
        assert len(pmf[0]) == 3
        np.testing.assert_allclose(
            pmf[-1][1] / BSS.Units.Energy.kcal_per_mol, energy, atol=0.1
        )

    def test_difference(self, gmx_complex, gmx_ligand):
        dG, error = BSS.FreeEnergy.Relative.difference(gmx_complex, gmx_ligand)
        np.testing.assert_allclose(
            dG / BSS.Units.Energy.kcal_per_mol, 14.216101, atol=0.1
        )
