try:
    from alchemlyb.parsing.gmx import extract_u_nk

    is_alchemlyb = True
except ModuleNotFoundError:
    is_alchemlyb = False

from matplotlib import cm
import numpy as np
import os
import pandas as pd
import pathlib
import pytest

import time

from BioSimSpace.Sandpit.Exscientia.Align._decouple import decouple
from BioSimSpace.Sandpit.Exscientia.FreeEnergy import Restraint
from BioSimSpace.Sandpit.Exscientia.Protocol import FreeEnergyEquilibration
from BioSimSpace.Sandpit.Exscientia.Units.Angle import radian, degree
from BioSimSpace.Sandpit.Exscientia.Units.Energy import kcal_per_mol
from BioSimSpace.Sandpit.Exscientia.Units.Length import angstrom
from BioSimSpace.Sandpit.Exscientia.Units.Temperature import kelvin
from BioSimSpace.Sandpit.Exscientia import Types as _Types
import BioSimSpace.Sandpit.Exscientia as BSS


# Make sure GROMSCS is installed.
has_gromacs = BSS._gmx_exe is not None

# Get the input files
url = BSS.tutorialUrl()


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
@pytest.mark.skipif(is_alchemlyb is False, reason="Requires alchemlyb to be installed.")
class Test_gmx_ABFE:
    @staticmethod
    @pytest.fixture(scope="class")
    def freenrg():
        m = BSS.IO.readMolecules(
            [
                f"{url}/crd.gro.bz2",
                f"{url}/complex.top.bz2",
            ]
        ).getMolecule(
            1
        )  # Molecule 1 is the ligand

        # Assign atoms for restraint
        atom_1 = m.getAtoms()[0]
        atom_2 = m.getAtoms()[1]
        atom_3 = m.getAtoms()[2]
        atom_4 = m.getAtoms()[3]
        atom_5 = m.getAtoms()[4]
        atom_6 = m.getAtoms()[5]

        decoupled_m = decouple(m)
        solvated = BSS.Solvent.tip3p(
            molecule=decoupled_m, box=3 * [3 * BSS.Units.Length.nanometer]
        )

        # Create random restraint dictionary
        restraint_dict = {
            "anchor_points": {
                "r1": atom_1,
                "r2": atom_2,
                "r3": atom_3,
                "l1": atom_4,
                "l2": atom_5,
                "l3": atom_6,
            },
            "equilibrium_values": {
                "r0": 5.08 * angstrom,
                "thetaA0": 64.051 * degree,
                "thetaB0": 39.618 * degree,
                "phiA0": 2.59 * radian,
                "phiB0": -1.20 * radian,
                "phiC0": 2.63 * radian,
            },
            "force_constants": {
                "kr": 10 * kcal_per_mol / angstrom**2,
                "kthetaA": 10 * kcal_per_mol / (radian * radian),
                "kthetaB": 10 * kcal_per_mol / (radian * radian),
                "kphiA": 10 * kcal_per_mol / (radian * radian),
                "kphiB": 10 * kcal_per_mol / (radian * radian),
                "kphiC": 10 * kcal_per_mol / (radian * radian),
            },
        }

        restraint = Restraint(
            solvated, restraint_dict, 298 * kelvin, restraint_type="Boresch"
        )
        protocol = FreeEnergyEquilibration(
            lam=pd.Series(data={"coul": 0.5, "vdw": 0.0}),
            lam_vals=pd.DataFrame(
                data={
                    "coul": [0.0, 0.5, 1.0, 1.0, 1.0],
                    "vdw": [0.0, 0.0, 0.0, 0.5, 1.0],
                }
            ),
            runtime=_Types.Time(0, "nanoseconds"),
        )
        freenrg = BSS.FreeEnergy.Absolute(
            solvated, protocol, engine="GROMACS", restraint=restraint
        )
        freenrg.run()
        freenrg.wait()
        return freenrg

    def test_file_exist(self, freenrg):
        """Test if all the files are there."""
        path = pathlib.Path(freenrg.workDir())
        for i in range(5):
            assert (path / f"lambda_{i}" / "gromacs.xvg").is_file()

    @pytest.mark.skipif(is_alchemlyb is False, reason="Need alchemlyb.")
    def test_lambda(self, freenrg):
        """Test if the xvg files contain the correct lambda."""
        path = pathlib.Path(freenrg.workDir())
        for i, (coul, vdw) in enumerate(
            zip([0.0, 0.5, 1.0, 1.0, 1.0], [0.0, 0.0, 0.0, 0.5, 1.0])
        ):
            u_nk = extract_u_nk(path / f"lambda_{i}" / "gromacs.xvg", 300)
            assert u_nk.index.names == ["time", "coul-lambda", "vdw-lambda"]
            assert np.isclose(u_nk.index.values[0][1], coul)
            assert np.isclose(u_nk.index.values[0][2], vdw)


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
class Test_Somd_ABFE:
    @staticmethod
    @pytest.fixture(scope="class")
    def freenrg():
        # Just use a single ligand with anchors based all on atoms in the ligand, for simplicity
        m = BSS.IO.readMolecules(
            [f"{url}/ligand01.prm7.bz2", f"{url}/ligand01.rst7.bz2"]
        ).getMolecule(0)

        # Assign atoms for restraint
        atom_1 = m.getAtoms()[0]
        atom_2 = m.getAtoms()[1]
        atom_3 = m.getAtoms()[2]
        atom_4 = m.getAtoms()[3]
        atom_5 = m.getAtoms()[4]
        atom_6 = m.getAtoms()[5]

        decoupled_m = decouple(m)
        solvated = BSS.Solvent.tip3p(
            molecule=decoupled_m, box=3 * [3 * BSS.Units.Length.nanometer]
        )

        # Create random restraint dictionary
        restraint_dict = {
            "anchor_points": {
                "r1": atom_1,
                "r2": atom_2,
                "r3": atom_3,
                "l1": atom_4,
                "l2": atom_5,
                "l3": atom_6,
            },
            "equilibrium_values": {
                "r0": 5.08 * angstrom,
                "thetaA0": 64.051 * degree,
                "thetaB0": 39.618 * degree,
                "phiA0": 2.59 * radian,
                "phiB0": -1.20 * radian,
                "phiC0": 2.63 * radian,
            },
            "force_constants": {
                "kr": 10 * kcal_per_mol / angstrom**2,
                "kthetaA": 10 * kcal_per_mol / (radian * radian),
                "kthetaB": 10 * kcal_per_mol / (radian * radian),
                "kphiA": 10 * kcal_per_mol / (radian * radian),
                "kphiB": 10 * kcal_per_mol / (radian * radian),
                "kphiC": 10 * kcal_per_mol / (radian * radian),
            },
        }

        restraint = Restraint(
            solvated, restraint_dict, 298 * kelvin, restraint_type="Boresch"
        )
        protocol = BSS.Protocol.FreeEnergy(
            lam_vals=[0.0, 0.5, 1.0], runtime=0.0001 * BSS.Units.Time.nanosecond
        )
        freenrg = BSS.FreeEnergy.Absolute(
            solvated, protocol, engine="SOMD", restraint=restraint
        )

        freenrg.run()
        freenrg.wait()
        # Sleep to allow files to be written
        time.sleep(25)
        return freenrg

    def test_files_exist(self, freenrg):
        """Test if the files have been created. Note that e.g. gradients.dat
        are not created until later in the simulation, so their presence is
        not tested for."""
        path = pathlib.Path(freenrg.workDir())
        for lam in ["0.0000", "0.5000", "1.0000"]:
            assert (path / f"lambda_{lam}" / "simfile.dat").is_file()
            assert (path / f"lambda_{lam}" / "SYSTEM.s3").is_file()
            assert (path / f"lambda_{lam}" / "somd.cfg").is_file()
            assert (path / f"lambda_{lam}" / "somd.rst7").is_file()
            assert (path / f"lambda_{lam}" / "somd.prm7").is_file()
            assert (path / f"lambda_{lam}" / "somd.err").is_file()
            assert (path / f"lambda_{lam}" / "somd.out").is_file()

    # @pytest.mark.xfail(reason="freenrg.wait() doesn't work for SOMD, so files aren't created in time.")
    def test_correct_conf_file(self, freenrg):
        """Check that lambda data is correct in somd.cfg"""
        path = pathlib.Path(freenrg.workDir())
        for lam in ["0.0000", "0.5000", "1.0000"]:
            with open(os.path.join(path, f"lambda_{lam}", "somd.cfg"), "rt") as f:
                lines = f.readlines()
                assert "lambda array = 0.0, 0.5, 1.0\n" in lines
                assert f"lambda_val = {float(lam):.1f}\n" in lines
