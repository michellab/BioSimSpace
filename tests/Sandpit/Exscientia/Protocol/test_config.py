import pandas as pd
import pytest

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia._Utils import _have_imported, _try_import
from BioSimSpace.Sandpit.Exscientia.Align._decouple import decouple
from BioSimSpace.Sandpit.Exscientia.Protocol import (
    ConfigFactory,
    Equilibration,
    FreeEnergy,
    FreeEnergyEquilibration,
    FreeEnergyMinimisation,
    Production,
)

# Make sure GROMACS is installed.
has_gromacs = BSS._gmx_exe is not None

# Make sure antechamber is installed.
has_antechamber = BSS.Parameters._Protocol._amber._antechamber_exe is not None

# Make sure openff is installed.
_openff = _try_import("openff")
has_openff = _have_imported(_openff)

# Store the tutorial URL.
url = BSS.tutorialUrl()


class TestAmber:
    @pytest.fixture(scope="class")
    def system(self):
        m0 = BSS.IO.readMolecules(
            [f"{url}/CAT-13c.prm7.bz2", f"{url}/CAT-13c.rst7.bz2"]
        ).getMolecule(0)
        return m0.toSystem()

    def test_NVT(self, system):
        config = ConfigFactory(system, Equilibration(pressure=None))
        res = [x.strip().strip(",") for x in config.generateAmberConfig()]
        assert "ntb=1" in res


class TestAmberRBFE:
    @pytest.fixture(scope="class")
    def system(self):
        m0 = BSS.IO.readMolecules(
            [f"{url}/CAT-13c.prm7.bz2", f"{url}/CAT-13c.rst7.bz2"]
        ).getMolecule(0)
        m1 = BSS.IO.readMolecules(
            [f"{url}/CAT-13a.prm7.bz2", f"{url}/CAT-13a.rst7.bz2"]
        ).getMolecule(0)
        atom_mapping = BSS.Align.matchAtoms(m0, m1)
        m0 = BSS.Align.rmsdAlign(m0, m1, atom_mapping)
        merged = BSS.Align.merge(m0, m1)
        return merged.toSystem()

    def test_generate_fep_masks(self, system):
        config = ConfigFactory(system, FreeEnergyMinimisation())
        res = config._generate_amber_fep_masks(0.004)
        expected_res = {
            "noshakemask": '""',
            "scmask1": '"@25-27,29-31"',
            "scmask2": '""',
            "timask1": '"@1-45"',
            "timask2": '"@46-84"',
        }
        assert res == expected_res

    def test_generate_restraint_masks(self, system):
        protocol = FreeEnergyEquilibration(restraint=[0, 24], force_constant=4.3)
        config = ConfigFactory(system, protocol)
        res = [x.strip().strip(",") for x in config.generateAmberConfig()]
        expected_res = {"ntr=1", 'restraintmask="@1,25,46"', "restraint_wt=4.3"}
        assert expected_res.issubset(res)

    @pytest.mark.parametrize(
        "protocol", [Equilibration, Production, FreeEnergyEquilibration, FreeEnergy]
    )
    def test_tau_t(self, system, protocol):
        config = ConfigFactory(system, protocol(tau_t=BSS.Types.Time(2, "picosecond")))
        res = [x.strip().strip(",") for x in config.generateAmberConfig()]
        expected_res = {"gamma_ln=0.50000"}
        assert expected_res.issubset(res)


class TestGromacsRBFE:
    @staticmethod
    @pytest.fixture(scope="class")
    def system():
        m0 = BSS.IO.readMolecules(
            [f"{url}/CAT-13a.prm7.bz2", f"{url}/CAT-13a.rst7.bz2"]
        ).getMolecule(0)
        m1 = BSS.IO.readMolecules(
            [f"{url}/CAT-13c.prm7.bz2", f"{url}/CAT-13c.rst7.bz2"]
        ).getMolecule(0)
        atom_mapping = BSS.Align.matchAtoms(m0, m1)
        m0 = BSS.Align.rmsdAlign(m0, m1, atom_mapping)
        merged = BSS.Align.merge(m0, m1)
        return merged.toSystem()

    @pytest.mark.parametrize(
        "protocol", [Equilibration, Production, FreeEnergyEquilibration, FreeEnergy]
    )
    def test_tau_t(self, system, protocol):
        config = ConfigFactory(system, protocol(tau_t=BSS.Types.Time(2, "picosecond")))
        res = config.generateGromacsConfig()
        expected_res = {"tau-t = 2.00000"}
        assert expected_res.issubset(res)

    @pytest.mark.skipif(
        has_gromacs is False, reason="Requires GROMACS to be installed."
    )
    def test_fep(self, system):
        """Test if the default config writer will write the expected thing."""
        protocol = FreeEnergyMinimisation(
            lam=0.0,
            lam_vals=None,
            min_lam=0.0,
            max_lam=1.0,
            num_lam=11,
            steps=10000,
            perturbation_type="full",
        )
        freenrg = BSS.FreeEnergy.Relative(
            system,
            protocol,
            engine="GROMACS",
        )
        with open(f"{freenrg._work_dir}/lambda_6/gromacs.mdp", "r") as f:
            mdp_text = f.read()
            assert (
                "fep-lambdas          = 0.00000 0.10000 0.20000 0.30000 0.40000 0.50000 0.60000 0.70000 0.80000 0.90000 1.00000"
                in mdp_text
            )
            assert "init-lambda-state = 6" in mdp_text

    @pytest.mark.skipif(
        has_gromacs is False, reason="Requires GROMACS to be installed."
    )
    def test_fep_df(self, system):
        """Test if the default config writer configured with the pd.DataFrame
        will write the expected thing.
        """
        protocol = FreeEnergyMinimisation(
            lam=pd.Series(data={"fep": 0.0}),
            lam_vals=None,
            min_lam=pd.Series(data={"fep": 0.0}),
            max_lam=pd.Series(data={"fep": 1.0}),
            num_lam=11,
            steps=10000,
            perturbation_type="full",
        )
        freenrg = BSS.FreeEnergy.Relative(
            system,
            protocol,
            engine="GROMACS",
        )
        with open(f"{freenrg._work_dir}/lambda_6/gromacs.mdp", "r") as f:
            mdp_text = f.read()
            assert (
                "fep-lambdas          = 0.00000 0.10000 0.20000 0.30000 0.40000 0.50000 0.60000 0.70000 0.80000 0.90000 1.00000"
                in mdp_text
            )
            assert "init-lambda-state = 6" in mdp_text

    @pytest.mark.skipif(
        has_gromacs is False, reason="Requires GROMACS to be installed."
    )
    def test_staged_fep_df(self, system):
        """Test if the multi-stage lambda will be written correctly."""
        protocol = FreeEnergyMinimisation(
            lam=pd.Series(data={"bonded": 0.0, "coul": 0.0, "vdw": 0.0}),
            lam_vals=pd.DataFrame(
                data={
                    "bonded": [0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0],
                    "coul": [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
                    "vdw": [0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0],
                }
            ),
            steps=10000,
            perturbation_type="full",
        )
        freenrg = BSS.FreeEnergy.Relative(
            system,
            protocol,
            engine="GROMACS",
        )
        with open(f"{freenrg._work_dir}/lambda_6/gromacs.mdp", "r") as f:
            mdp_text = f.read()
            assert (
                "bonded-lambdas       = 0.00000 0.25000 0.50000 0.75000 1.00000 1.00000 1.00000"
                in mdp_text
            )
            assert (
                "coul-lambdas         = 0.00000 0.00000 0.00000 0.50000 1.00000 1.00000 1.00000"
                in mdp_text
            )
            assert (
                "vdw-lambdas          = 0.00000 0.00000 0.00000 0.00000 0.00000 0.50000 1.00000"
                in mdp_text
            )
            assert "init-lambda-state = 6" in mdp_text


@pytest.mark.skipif(
    has_antechamber is False or has_openff is False,
    reason="Requires ambertools/antechamber and openff to be installed",
)
class TestGromacsABFE:
    @staticmethod
    @pytest.fixture(scope="class")
    def system():
        # Benzene.
        m0 = BSS.Parameters.openff_unconstrained_2_0_0("c1ccccc1").getMolecule()
        m0._sire_object = m0._sire_object.edit().rename("LIG").molecule().commit()
        protocol = FreeEnergyMinimisation(
            lam=0.0,
            lam_vals=None,
            min_lam=0.0,
            max_lam=1.0,
            num_lam=11,
            steps=10000,
            perturbation_type="full",
        )
        return m0, protocol

    @pytest.mark.skipif(
        has_gromacs is False, reason="Requires GROMACS to be installed."
    )
    def test_decouple_vdw_q(self, system):
        m, protocol = system
        """Test the decoupling where lambda0 = vdw-q and lambda1=none."""
        mol = decouple(m)
        freenrg = BSS.FreeEnergy.Relative(
            mol.toSystem(),
            protocol,
            engine="GROMACS",
        )
        with open(f"{freenrg._work_dir}/lambda_6/gromacs.mdp", "r") as f:
            mdp_text = f.read()
            assert "couple-moltype = LIG" in mdp_text
            assert "couple-lambda0 = vdw-q" in mdp_text
            assert "couple-lambda1 = none" in mdp_text
            assert "couple-intramol = yes" in mdp_text

    @pytest.mark.skipif(
        has_gromacs is False, reason="Requires GROMACS to be installed."
    )
    def test_annihilate_vdw2q(self, system):
        """Test the annihilation where lambda0 = vdw-q and lambda1=none."""
        m, protocol = system
        mol = decouple(m, charge=(False, True), LJ=(True, False), intramol=False)
        freenrg = BSS.FreeEnergy.Relative(
            mol.toSystem(),
            protocol,
            engine="GROMACS",
        )
        with open(f"{freenrg._work_dir}/lambda_6/gromacs.mdp", "r") as f:
            mdp_text = f.read()
            assert "couple-moltype = LIG" in mdp_text
            assert "couple-lambda0 = vdw" in mdp_text
            assert "couple-lambda1 = q" in mdp_text
            assert "couple-intramol = no" in mdp_text

    @pytest.mark.skipif(
        has_gromacs is False, reason="Requires GROMACS to be installed."
    )
    def test_sc_parameters(self, system):
        """Test if the soft core parameters have been written.
        The default sc-alpha is 0, which means the soft-core of the vdw is not
        turned on by default. This checks if the this value has been changed to
        0.5.
        """
        m, protocol = system
        mol = decouple(m)
        freenrg = BSS.FreeEnergy.Relative(
            mol.toSystem(),
            protocol,
            engine="GROMACS",
        )

        with open(f"{freenrg._work_dir}/lambda_6/gromacs.mdp", "r") as f:
            mdp_text = f.read()
            assert "sc-alpha = 0.5" in mdp_text


class TestAmberASFE:
    @staticmethod
    @pytest.fixture(scope="class")
    def system():
        mol = BSS.IO.readMolecules(
            [f"{url}/CAT-13c.prm7.bz2", f"{url}/CAT-13c.rst7.bz2"]
        ).getMolecule(0)
        mol = decouple(mol)
        return mol.toSystem()

    def test_generate_fep_masks(self, system):
        config = ConfigFactory(system, FreeEnergyMinimisation())
        res = config._generate_amber_fep_masks(0.001)
        expected_res = {
            "noshakemask": '"@1-45"',
            "scmask1": '"@1-45"',
            "scmask2": '""',
            "timask1": '"@1-45"',
            "timask2": '""',
        }
        for key in expected_res:
            assert expected_res[key] == res[key]
