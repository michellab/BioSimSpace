import itertools
import pytest

import pandas as pd

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Units.Length import angstrom
from BioSimSpace.Sandpit.Exscientia.Units.Energy import kj_per_mol


# Make sure GROMSCS is installed.
has_gromacs = BSS._gmx_exe is not None


@pytest.fixture
def system():
    ff = "openff_unconstrained-2.0.0"
    mol0 = BSS.Parameters.parameterise("c1ccccc1C", ff).getMolecule()
    mol1 = BSS.Parameters.parameterise("c1ccccc1", ff).getMolecule()
    return BSS.Align.merge(mol0, mol1).toSystem()


@pytest.fixture
def restraint():
    return {
        "restraint": "heavy",
        "force_constant": 1000 * kj_per_mol / angstrom**2,
    }


@pytest.fixture
def free_energy():
    return {
        "lam": pd.Series(data={"fep": 0.5}),
        "lam_vals": pd.DataFrame(
            data={"fep": [0.5, 0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0]}
        ),
    }


@pytest.fixture()
def minimisation():
    return {"steps": 100}


@pytest.fixture()
def equilibration(production):
    return {
        **production,
        "temperature_start": BSS.Types.Temperature(300, "kelvin"),
        "temperature_end": BSS.Types.Temperature(310, "kelvin"),
    }


@pytest.fixture()
def production():
    return {"runtime": BSS.Types.Time(0.001, "nanoseconds")}


@pytest.fixture(
    params=itertools.product(
        [True, False], ["Minimisation", "Equilibration", "Production"]
    )
)
def protocol(request, restraint, free_energy, minimisation, equilibration, production):
    FE, name = request.param
    if FE:
        if name == "Minimisation":
            return BSS.Protocol.FreeEnergyMinimisation(
                **minimisation, **restraint, **free_energy
            )
        elif name == "Equilibration":
            return BSS.Protocol.FreeEnergyEquilibration(
                **equilibration, **restraint, **free_energy
            )
        else:
            return BSS.Protocol.FreeEnergy(**production, **restraint, **free_energy)
    else:
        if name == "Minimisation":
            return BSS.Protocol.Minimisation(**minimisation, **restraint)
        elif name == "Equilibration":
            return BSS.Protocol.Equilibration(**equilibration, **restraint)
        else:
            return BSS.Protocol.Production(**production, **restraint)


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_gromacs(protocol, system, tmp_path):
    BSS.Process.Gromacs(system, protocol, work_dir=str(tmp_path), ignore_warnings=True)
    assert (tmp_path / "posre_0001.itp").is_file()
    with open(tmp_path / "gromacs.top", "r") as f:
        assert "posre_0001.itp" in f.read()


def test_amber(protocol, system, tmp_path):
    BSS.Process.Amber(system, protocol, work_dir=str(tmp_path))
    with open(tmp_path / "amber.cfg", "r") as f:
        cfg = f.read()
        assert "restraint_wt" in cfg
        assert "restraintmask" in cfg
