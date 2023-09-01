from difflib import unified_diff

import itertools
import os

import pandas as pd
import pytest
from sire.legacy.IO import AmberRst

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Units.Energy import kj_per_mol
from BioSimSpace.Sandpit.Exscientia.Units.Length import angstrom

from tests.Sandpit.Exscientia.conftest import has_amber, has_gromacs, has_openff


@pytest.fixture
def system():
    ff = "openff_unconstrained-2.0.0"
    mol0 = BSS.Parameters.parameterise("c1ccccc1C", ff).getMolecule()
    mol1 = BSS.Parameters.parameterise("c1ccccc1", ff).getMolecule()
    return BSS.Align.merge(mol0, mol1).toSystem()


@pytest.fixture
def ref_system(system):
    mol = system[0]
    mol.translate(3 * [BSS.Units.Length.nanometer * 1])
    return mol.toSystem()


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


@pytest.mark.skipif(
    has_gromacs is False or has_openff is False,
    reason="Requires GROMACS and openff to be installed",
)
def test_gromacs(protocol, system, ref_system, tmp_path):
    proc = BSS.Process.Gromacs(
        system,
        protocol,
        reference_system=ref_system,
        work_dir=str(tmp_path),
        ignore_warnings=True,
    )

    # We have the restraints
    assert (tmp_path / "posre_0001.itp").is_file()
    with open(tmp_path / "gromacs.top", "r") as f:
        assert "posre_0001.itp" in f.read()

    # We have generated a separate restraint reference
    assert os.path.exists(proc._ref_file)
    contents_ref, contents_gro = (
        open(proc._ref_file).readlines(),
        open(proc._gro_file).readlines(),
    )
    diff = list(unified_diff(contents_ref, contents_gro))
    assert len(diff)


@pytest.mark.skipif(
    has_openff is False,
    reason="Requires openff to be installed",
)
def test_amber(protocol, system, ref_system, tmp_path):
    proc = BSS.Process.Amber(
        system, protocol, reference_system=ref_system, work_dir=str(tmp_path)
    )

    # We have the restraints
    with open(tmp_path / "amber.cfg", "r") as f:
        cfg = f.read()
        assert "restraint_wt" in cfg
        assert "restraintmask" in cfg

    # We have generated a separate restraint reference
    assert os.path.exists(proc._ref_file)

    ref = AmberRst(proc._ref_file).getFrame(0)
    rst = AmberRst(proc._rst_file).getFrame(0)

    assert ref != rst

    # We are pointing the reference to the correct file
    assert f"{proc._work_dir}/{proc.getArgs()['-ref']}" == proc._ref_file
