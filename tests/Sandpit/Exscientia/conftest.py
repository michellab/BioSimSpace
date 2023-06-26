from typing import Dict, Optional

import pytest
from BioSimSpace.Sandpit import Exscientia as BSS

from BioSimSpace.Sandpit.Exscientia.Process import Gromacs
from BioSimSpace.Sandpit.Exscientia.Protocol import Production
from BioSimSpace.Sandpit.Exscientia.Units.Temperature import kelvin
from BioSimSpace.Sandpit.Exscientia.Units.Time import femtosecond
from BioSimSpace.Sandpit.Exscientia._SireWrappers import Molecule


def get_energy(
    mol: Molecule,
    coord_mol: Optional[Molecule] = None,
    coord_mol_mapping: Optional[Dict[int, int]] = None,
) -> float:
    if coord_mol is not None:
        c = mol._sire_object.cursor()
        c_coord = coord_mol._sire_object.cursor()
        coordinates = c_coord["coordinates"].toVector()
        if coord_mol_mapping is not None:
            sorted_keys = sorted(coord_mol_mapping)
            coordinates = [coordinates[coord_mol_mapping[k]] for k in sorted_keys]
        for atom, coord in zip(c.atoms(), coordinates):
            atom["coordinates"] = coord
        mol = Molecule(c.commit())

    system = mol.toSystem()

    protocol = Production(
        timestep=0.1 * femtosecond,
        runtime=0.1 * femtosecond,
        temperature=0 * kelvin,
    )

    proc = Gromacs(system, protocol)
    proc.start()
    proc.wait()

    res = proc.getPotentialEnergy().value()

    return res


@pytest.fixture
def benzene():
    return BSS.Parameters.parameterise(
        "c1ccccc1", "openff_unconstrained-2.0.0"
    ).getMolecule()


@pytest.fixture
def pyrrole():
    return BSS.Parameters.parameterise(
        "C1=CNC=C1", "openff_unconstrained-2.0.0"
    ).getMolecule()


@pytest.fixture
def mapping_benzene_pyrrole(benzene, pyrrole):
    return BSS.Align.matchAtoms(benzene, pyrrole, complete_rings_only=False)


@pytest.fixture
def merged_benzene_pyrrole(benzene, pyrrole, mapping_benzene_pyrrole):
    return BSS.Align.merge(
        benzene,
        pyrrole,
        mapping=mapping_benzene_pyrrole,
        allow_ring_breaking=True,
        allow_ring_size_change=True,
    )
