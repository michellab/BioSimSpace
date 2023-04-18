import pytest
from BioSimSpace.Sandpit.Exscientia.Align._merge import _removeDummies

from test.Sandpit.Exscientia.conftest import get_energy


def test__removeDummies(
    benzene, pyrrole, mapping_benzene_pyrrole, merged_benzene_pyrrole
):
    benzene_merged = _removeDummies(merged_benzene_pyrrole, is_lambda1=False)
    pyrrole_merged = _removeDummies(merged_benzene_pyrrole, is_lambda1=True)
    energy0_unmerged = get_energy(benzene)
    energy1_unmerged = get_energy(pyrrole)
    energy0_merged = get_energy(benzene_merged)
    energy1_merged = get_energy(
        pyrrole_merged, coord_mol=pyrrole, coord_mol_mapping=mapping_benzene_pyrrole
    )
    assert energy0_unmerged == pytest.approx(energy0_merged)
    assert energy1_unmerged == pytest.approx(energy1_merged)
