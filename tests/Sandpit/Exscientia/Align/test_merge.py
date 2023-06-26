import os
import pytest

import BioSimSpace.Sandpit.Exscientia as BSS

from BioSimSpace.Sandpit.Exscientia.Align._merge import _removeDummies
from BioSimSpace.Sandpit.Exscientia._Utils import _try_import, _have_imported

from tests.Sandpit.Exscientia.conftest import get_energy

# Check whether AMBER is installed.
if BSS._amber_home is not None:
    exe = "%s/bin/sander" % BSS._amber_home
    if os.path.isfile(exe):
        has_amber = True
    else:
        has_amber = False
else:
    has_amber = False

# Check whether GROMACS is installed.
has_gromacs = BSS._gmx_exe is not None

# Make sure openff is installed.
_openff = _try_import("openff")
has_openff = _have_imported(_openff)

try:
    import alchemlyb

    has_alchemlyb = True
except ModuleNotFoundError:
    has_alchemlyb = False


@pytest.mark.skipif(
    has_amber is False or has_gromacs is False or has_openff is False or has_alchemlyb,
    reason="Requires that AMBER, GROMACS, OpenFF, and alchemlyb are installed.",
)
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
