import BioSimSpace as BSS

import os
import pytest

# Make sure required AMBER executables are present.
if BSS._amber_home is not None:
    tleap = "%s/bin/tleap" % BSS._amber_home
    if os.path.isfile(tleap):
        has_tleap = True
    else:
        has_tleap = False
else:
    has_tleap = False


@pytest.mark.skipif(has_tleap is False, reason="Requires tLEaP to be installed.")
@pytest.mark.parametrize("ff", BSS.Parameters.amberProteinForceFields())
def test_disulphide(ff):
    """Test parameterisation in the presence of disulphide bridges."""

    # Load the example molecule.
    molecule = BSS.IO.readMolecules("test/input/molecules/4LYT_Fixed.pdb")[0]

    # Try to parameterise with the named force field. If working, this should
    # auto-detect disulphide bonds and add the appropriate bond records to the
    # tLEaP input script.
    molecule = getattr(BSS.Parameters, ff)(molecule).getMolecule()

    # Check that we actually generate records for four disulphide bonds.
    bonds = BSS.Parameters._Protocol.AmberProtein._get_disulphide_bonds(
        molecule._sire_object
    )
    assert len(bonds) == 4

    # Check that the bond parameters are present in the molecule.
    bonds = molecule.search("bonds from element S to element S")
    assert len(bonds) == 4
