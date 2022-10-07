import BioSimSpace.Sandpit.Exscientia as BSS

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
def test_disulphide():
    """Test parameterisation in the presence of disulphide bridges."""

    # Load the example molecule.
    molecule = BSS.IO.readMolecules("test/Sandpit/Exscientia/input/molecules/4LYT_Fixed.pdb")[0]

    # Try to parameterise with ff14SB. If working, this should auto-detect
    # disulphide bonds and add the appropriate bond records to the tLEaP input
    # script.
    molecule = BSS.Parameters.ff14SB(molecule).getMolecule()

    # Check that we actually generate records for four disulphide bonds.
    bonds = BSS.Parameters.Protocol._protocol._get_disulphide_bonds(molecule._sire_object)
    assert len(bonds) == 4
