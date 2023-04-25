import BioSimSpace as BSS
from BioSimSpace._Utils import _try_import, _have_imported

import os
import pytest

# Make sure openff is installed.
_openff = _try_import("openff")
has_openff = _have_imported(_openff)

# Store the tutorial URL.
url = BSS.tutorialUrl()

# Make sure required AMBER executables are present.
if BSS._amber_home is not None:
    tleap = "%s/bin/tleap" % BSS._amber_home
    if os.path.isfile(tleap):
        has_tleap = True
    else:
        has_tleap = False
else:
    has_tleap = False

# Make sure antechamber is installed.
has_antechamber = BSS.Parameters._Protocol._amber._antechamber_exe is not None


@pytest.fixture(scope="session")
def molecule0():
    return BSS.IO.readMolecules(f"{url}/4LYT_Fixed.pdb.bz2")[0]


@pytest.fixture(scope="session")
def molecule1():
    return BSS.IO.readMolecules(f"{url}/3G8K_Fixed.pdb.bz2")[0]


@pytest.mark.skipif(has_tleap is False, reason="Requires tLEaP to be installed.")
@pytest.mark.parametrize("ff", BSS.Parameters.amberProteinForceFields())
def test_disulphide(molecule0, ff):
    """Test parameterisation in the presence of disulphide bridges."""

    # Try to parameterise with the named force field. If working, this should
    # auto-detect disulphide bonds and add the appropriate bond records to the
    # tLEaP input script.
    molecule = getattr(BSS.Parameters, ff)(molecule0).getMolecule()

    # Check that we actually generate records for four disulphide bonds.
    bonds = BSS.Parameters._Protocol.AmberProtein._get_disulphide_bonds(
        molecule0._sire_object
    )
    assert len(bonds) == 4

    # Check that the bond parameters are present in the molecule.
    bonds = molecule.search("bonds from element S to element S")
    assert len(bonds) == 4


@pytest.mark.skipif(has_tleap is False, reason="Requires tLEaP to be installed.")
def test_disulphide_renumber(molecule1, ff="ff14SB"):
    """
    Test parameterisation in the presence of disulphide bridges using a
    multi-chain PDB with duplicate residue numbering
    """

    # Try to parameterise with the named force field. If working, this should
    # auto-detect disulphide bonds and add the appropriate bond records to the
    # tLEaP input script.
    molecule = getattr(BSS.Parameters, ff)(molecule1).getMolecule()

    # Check that we actually generate records for eight disulphide bonds.
    bonds = BSS.Parameters._Protocol.AmberProtein._get_disulphide_bonds(
        molecule1._sire_object
    )
    assert len(bonds) == 8

    # Check that the bond parameters are present in the molecule.
    bonds = molecule.search("bonds from element S to element S")
    assert len(bonds) == 8


@pytest.mark.skipif(
    has_antechamber is False or has_openff is False,
    reason="Requires AmberTools/antechamber and OpenFF to be installed.",
)
def test_molecule_rename():
    """
    Test that a parameterised molecule generated from a SMILES string
    starting with the "[" character is renamed with a "smiles:" prefix
    so that it can be parsed by gmx when used in a GROMACS topology file.
    """

    # Create the parameterised molecule.
    mol = BSS.Parameters.openff_unconstrained_2_0_0(
        "[C@@H](C(F)(F)F)(OC(F)F)Cl"
    ).getMolecule()

    # Make sure the name is correct.
    assert mol._sire_object.name().value() == "smiles:[C@@H](C(F)(F)F)(OC(F)F)Cl"
