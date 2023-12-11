import os
import pytest
import tempfile

import BioSimSpace.Sandpit.Exscientia as BSS

from tests.Sandpit.Exscientia.conftest import (
    url,
    has_openff,
    has_tleap,
    has_antechamber,
)


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


@pytest.mark.skipif(
    has_antechamber is False or has_tleap is False,
    reason="Requires AmberTools/antechamber and tLEaP to be installed.",
)
def test_custom_parameters():
    """
    Test that custom parameters are correctly inserted into the LEaP input
    script.
    """

    # First parameterise a molecule with the GAFF force field.
    mol = BSS.Parameters.gaff("C").getMolecule()

    # Now try to parameterise with a protein force field, also including the
    # GAFF parameters and a custom parameter file. Note that we only validate
    # that the file exists and LEaP won't error if it doesn't contain relevant
    # parameters.

    # The path to this file.
    file_path = os.path.realpath(__file__)

    with tempfile.TemporaryDirectory() as tmp:
        # This will fail, but we only care about the LEaP input script.
        try:
            mol = BSS.Parameters.ff14SB(
                mol,
                work_dir=tmp,
                custom_parameters=["leaprc.gaff", file_path],
            ).getMolecule()
        except:
            pass

        # Load the script and check that the custom parameters are present and
        # in the correct order.
        with open(f"{tmp}/leap.txt", "r") as f:
            script = f.readlines()
            line_gaff = 0
            line_file = 0
            line_mol = 0
            for x, line in enumerate(script):
                line = line.strip()
                if line == "source leaprc.gaff":
                    line_gaff = x
                elif line == f"loadAmberParams {file_path}":
                    line_file = x
                elif line == "mol = loadPdb leap.pdb":
                    line_mol = x

            # Make sure the lines are found.
            assert line_gaff != 0
            assert line_file != 0
            assert line_mol != 0

            # Make sure the lines are in the correct order.
            assert line_gaff < line_file < line_mol
