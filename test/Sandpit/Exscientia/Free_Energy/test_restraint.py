import pytest

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Align._decouple import decouple
from BioSimSpace.Sandpit.Exscientia.Protocol import Production
from BioSimSpace.Sandpit.Exscientia.FreeEnergy import Restraint
# from BioSimSpace.Sandpit.Exscientia import Process as _Process
from BioSimSpace.Sandpit.Exscientia.Units.Length import angstrom
from BioSimSpace.Sandpit.Exscientia.Units.Angle import radian
from BioSimSpace.Sandpit.Exscientia.Units.Energy import kcal_per_mol
from BioSimSpace.Sandpit.Exscientia._SireWrappers import Atom
from BioSimSpace.Sandpit.Exscientia.Types import Length, Angle, Energy

@pytest.fixture(scope='session')
def restraint():
    '''Generate the restaint object with _Boresch already filled.'''
    ligand = BSS.IO.readMolecules(BSS.IO.glob("test/input/ligands/ligand01*")).getMolecule(0)
    decoupled_ligand = decouple(ligand)

    protein = BSS.IO.readMolecules(['test/input/molecules/1jr5.crd',
                                    'test/input/molecules/1jr5.top']).getMolecule(0)

    system = (protein + decoupled_ligand).toSystem()

    protocol = Production()

    restraint = Restraint(system, protocol=protocol, engine='GROMACS')

    # Assign three atoms from the ligand
    ligand_2 = decoupled_ligand.getAtoms()[2]
    ligand_1 = decoupled_ligand.getAtoms()[1]
    ligand_0 = decoupled_ligand.getAtoms()[0]

    # And three atoms from the protein
    protein_0 = protein.getAtoms()[0]
    protein_1 = protein.getAtoms()[1]
    protein_2 = protein.getAtoms()[2]

    # To form the Boresch with arbitrary values
    restraint._Boresch = {'Bonds': [
        ((ligand_0, protein_0), 1 * angstrom, 2 * kcal_per_mol / angstrom ** 2),
    ],
        'Angles': [
            ((ligand_1, ligand_0, protein_0), 3 * radian,
             4 * kcal_per_mol / (radian * radian)),
            ((ligand_0, protein_1, protein_2), 5 * radian,
             6 * kcal_per_mol / (radian * radian)),
        ],
        'Dihedrals': [
            ((ligand_2, ligand_1, ligand_0, protein_0), 7 * radian,
             8 * kcal_per_mol / (radian * radian)),
            ((ligand_1, ligand_0, protein_0, protein_1), 9 * radian,
             10 * kcal_per_mol / (radian * radian)),
            ((ligand_0, protein_0, protein_1, protein_2), 11 * radian,
             12 * kcal_per_mol / (radian * radian)),
        ]}
    return restraint

class TestBoreschType():
    '''Test the output type of the Boresch.'''
    @staticmethod
    @pytest.fixture(scope='class')
    def Boresch(restraint):
        return restraint._Boresch

    @pytest.mark.parametrize("name", ['Bonds', 'Angles', 'Dihedrals'])
    def test_type_atom(self, Boresch, name):
        'Check if the atoms are the type atom'
        fields = Boresch[name]
        for field in fields:
            for atom in field[0]:
                assert isinstance(atom, Atom)

    @pytest.mark.parametrize("name,length", [
        ('Bonds', 2), ('Angles', 3), ('Dihedrals', 4)])
    def test_length_atom(self, Boresch, name, length):
        'Check if the bond, angle, dihedral have the correct number of atoms.'
        fields = Boresch[name]
        for field in fields:
            assert len(field[0]) == length

    @pytest.mark.parametrize("name,unit", [
        ('Bonds', Length), ('Angles', Angle), ('Dihedrals', Angle)])
    def test_unit_value(self, Boresch, name, unit):
        '''Check if the the equilibrium value of bond, angle, dihedral has the
        correct unit.'''
        fields = Boresch[name]
        for field in fields:
            assert isinstance(field[1], unit)

    @pytest.mark.parametrize("name,angle,length", [
        ('Bonds', 0, 0), ('Angles', -2, 2), ('Dihedrals', -2, 2)])
    def test_unit_fc(self, Boresch, name, angle, length):
        '''Check if the the force constant of bond, angle, dihedral has the
        correct unit.'''
        fields = Boresch[name]
        for field in fields:
            assert field[2].dimensions() == (angle, 0, length, 1, -1, 0, -2)
