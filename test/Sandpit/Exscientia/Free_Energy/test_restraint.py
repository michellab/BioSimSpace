import pytest

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Align._decouple import decouple
from BioSimSpace.Sandpit.Exscientia.Protocol import Equilibration, Production, Minimisation
from BioSimSpace.Sandpit.Exscientia.FreeEnergy import Restraint
# from BioSimSpace.Sandpit.Exscientia import Process as _Process
from BioSimSpace.Sandpit.Exscientia.Units.Length import angstrom
from BioSimSpace.Sandpit.Exscientia.Units.Angle import radian
from BioSimSpace.Sandpit.Exscientia.Units.Energy import kcal_per_mol

@pytest.fixture(scope='session')
def restraint():
    '''Generate the restaint object with _Boresch already filled.'''
    ligand = BSS.IO.readMolecules('../input/ligands/ligand01*').getMolecule(0)
    decoupled_ligand = decouple(ligand)

    protein = BSS.IO.readMolecules(['../input/molecules/1jr5.crd',
                                    '../input/molecules/1jr5.top']).getMolecule(0)

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
             4 * kcal_per_mol / radian ** 2),
            ((ligand_0, protein_1, protein_2), 5 * radian,
             6 * kcal_per_mol / radian ** 2),
        ],
        'Dihedrals': [
            ((ligand_2, ligand_1, ligand_0, protein_0), 7 * radian,
             8 * kcal_per_mol / radian ** 2),
            ((ligand_1, ligand_0, protein_0, protein_1), 9 * radian,
             10 * kcal_per_mol / radian ** 2),
            ((ligand_0, protein_0, protein_1, protein_2), 11 * radian,
             12 * kcal_per_mol / radian ** 2),
        ]}
    return restraint
