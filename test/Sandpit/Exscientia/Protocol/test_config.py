import pandas as pd
import pytest

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Protocol import FreeEnergyMinimisation
from BioSimSpace.Sandpit.Exscientia.Align._decouple import decouple

class TestGromacsRBFE():
    @staticmethod
    @pytest.fixture(scope='class')
    def system():
        m0 = BSS.IO.readMolecules("test/input/ligands/CAT-13a*").getMolecule(0)
        m1 = BSS.IO.readMolecules("test/input/ligands/CAT-13c*").getMolecule(0)
        atom_mapping = BSS.Align.matchAtoms(m0, m1)
        m0 = BSS.Align.rmsdAlign(m0, m1, atom_mapping)
        merged = BSS.Align.merge(m0, m1)
        return merged.toSystem()

    def test_fep(self, system):
        '''Test if the default config writer will write the expected thing'''
        protocol = FreeEnergyMinimisation(lam=0.0,
                 lam_vals=None,
                 min_lam=0.0,
                 max_lam=1.0,
                 num_lam=11,
                 steps=10000,
                 perturbation_type="full")
        freenrg = BSS.FreeEnergy.Relative(system, protocol, engine='GROMACS', )
        with open(f"{freenrg._work_dir}/lambda_6/gromacs.mdp", 'r') as f:
            mdp_text = f.read()
            assert 'fep-lambdas          = 0.00000 0.10000 0.20000 0.30000 0.40000 0.50000 0.60000 0.70000 0.80000 0.90000 1.00000' in mdp_text
            assert 'init-lambda-state = 6' in mdp_text

    def test_fep_df(self, system):
        '''Test if the default config writer configured with the pd.DataFrame
        will write the expected thing.'''
        protocol = FreeEnergyMinimisation(lam=pd.Series(data={'fep': 0.0}),
                 lam_vals=None,
                 min_lam=pd.Series(data={'fep': 0.0}),
                 max_lam=pd.Series(data={'fep': 1.0}),
                 num_lam=11,
                 steps=10000,
                 perturbation_type="full")
        freenrg = BSS.FreeEnergy.Relative(system, protocol, engine='GROMACS', )
        with open(f"{freenrg._work_dir}/lambda_6/gromacs.mdp", 'r') as f:
            mdp_text = f.read()
            assert 'fep-lambdas          = 0.00000 0.10000 0.20000 0.30000 0.40000 0.50000 0.60000 0.70000 0.80000 0.90000 1.00000' in mdp_text
            assert 'init-lambda-state = 6' in mdp_text

    def test_staged_fep_df(self, system):
        '''Test if the multi-stage lambda will be written correctly'''
        protocol = FreeEnergyMinimisation(
            lam=pd.Series(data={'bonded': 0.0, 'coul': 0.0, 'vdw': 0.0}),
            lam_vals=pd.DataFrame(
                data={'bonded': [0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0],
                      'coul': [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
                      'vdw': [0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0]}),
            steps=10000, perturbation_type="full")
        freenrg = BSS.FreeEnergy.Relative(system, protocol, engine='GROMACS', )
        with open(f"{freenrg._work_dir}/lambda_6/gromacs.mdp", 'r') as f:
            mdp_text = f.read()
            assert 'bonded-lambdas       = 0.00000 0.25000 0.50000 0.75000 1.00000 1.00000 1.00000' in mdp_text
            assert 'coul-lambdas         = 0.00000 0.00000 0.00000 0.50000 1.00000 1.00000 1.00000' in mdp_text
            assert 'vdw-lambdas          = 0.00000 0.00000 0.00000 0.00000 0.00000 0.50000 1.00000' in mdp_text
            assert 'init-lambda-state = 6' in mdp_text

class TestGromacsABFE():
    @staticmethod
    @pytest.fixture(scope='class')
    def system():
        # Benzene.
        m0 = BSS.Parameters.openff_unconstrained_2_0_0(
            "c1ccccc1").getMolecule()
        protocol = FreeEnergyMinimisation(lam=0.0,
                 lam_vals=None,
                 min_lam=0.0,
                 max_lam=1.0,
                 num_lam=11,
                 steps=10000,
                 perturbation_type="full")
        return m0, protocol

    def test_decouple_vdw_q(self, system):
        m, protocol = system
        '''Test the decoupling where lambda0 = vdw-q and lambda1=none'''
        mol = decouple(m)
        freenrg = BSS.FreeEnergy.Relative(mol.toSystem(), protocol, engine='GROMACS', )
        with open(f"{freenrg._work_dir}/lambda_6/gromacs.mdp", 'r') as f:
            mdp_text = f.read()
            assert 'couple-moltype = c1ccccc1' in mdp_text
            assert 'couple-lambda0 = vdw-q' in mdp_text
            assert 'couple-lambda1 = none' in mdp_text
            assert 'couple-intramol = yes' in mdp_text

    def test_annihilate_vdw2q(self, system):
        '''Test the annihilation where lambda0 = vdw-q and lambda1=none'''
        m, protocol = system
        mol = decouple(m, charge=(False, True), LJ=(True, False),
                       intramol=False)
        freenrg = BSS.FreeEnergy.Relative(mol.toSystem(), protocol, engine='GROMACS', )
        with open(f"{freenrg._work_dir}/lambda_6/gromacs.mdp", 'r') as f:
            mdp_text = f.read()
            assert 'couple-moltype = c1ccccc1' in mdp_text
            assert 'couple-lambda0 = vdw' in mdp_text
            assert 'couple-lambda1 = q' in mdp_text
            assert 'couple-intramol = no' in mdp_text

    def test_sc_parameters(self, system):
        '''Test if the soft core parameters have been written.
        The default sc-alpha is 0, which means the soft-core of the vdw is not
        turned on by default. This checks if the this value has been changed to
        0.5.'''
        m, protocol = system
        mol = decouple(m)
        freenrg = BSS.FreeEnergy.Relative(mol.toSystem(), protocol, engine='GROMACS', )

        with open(f"{freenrg._work_dir}/lambda_6/gromacs.mdp", 'r') as f:
            mdp_text = f.read()
            assert 'sc-alpha = 0.5' in mdp_text
