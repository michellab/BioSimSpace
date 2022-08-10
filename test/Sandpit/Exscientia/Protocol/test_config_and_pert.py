import pandas as pd
import pytest

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Protocol import FreeEnergyMinimisation, FreeEnergy
from BioSimSpace.Sandpit.Exscientia.Align._decouple import decouple
from BioSimSpace.Sandpit.Exscientia.Units.Length import angstrom
from BioSimSpace.Sandpit.Exscientia.Units.Angle import radian, degree
from BioSimSpace.Sandpit.Exscientia.Units.Energy import kcal_per_mol
from BioSimSpace.Sandpit.Exscientia.Units.Temperature import kelvin
from BioSimSpace.Sandpit.Exscientia.FreeEnergy import Restraint


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

class TestSomdABFE():
    @staticmethod
    @pytest.fixture(scope='class')
    def system_and_restraint():
        # Benzene.
        m = BSS.Parameters.openff_unconstrained_2_0_0(
                    "c1ccccc1").getMolecule()

        # Assign atoms for restraint
        atom_1 = m.getAtoms()[0]
        atom_2 = m.getAtoms()[1]
        atom_3 = m.getAtoms()[2]
        atom_4 = m.getAtoms()[3]
        atom_5 = m.getAtoms()[4]
        atom_6 = m.getAtoms()[5]

        mol = decouple(m)
        system = mol.toSystem()

        # Create random restraint dictionary
        restraint_dict = {
                "anchor_points":{"r1":atom_1, "r2":atom_2, "r3":atom_3,
                                "l1":atom_4, "l2":atom_5, "l3":atom_6},
                "equilibrium_values":{"r0": 5.08 * angstrom,
                                    "thetaA0": 64.051 * degree,
                                    "thetaB0": 39.618 * degree,
                                    "phiA0":2.59 * radian,
                                    "phiB0":-1.20 * radian,
                                    "phiC0":2.63 * radian},
                "force_constants":{"kr":10 * kcal_per_mol / angstrom ** 2,
                                "kthetaA":10 * kcal_per_mol / (radian * radian),
                                "kthetaB":10 * kcal_per_mol / (radian * radian),
                                "kphiA":10 * kcal_per_mol / (radian * radian),
                                "kphiB":10 * kcal_per_mol / (radian * radian),
                                "kphiC":10 * kcal_per_mol / (radian * radian)}}

        restraint = Restraint(system, restraint_dict, 298 * kelvin, rest_type='Boresch')

        return system, restraint

    def test_turn_on_restraint(self, system_and_restraint):
        '''Test for turning on the restraint'''
        system, restraint = system_and_restraint
        protocol = FreeEnergy(perturbation_type="restraint")
        freenrg = BSS.FreeEnergy.Absolute(system, protocol, engine='SOMD', restraint=restraint)

        # Test .cfg file
        with open(f"{freenrg._work_dir}/lambda_0.000/somd.cfg", 'r') as f:
            cfg_text = f.read()
            assert 'use boresch restraints = True' in cfg_text
            assert 'boresch restraints dictionary = {"anchor_points":{"r1":1,'
            ' "r2":2, "r3":3, "l1":4, "l2":5, "l3":6}, "equilibrium_values"'
            ':{"r0":5.08, "thetaA0":1.12, "thetaB0":0.69,"phiA0":2.59, "phiB0"'
            ':-1.20, "phiC0":2.63}, "force_constants":{"kr":5.00, "kthetaA":5.00,'
            ' "kthetaB":5.00, "kphiA":5.00, "kphiB":5.00, "kphiC":5.00}}' in cfg_text
            assert 'turn on receptor-ligand restraints mode = True' in cfg_text

        # Test .pert file
        with open(f"{freenrg._work_dir}/lambda_0.000/somd.pert", 'r') as f:
            pert_text = f.read()
            # Check all atoms present
            carbons = [f"C{i}x" for i in range(1, 7)]
            hydrogens = [f"H{i}x" for i in range(1, 7)]
            atoms = carbons + hydrogens
            for atom in atoms:
                assert atom in pert_text
            # Check perturbations are correct
            lines = ["initial_type   C1",
                    "final_type     C1",
                    "initial_LJ     3.48065 0.08688",
                    "final_LJ       3.48065 0.08688",
                    "initial_charge -0.13000",
                    "final_charge   -0.13000"]
            for line in lines:
                assert line in pert_text
            
    def test_discharge(self, system_and_restraint):
        '''Test for discharging the ligand'''
        system, restraint = system_and_restraint
        protocol = FreeEnergy(perturbation_type="discharge_soft")
        freenrg = BSS.FreeEnergy.Absolute(system, protocol, engine='SOMD', restraint=restraint)

        # Test .cfg file
        with open(f"{freenrg._work_dir}/lambda_0.000/somd.cfg", 'r') as f:
            cfg_text = f.read()
            assert 'use boresch restraints = True' in cfg_text
            assert 'boresch restraints dictionary = {"anchor_points":{"r1":1,'
            ' "r2":2, "r3":3, "l1":4, "l2":5, "l3":6}, "equilibrium_values"'
            ':{"r0":5.08, "thetaA0":1.12, "thetaB0":0.69,"phiA0":2.59, "phiB0"'
            ':-1.20, "phiC0":2.63}, "force_constants":{"kr":5.00, "kthetaA":5.00,'
            ' "kthetaB":5.00, "kphiA":5.00, "kphiB":5.00, "kphiC":5.00}}' in cfg_text
            assert 'turn on receptor-ligand restraints mode = True' not in cfg_text

        # Test .pert file
        with open(f"{freenrg._work_dir}/lambda_0.000/somd.pert", 'r') as f:
            pert_text = f.read()
            # Check all atoms present
            carbons = [f"C{i}x" for i in range(1, 7)]
            hydrogens = [f"H{i}x" for i in range(1, 7)]
            atoms = carbons + hydrogens
            for atom in atoms:
                assert atom in pert_text
            # Check perturbations are correct
            lines = ["initial_type   C1",
                    "final_type     C1",
                    "initial_LJ     3.48065 0.08688",
                    "final_LJ       3.48065 0.08688",
                    "initial_charge -0.13000",
                    "final_charge   -0.00000"]
            for line in lines:
                assert line in pert_text

    def test_vanish(self, system_and_restraint):
        '''Test for vanishing the ligand'''
        system, restraint = system_and_restraint
        protocol = FreeEnergy(perturbation_type="vanish_soft")
        freenrg = BSS.FreeEnergy.Absolute(system, protocol, engine='SOMD', restraint=restraint)

        # Test .cfg file
        with open(f"{freenrg._work_dir}/lambda_0.000/somd.cfg", 'r') as f:
            cfg_text = f.read()
            assert 'use boresch restraints = True' in cfg_text
            assert 'boresch restraints dictionary = {"anchor_points":{"r1":1,'
            ' "r2":2, "r3":3, "l1":4, "l2":5, "l3":6}, "equilibrium_values"'
            ':{"r0":5.08, "thetaA0":1.12, "thetaB0":0.69,"phiA0":2.59, "phiB0"'
            ':-1.20, "phiC0":2.63}, "force_constants":{"kr":5.00, "kthetaA":5.00,'
            ' "kthetaB":5.00, "kphiA":5.00, "kphiB":5.00, "kphiC":5.00}}' in cfg_text
            assert 'turn on receptor-ligand restraints mode = True' not in cfg_text

        # Test .pert file
        with open(f"{freenrg._work_dir}/lambda_0.000/somd.pert", 'r') as f:
            pert_text = f.read()
            # Check all atoms present
            carbons = [f"C{i}x" for i in range(1, 7)]
            hydrogens = [f"H{i}x" for i in range(1, 7)]
            atoms = carbons + hydrogens
            for atom in atoms:
                assert atom in pert_text
            # Check perturbations are correct
            lines = ["initial_type   C1",
                    "final_type     du",
                    "initial_LJ     3.48065 0.08688",
                    "final_LJ       0.00000 0.00000",
                    "initial_charge -0.00000",
                    "final_charge   -0.00000"]
            for line in lines:
                assert line in pert_text

    def test_discharge_and_vanish(self, system_and_restraint):
        '''Test for simultaneously discharging and vanishing the ligand'''
        system, restraint = system_and_restraint
        protocol = FreeEnergy(perturbation_type="full")
        freenrg = BSS.FreeEnergy.Absolute(system, protocol, engine='SOMD', restraint=restraint)

        # Test .cfg file
        with open(f"{freenrg._work_dir}/lambda_0.000/somd.cfg", 'r') as f:
            cfg_text = f.read()
            assert 'use boresch restraints = True' in cfg_text
            assert 'boresch restraints dictionary = {"anchor_points":{"r1":1,'
            ' "r2":2, "r3":3, "l1":4, "l2":5, "l3":6}, "equilibrium_values"'
            ':{"r0":5.08, "thetaA0":1.12, "thetaB0":0.69,"phiA0":2.59, "phiB0"'
            ':-1.20, "phiC0":2.63}, "force_constants":{"kr":5.00, "kthetaA":5.00,'
            ' "kthetaB":5.00, "kphiA":5.00, "kphiB":5.00, "kphiC":5.00}}' in cfg_text
            assert 'turn on receptor-ligand restraints mode = True' not in cfg_text

        # Test .pert file
        with open(f"{freenrg._work_dir}/lambda_0.000/somd.pert", 'r') as f:
            pert_text = f.read()
            # Check all atoms present
            carbons = [f"C{i}x" for i in range(1, 7)]
            hydrogens = [f"H{i}x" for i in range(1, 7)]
            atoms = carbons + hydrogens
            for atom in atoms:
                assert atom in pert_text
            # Check perturbations are correct
            lines = ["initial_type   C1",
                    "final_type     du",
                    "initial_LJ     3.48065 0.08688",
                    "final_LJ       0.00000 0.00000",
                    "initial_charge -0.13000",
                    "final_charge   -0.00000"]
            for line in lines:
                assert line in pert_text