import pathlib
from os.path import exists

import pytest
import pandas as pd
import numpy as np
try:
    from alchemlyb.parsing.gmx import extract_u_nk
    is_alchemlyb = True
except ModuleNotFoundError:
    is_alchemlyb = False

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Protocol import FreeEnergyEquilibration, FreeEnergy
from BioSimSpace.Sandpit.Exscientia.Align._decouple import decouple
from BioSimSpace.Sandpit.Exscientia import Types as _Types

# Make sure GROMSCS is installed.
has_gromacs = BSS._gmx_exe is not None

@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
class Test_gmx_ABFE():
    @staticmethod
    @pytest.fixture(scope='class')
    def freenrg():
        m = BSS.IO.readMolecules("test/input/amber/ala/*").getMolecule(0)
        decouple_m = decouple(m)
        solvated = BSS.Solvent.tip3p(molecule=decouple_m,
                                     box=3 * [3 * BSS.Units.Length.nanometer])
        protocol = FreeEnergyEquilibration(
            lam=pd.Series(data={'coul': 0.5, 'vdw': 0.0}),
            lam_vals=pd.DataFrame(
                data={'coul': [0.0, 0.5, 1.0, 1.0, 1.0],
                      'vdw':  [0.0, 0.0, 0.0, 0.5, 1.0]}),
            runtime=_Types.Time(0, "nanoseconds"),)
        freenrg = BSS.FreeEnergy.Relative(solvated, protocol,
                                          engine='GROMACS', )
        freenrg.run()
        freenrg.wait()
        return freenrg

    def test_file_exist(self, freenrg):
        '''Test if all the files are there.'''
        path = pathlib.Path(freenrg._work_dir)
        for i in range(5):
            assert (path / f'lambda_{i}' / 'gromacs.xvg').is_file()

    @pytest.mark.skipif(is_alchemlyb is False, reason='Need alchemlyb.')
    def test_lambda(self, freenrg):
        '''Test if the xvg files contain the correct lambda.'''
        path = pathlib.Path(freenrg._work_dir)
        for i, (coul, vdw) in enumerate(zip([0.0, 0.5, 1.0, 1.0, 1.0],
                                            [0.0, 0.0, 0.0, 0.5, 1.0])):
            u_nk = extract_u_nk(path / f'lambda_{i}' / 'gromacs.xvg', 300)
            assert u_nk.index.names == ['time', 'coul-lambda', 'vdw-lambda']
            assert np.isclose(u_nk.index.values[0][1], coul)
            assert np.isclose(u_nk.index.values[0][2], vdw)

class Test_Somd_ABFE():
    @staticmethod
    @pytest.fixture(scope='class')
    def freenrg():
        m = BSS.Parameters.openff_unconstrained_2_0_0(
                    "c1ccccc1").getMolecule()
        mol = decouple(m)
        system = mol.toSystem()

        decouple_m = decouple(m)
        solvated = BSS.Solvent.tip3p(molecule=decouple_m,
                                     box=3 * [3 * BSS.Units.Length.nanometer])
        protocol = BSS.Protocol.FreeEnergy(lam_vals = [0.0, 0.5, 1.0],
                              runtime=0.0001*BSS.Units.Time.nanosecond)
        freenrg = BSS.FreeEnergy.Relative(solvated, protocol,
                                          engine='SOMD')

        freenrg.run()
        freenrg.wait()
        return freenrg

    def test_files_exist(self, freenrg):
        '''Test if all the files have been created.'''
        path = pathlib.Path(freenrg._work_dir)
        # TODO: Fix _process_runner.wait so this works with range(3)
        for i in range(1):
            assert (path / f'lambda_{i}' / 'simfile.dat').is_file()
            assert (path / f'lambda_{i}' / 'sim_restart.s3').is_file()
            assert (path / f'lambda_{i}' / 'traj000000001.dcd').is_file()
            assert (path / f'lambda_{i}' / 'gradients.dat').is_file()

    def test_correct_simfile(self, freenrg):
        '''Check if correct numbers of entries in simfile.dat'''
        path = pathlib.Path(freenrg._work_dir)
        with open(path / 'lambda_0' / 'simfile.dat', 'rt') as f:
            assert "#Alchemical array is		 (0.0, 0.5, 1.0)" in f.read()