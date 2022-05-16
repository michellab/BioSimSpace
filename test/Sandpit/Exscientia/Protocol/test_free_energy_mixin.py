import pandas as pd
import pytest
import numpy as np

from BioSimSpace.Sandpit.Exscientia.Protocol import _FreeEnergyMixin

class Test_List_linspace():
    '''Test if the new _FreeEnergyMixin class could handle the num_lam as input
    with input as float'''
    @staticmethod
    @pytest.fixture(scope='class')
    def protocol():
        protocol = _FreeEnergyMixin(lam=0.5, lam_vals=None, min_lam=0.0,
                                    max_lam=1.0, num_lam=11)
        return protocol

    def test_getLambdaValues(self, protocol):
        value = protocol.getLambdaValues(type='list')
        assert isinstance(value, list)
        assert isinstance(value[5], float)
        assert np.isclose(value[5], 0.5)

    def test_getLambda(self, protocol):
        value = protocol.getLambda(type='float')
        assert isinstance(value, float)
        assert np.isclose(value, 0.5)

    def test_getLambdaIndex(self, protocol):
        value = protocol.getLambdaIndex()
        assert isinstance(value, int)
        assert np.isclose(value, 5)

    def test_getLambdaValues_df(self, protocol):
        value = protocol.getLambdaValues(type='dataframe')
        assert isinstance(value, pd.DataFrame)
        assert isinstance(value.loc[5], pd.Series)
        assert np.isclose(value.loc[5], 0.5)

    def test_getLambda_Series(self, protocol):
        value = protocol.getLambda(type='series')
        assert isinstance(value, pd.Series)
        assert np.isclose(value, 0.5)

class Test_List_explict(Test_List_linspace):
    '''Test if the new _FreeEnergyMixin class could handle the lam_vals as input
    with input as float'''
    @staticmethod
    @pytest.fixture(scope='class')
    def protocol():
        protocol = _FreeEnergyMixin(lam=0.5,
                                    # Note that 0.5 is the first to check if being sorted
                                    lam_vals=[0.5, 0.0, 0.1, 0.2, 0.3, 0.4,
                                              0.6, 0.7, 0.8, 0.9, 1.0], )
        return protocol

class Test_List_explict_df(Test_List_explict):
    '''Test if the new _FreeEnergyMixin class could handle the lam_vals as input
    with input as pd.DataFrame'''
    @staticmethod
    @pytest.fixture(scope='class')
    def protocol():
        protocol = _FreeEnergyMixin(lam=pd.Series(data={'fep': 0.5}),
                                    # Note that 0.5 is the first to check if being sorted
                                    lam_vals=pd.DataFrame(data={'fep':
                                        [0.5, 0.0, 0.1, 0.2, 0.3, 0.4, 0.6,
                                         0.7, 0.8, 0.9, 1.0]}))
        return protocol

class Test_List_linspace_df():
    '''Test if the new _FreeEnergyMixin class could handle the num_lam as input
    with input as pd.DataFrame'''
    @staticmethod
    @pytest.fixture(scope='class')
    def protocol():
        protocol = _FreeEnergyMixin(lam=pd.Series(data={'fep': 0.5}),
                                    lam_vals=None,
                                    min_lam=pd.Series(data={'fep': 0.0}),
                                    max_lam=pd.Series(data={'fep': 1.0}),
                                    num_lam=11)
        return protocol

def test_df_num_lam():
    '''Test use multiple column df as input to num_lam'''
    protocol = _FreeEnergyMixin(lam=pd.Series(data={'bonded': 0.5,
                                                    'coul': 0.5,
                                                    'vdw': 0.5}),
                                lam_vals=None,
                                min_lam=pd.Series(data={'bonded': 0.0,
                                                        'coul': 0.0,
                                                        'vdw': 0.0}),
                                max_lam=pd.Series(data={'bonded': 1.0,
                                                        'coul': 1.0,
                                                        'vdw': 1.0}),
                                num_lam=11)
    assert protocol.getLambdaIndex() == 5

def test_df_lam_vals():
    '''Test use multiple column df as input to lam_vals'''
    protocol = _FreeEnergyMixin(
        lam=pd.Series(data={'bonded': 0.75, 'coul': 0.5, 'vdw': 0.0}),
        lam_vals=pd.DataFrame(data={'bonded': [0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0],
                                    'coul':   [0.0, 0.0,  0.0, 0.5,  1.0, 1.0, 1.0],
                                    'vdw':    [0.0, 0.0,  0.0, 0.0,  0.0, 0.5, 1.0]}))
    assert protocol.getLambdaIndex() == 3