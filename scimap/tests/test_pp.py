#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 22:16:06 2024
@author: aj
Pre processing tools tests
"""

import pytest
import os
import anndata as ad
import pandas as pd
import numpy as np


#os.chdir ('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/softwares/scimap')
# np.savez('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/softwares/scimap/scimap/tests/expected_test_values/test_log1p.npz', data=adata.layers['log_test'])

# loading via mcmicro
def test_mcmicro_to_scimap():
    from scimap.preprocessing.mcmicro_to_scimap import mcmicro_to_scimap
    feature_table_path = [os.getcwd() + '/scimap/tests/scimapExampleData/exemplar-001--unmicst_cell.csv']
    adata = mcmicro_to_scimap (feature_table_path)
    assert adata.shape == (11201, 9)
    
    
@pytest.fixture
def adata():
    image_path = os.getcwd() + '/scimap/tests/scimapExampleData/scimapExampleData.h5ad'
    adata = ad.read_h5ad(image_path)
    return adata


# log1p
def test_log1p (adata):
    from scimap.preprocessing.log1p import log1p
    adata = log1p (adata, layer='log_test')
    # load expected data
    loaded_data = np.load( os.getcwd() + '/scimap/tests/expected_test_values/test_log1p.npz')['data']
    assert np.allclose(loaded_data, adata.layers['log_test']), "The arrays do not match."


# rescale
def test_rescale (adata):
    from scimap.preprocessing.rescale import rescale
    manual_gate = pd.read_csv(os.getcwd() + '/scimap/tests/scimapExampleData/manual_gates.csv')
    adata = rescale (adata, gate=manual_gate)
    # load expected data
    loaded_data = np.load( os.getcwd() + '/scimap/tests/expected_test_values/test_rescale.npz')['data']
    assert np.allclose(loaded_data, adata.X), "The arrays do not match."
    

# combat
def test_combat (adata):
    from scimap.preprocessing.combat import combat   
    adata = combat (adata, batch='ROI')
    assert adata.layers['combat'].shape == (11201, 9)
    adata = combat (adata, batch='ROI', layer='raw', label='combat_raw')
    assert adata.layers['combat_raw'].shape == (11201, 9)
    adata = combat (adata, batch='ROI', log=True, label='combat_log')
    assert adata.layers['combat_log'].shape == (11201, 9)
    adata = combat (adata, batch='ROI', layer='log', label='combat_log_layer')
    assert adata.layers['combat_log_layer'].shape == (11201, 9)
    adata = combat (adata, batch='ROI', replaceOriginal=True)
    # load expected data
    loaded_data = np.load( os.getcwd() + '/scimap/tests/expected_test_values/test_combat.npz')['data']
    assert np.allclose(loaded_data, adata.X), "The arrays do not match."
    
    

    
