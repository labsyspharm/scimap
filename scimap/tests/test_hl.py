#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 20:14:55 2024
@author: aj
Test helper function
"""

import pytest
import os
import anndata as ad
import pandas as pd
import pickle


load_pickle = lambda filename: pickle.load(open(filename, 'rb'))
#save_pickle(list(adata.obs['phenotype_renamed']), '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/softwares/scimap/scimap/tests/expected_test_values/test_rename.pkl')

# load data
@pytest.fixture
def adata():
    image_path = os.getcwd() + '/scimap/tests/scimapExampleData/scimapExampleData.h5ad'
    adata = ad.read_h5ad(image_path)
    return adata


# classify
def test_classify (adata):
    from scimap.helpers.classify import classify
    adata = classify(adata, pos=['FOXP3'], neg=['ECAD'], phenotype='phenotype')
    # load expected data
    loaded_data = load_pickle(os.getcwd() + '/scimap/tests/expected_test_values/test_classify.pkl')
    assert loaded_data == list(adata.obs['classify']), "The lists do not match."

# rename
def test_rename (adata):
    from scimap.helpers.rename import rename
    name= {'renamed': ['Dendritic cells', 'NK cells']}
    adata = rename (adata, name, from_column='phenotype', to_column='phenotype_renamed')
    # load expected data
    loaded_data = load_pickle(os.getcwd() + '/scimap/tests/expected_test_values/test_rename.pkl')
    assert loaded_data == list(adata.obs['phenotype_renamed']), "The lists do not match."
    
# dropFeatures
def test_dropFeatures (adata):
    from scimap.helpers.dropFeatures import dropFeatures
    adata = dropFeatures(adata, drop_markers=['ELANE', 'NCAM'])
    assert len(adata.var.index) == 7

# merge_adata_obs
def test_merge_adata_obs (adata):
    from scimap.helpers.merge_adata_obs import merge_adata_obs
    bdata = adata.copy()
    bdata.obs['new_col'] = bdata.obs['imageid']
    combined_adata = merge_adata_obs(adata=[adata, bdata])
    assert len(combined_adata.obs.columns) == 14

# scimap_to_csv
def test_scimap_to_csv (adata):
    from scimap.helpers.scimap_to_csv import scimap_to_csv
    data = scimap_to_csv(adata)
    assert data.shape == (11201, 22)
    
# downloadDemoData
def downloadDemoData ():
    from scimap.helpers.downloadDemoData import downloadDemoData
    download_directory = os.getcwd() + '/demodata'
    downloadDemoData (download_directory, api_url="https://zenodo.org/api/records/11054442")
    downloaded_data = pd.read_csv(os.getcwd() + '/demodata/manual_gates.csv')
    assert downloaded_data.shape == (9, 2)
        
# addROI_omero

# animate