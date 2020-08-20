#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 09:32:49 2020
@author: Ajit Johnson Nirmal
Tests
"""

import pytest
import sys, os
import anndata as ad

#os.chdir ("/Users/aj/Dropbox (Partners HealthCare)/packages/scimap")

@pytest.fixture
def adata():
    image_path = os. getcwd() + '/scimap/tests/_data/example_data.h5ad'
    adata = ad.read(image_path)
    return adata

# Testing the phenotyping function
def test_phenotype(adata):
    import pandas as pd
    import numpy as np
    from scimap.tools._phenotype_cells import phenotype_cells

    # Load phenotype and test phenotyping
    phenotype = pd.read_csv(os. getcwd() + '/scimap/tests/_data/phenotype_workflow.csv')
    adata = phenotype_cells (adata, phenotype=phenotype, gate = 0.5, label="phenotype")
    a = adata.obs['phenotype'][0]

    # test
    assert a == 'M2 Macrophages'

# Testing the spatial_count function
def test_spatial_count(adata):
    from scimap.tools._spatial_count import spatial_count
    
    adata = spatial_count (adata,x_coordinate='X_position',y_coordinate='Y_position',
                           phenotype='phenotype',method='radius',radius=30,
                           imageid='ImageId',subset=None,label='spatial_count_radius')
    a = round(adata.uns['spatial_count_radius']['ASMA+ cells'].sum(),0)
    
    # test
    assert a == 330
