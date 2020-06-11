#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 09:32:49 2020
@author: Ajit Johnson Nirmal
Tests
"""

import pytest
import sys, os

@pytest.fixture
def adata():
    from scimap.preprocessing._mcmicro_to_scimap import mcmicro_to_scimap
    image_path = [os. getcwd() + '/scimap/tests/_data/example_data.csv']
    adata = mcmicro_to_scimap (image_path, drop_markers=['BG1', 'BG2', 'BG3'])
    return adata

    
def test_mcmicro_to_scimap(adata):
    assert adata.shape == (3029, 33)
    
def test_rescale_phenotype(adata):
    import pandas as pd
    import numpy as np
    from scimap.preprocessing._rescale import rescale
    from scimap.tools._phenotype_cells import phenotype_cells
    
    # test rescaling data
    manual_gate = pd.DataFrame({'marker': ['CD3D', 'KI67'], 'gate': [7, 8]})
    adata = rescale (adata, gate=manual_gate, failed_markers=['CD20', 'CD21'])
    #assert np.round(adata[:,'CD3D'].X[0],2) == 0.37
    a = np.round(adata[:,'CD3D'].X[0],2)
    
    # Load phenotype and test phenotyping
    phenotype = pd.read_csv(os. getcwd() + '/scimap/tests/_data/phenotype_workflow.csv')
    adata = phenotype_cells (adata, phenotype=phenotype, gate = 0.5, label="phenotype") 
    #assert adata.obs['phenotype'][0] == 'M2 Macrophages'
    b = adata.obs['phenotype'][0]
    
    # test
    assert (a, b) == (0.37, 'M2 Macrophages')
    

    