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
    
def test_rescale(adata):
    import pandas as pd
    import numpy as np
    from scimap.preprocessing._rescale import rescale
    manual_gate = pd.DataFrame({'marker': ['CD3D', 'KI67'], 'gate': [7, 8]})
    adata = rescale (adata, gate=manual_gate, failed_markers=['CD20', 'CD21'])
    assert np.round(adata[:,'CD3D'].X[0],2) == 0.37
    
    