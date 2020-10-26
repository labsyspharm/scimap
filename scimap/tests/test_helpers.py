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
    import anndata as ad    
    adata = ad.read(os. getcwd() + '/scimap/tests/_data/example_data.h5ad')
    return adata

# test rename function 
def test_rename (adata):
    from scimap.helpers._rename import rename
    name= {'T cells': ['CD8 T cells', 'CD4 T cells']}
    adata = rename (adata, name, from_column='phenotype', to_column='phenotype_renamed')
    assert adata.obs['phenotype_renamed'].value_counts()['T cells'] == 97
    

# test classify function
def test_classify (adata):
    from scimap.helpers._classify import classify
    adata = classify(adata, pos=['CD3D','FOXP3'], neg=['ASMA'], 
                     subclassify_phenotype=['T cells','Regulatory T cells'])
    
    a = adata.obs['classify'].value_counts()['passed_classify']
    assert a == 16