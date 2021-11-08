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
    
# Testing the spatial_expression function
def test_spatial_expression(adata):
    from scimap.tools._spatial_expression import spatial_expression
    
    adata = spatial_expression (adata, x_coordinate='X_position',y_coordinate='Y_position',
                                method='radius', radius=30, imageid='ImageId', 
                                use_raw=True,subset=None,label='spatial_expression')
    a = round(adata.uns['spatial_expression']['CD3D'].sum(),0)
    
    # test
    assert a == 19920
    
# Testing the spatial_aggregate function
def test_spatial_aggregate(adata):
    from scimap.tools._spatial_aggregate import spatial_aggregate
    
    adata = spatial_aggregate (adata, x_coordinate='X_position',y_coordinate='Y_position',
                           purity = 60, phenotype='phenotype', method='radius', radius=30,
                           imageid='ImageId',subset=None,label='spatial_aggregate')
    a = adata.obs['spatial_aggregate'].value_counts()['M2 Macrophages']
    
    # test
    assert a == 20
    
# Testing cluster function
def test_cluster(adata):
    from scimap.tools._cluster import cluster
    adata = cluster (adata,  method = 'kmeans', k= 5, use_raw = True)
    a = adata.obs['kmeans'].value_counts()[4]
    
    #test
    assert a == 252

# Testing spatial_interaction function
def test_spatial_interaction(adata):
    from scimap.tools._spatial_interaction import spatial_interaction
    adata = spatial_interaction (adata,  method = 'knn', knn= 5, permutation = 10, imageid='ImageId',x_coordinate='X_position',y_coordinate='Y_position')
    a = adata.uns['spatial_interaction']
    
    #test
    assert a is not None
    
# Testing spatial_distance function
def test_spatial_distance(adata):
    from scimap.tools._spatial_distance import spatial_distance
    adata = spatial_distance (adata, imageid='ImageId',x_coordinate='X_position',y_coordinate='Y_position')
    a = adata.uns['spatial_distance']
    
    #test
    assert a is not None
    
# Testing spatial_pscore function
def test_spatial_pscore(adata):
    from scimap.tools._spatial_pscore import spatial_pscore
    adata = spatial_pscore (adata, imageid='ImageId',x_coordinate='X_position',y_coordinate='Y_position', 
                            score_by='ImageId', proximity= ['Tumor CD30+', 'M2 Macrophages'])
    a = adata.uns['spatial_pscore']['All Cells'].values
    
    # test
    assert a == 3029

# Testing spatial_lda function
def test_spatial_lda (adata):
    from scimap.tools._spatial_lda import spatial_lda
    adata = spatial_lda (adata, num_motifs=10, radius=30,imageid='ImageId',
                         x_coordinate='X_position',y_coordinate='Y_position')
    a = round(adata.uns['spatial_lda']['Motif_0'][0], 3)
    
    # test
    assert a == 0.775
    
# Testing foldchange function
def test_foldchange (adata):
    from scimap.tools._foldchange import foldchange
    import numpy as np
    # prepare
    x = np.repeat('ROI1', round(adata.shape[0]/2))
    y = np.repeat('ROI2', adata.shape[0] - round(adata.shape[0]/2))
    z = np.concatenate((x, y))    
    adata.obs['ROI'] = z
    # start test
    adata = foldchange(adata, from_group = 'ROI1', imageid='ROI')
    a = round(adata.uns['foldchange_fc']['ASMA+ cells'], 2).values
    
    # test
    assert a == 1.27





























