#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 10:51:50 2024
@author: aj
Test for tools
"""

import pytest
import os
import anndata as ad
import pandas as pd
import numpy as np
import pickle

# def save_pickle(obj, filename): pickle.dump(obj, open(filename, 'wb'))
#save_pickle(list(adata.obs['spatial_similarity_search_ROI2']), '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/softwares/scimap/scimap/tests/expected_test_values/test_spatial_similarity_search.pkl')
load_pickle = lambda filename: pickle.load(open(filename, 'rb'))

# load data
@pytest.fixture
def adata():
    image_path = os.getcwd() + '/scimap/tests/scimapExampleData/scimapExampleData.h5ad'
    adata = ad.read_h5ad(image_path)
    return adata


#phenotype_cells
def test_phenotype (adata):
    from scimap.tools.phenotype_cells import phenotype_cells
    # Load phenotype and test phenotyping
    phenotype = pd.read_csv(os. getcwd() + '/scimap/tests/scimapExampleData/phenotype_workflow.csv')
    adata = phenotype_cells (adata, phenotype=phenotype, gate = 0.5, label="phenotype")
    # load expected data
    loaded_data = load_pickle(os.getcwd() + '/scimap/tests/expected_test_values/test_phenotype.pkl')
    assert loaded_data == list(adata.obs['phenotype']), "The lists do not match."


#cluster
def test_cluster (adata):
    from scimap.tools.cluster import cluster
    adata = cluster (adata,  method = 'kmeans', k= 5, use_raw = True)
    # load expected data
    loaded_data = load_pickle(os.getcwd() + '/scimap/tests/expected_test_values/test_cluster.pkl')
    assert loaded_data == list(adata.obs['kmeans']), "The lists do not match."
    
#umap
def test_umap (adata):
    from scimap.tools.umap import umap
    adata = umap(adata, label='umap_test')
    assert adata.obsm['umap_test'].shape == (11201, 2)


#foldchange
def test_foldchange (adata):
    from scimap.tools.foldchange import foldchange
    adata = foldchange(adata, from_group = 'ROI1', imageid='ROI')
    # load expected data
    loaded_data = np.load( os.getcwd() + '/scimap/tests/expected_test_values/test_foldchange.npz')['data']
    assert np.allclose(loaded_data, adata.uns['foldchange_fc'].to_numpy() ), "The arrays do not match."


#spatial_distance
def test_spatial_distance (adata):
    from scimap.tools.spatial_distance import spatial_distance
    adata = spatial_distance (adata)
    # load expected data
    loaded_data = np.load( os.getcwd() + '/scimap/tests/expected_test_values/test_spatial_distance.npz')['data']
    assert np.allclose(loaded_data, adata.uns['spatial_distance'].to_numpy()), "The arrays do not match."
    

#spatial_interaction
def test_spatial_interaction (adata):
    from scimap.tools.spatial_interaction import spatial_interaction
    adata = spatial_interaction (adata,  method = 'knn', knn= 5, permutation = 10)
    assert adata.uns['spatial_interaction'] is not None
    

#spatial_count
def test_spatial_count (adata):
    from scimap.tools.spatial_count import spatial_count
    adata = spatial_count (adata, phenotype='phenotype',method='knn',radius=5)
    # load expected data
    loaded_data = np.load( os.getcwd() + '/scimap/tests/expected_test_values/test_spatial_count.npz')['data']
    assert np.allclose(loaded_data, adata.uns['spatial_count'].to_numpy()), "The arrays do not match."
    
    
#spatial_cluster
def test_spatial_cluster (adata):
    from scimap.tools.spatial_cluster import spatial_cluster
    adata = spatial_cluster (adata, df_name='spatial_count_test')
    # load expected data
    #loaded_data = load_pickle(os.getcwd() + '/scimap/tests/expected_test_values/test_spatial_cluster.pkl')
    #assert loaded_data == list(adata.obs['spatial_kmeans']), "The lists do not match."
    

#spatial_lda
def test_spatial_lda (adata):
    from scimap.tools.spatial_lda import spatial_lda
    adata = spatial_lda (adata, num_motifs=10, radius=30)
    assert adata.uns['spatial_lda'] is not None


#spatial_expression
def test_spatial_expression (adata):
    from scimap.tools.spatial_expression import spatial_expression
    adata = spatial_expression (adata, method='radius', radius=30, use_raw=True, label='spatial_expression')
    # load expected data
    loaded_data = np.load( os.getcwd() + '/scimap/tests/expected_test_values/test_spatial_expression.npz')['data']
    assert np.allclose(loaded_data, adata.uns['spatial_expression'].to_numpy()), "The arrays do not match."


#spatial_pscore
def test_spatial_pscore (adata):
    from scimap.tools.spatial_pscore import spatial_pscore
    adata = spatial_pscore (adata,  proximity= ['Immune', 'ECAD+'])
    # load expected data
    loaded_data = np.load( os.getcwd() + '/scimap/tests/expected_test_values/test_spatial_pscore.npz')['data']
    assert np.allclose(loaded_data, adata.uns['spatial_pscore'].to_numpy()), "The arrays do not match."


#spatial_aggregate
def test_spatial_aggregate (adata):
    from scimap.tools.spatial_aggregate import spatial_aggregate
    adata = spatial_aggregate (adata, purity = 60, phenotype='phenotype', method='knn', radius=10)
    # load expected data
    loaded_data = load_pickle(os.getcwd() + '/scimap/tests/expected_test_values/test_spatial_aggregate.pkl')
    assert loaded_data == list(adata.obs['spatial_aggregate']), "The lists do not match."
    

#spatial_similarity_search
def test_spatial_similarity_search (adata):
    from scimap.tools.spatial_similarity_search import spatial_similarity_search
    adata = spatial_similarity_search (adata, ROI_column='ROI', similarity_threshold=0.6,  method='knn', radius=10)
    # load expected data
    loaded_data = load_pickle(os.getcwd() + '/scimap/tests/expected_test_values/test_spatial_similarity_search.pkl')
    assert loaded_data == list(adata.obs['spatial_similarity_search_ROI2']), "The lists do not match."
    











