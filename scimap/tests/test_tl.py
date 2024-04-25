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
    assert adata.obs['phenotype'].iloc[0] == 'Immune'

#cluster
def test_cluster (adata):
    from scimap.tools.cluster import cluster
    adata = cluster (adata,  method = 'kmeans', k= 5, use_raw = True)
    assert adata.obs['kmeans'].value_counts().iloc[4] == 81
    

#umap
def test_umap (adata):
    from scimap.tools.umap import umap
    adata = umap(adata, label='umap_test')
    assert adata.obsm['umap_test'].shape == (11201, 2)


#foldchange
def test_foldchange (adata):
    from scimap.tools.foldchange import foldchange
    adata = foldchange(adata, from_group = 'ROI1', imageid='ROI')
    assert round(adata.uns['foldchange_fc']['Immune'], 2).values[0] == 1.87


#spatial_distance
def test_spatial_distance (adata):
    from scimap.tools.spatial_distance import spatial_distance
    adata = spatial_distance (adata)
    assert round(adata.uns['spatial_distance']['Immune'].iloc[0], 2) == 521.27

    
#spatial_interaction
def test_spatial_interaction (adata):
    from scimap.tools.spatial_interaction import spatial_interaction
    adata = spatial_interaction (adata,  method = 'knn', knn= 5, permutation = 10)
    assert adata.uns['spatial_interaction'] is not None
    

#spatial_count &
#spatial_cluster
def test_spatial_count (adata):
    from scimap.tools.spatial_count import spatial_count
    adata = spatial_count (adata, phenotype='phenotype',method='knn',radius=5)
    assert adata.uns['spatial_count'] is not None
    
    # test spatial cluster
    from scimap.tools.spatial_cluster import spatial_cluster
    adata = spatial_cluster(adata, df_name='spatial_count')
    assert adata.obs['spatial_kmeans'].value_counts().iloc[3] == 1213


#spatial_lda
def test_spatial_lda (adata):
    from scimap.tools.spatial_lda import spatial_lda
    adata = spatial_lda (adata, num_motifs=10, radius=30)
    assert round(adata.uns['spatial_lda']['Motif_0'].iloc[0], 3) == 0.05


#spatial_expression
def test_spatial_expression (adata):
    from scimap.tools.spatial_expression import spatial_expression
    adata = spatial_expression (adata, method='radius', radius=30, use_raw=True, label='spatial_expression')
    assert round(adata.uns['spatial_expression']['ECAD'].sum(),0) == 81404.0


#spatial_pscore
def test_spatial_pscore (adata):
    from scimap.tools.spatial_pscore import spatial_pscore
    adata = spatial_pscore (adata,  proximity= ['Immune', 'ECAD+'])
    assert adata.uns['spatial_pscore']['Immune_ECAD+'].values[0] == 2095


#spatial_aggregate
def test_spatial_aggregate (adata):
    from scimap.tools.spatial_aggregate import spatial_aggregate
    adata = spatial_aggregate (adata, purity = 60, phenotype='phenotype', method='knn', radius=10)
    assert adata.obs['spatial_aggregate'].value_counts()['Immune'] == 3439


#spatial_similarity_search
def test_spatial_similarity_search (adata):
    from scimap.tools.spatial_similarity_search import spatial_similarity_search
    adata = spatial_similarity_search (adata, ROI_column='ROI', similarity_threshold=0.6,  method='knn', radius=10)
    assert adata.obs['spatial_similarity_search_ROI2'].value_counts()['similar_to_ROI'] == 4786

