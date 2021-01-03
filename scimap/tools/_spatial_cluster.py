#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 17:45:14 2020
@author: Ajit Johnson Nirmal 
Function to cluster the spatial motif data genereated by either `sm.tl.spatial_expression` 
or `sm.tl.spatial_count` methods.
"""

# Import library
from scimap.tools._cluster import cluster
import anndata as ad


# Function
def spatial_cluster (adata, df_name='spatial_count', method = 'kmeans',k=10, phenotype='phenotype',
                     n_pcs=None, resolution=1, phenograph_clustering_metric='euclidean', 
                     nearest_neighbors=30, random_state=0,label=None):
    """
    

    Parameters
    ----------
    adata : AnnData Object
    
    df_name : string, required
        Label of the spatial analysis performed.
        By default if `sm.tl.spatial_count` was run the results will be saved under `spatial_count` and
        if `sm.tl.spatial_expression` was run, the results will be saved under `spatial_expression`.
        The default is 'spatial_count'.
    method : string, optional
        Clustering method to be used- Implemented methods- kmeans, phenograph and leiden. The default is 'kmeans'.
    k : int, optional
        Number of clusters to return when using K-Means clustering. The default is 10.
    phenotype : string, optional
        The column name that contains the cluster/phenotype information. The default is 'phenotype'.
    n_pcs : int, optional
        Number of PC's to be used in leiden clustering. By default it uses all PC's.
    resolution : float, optional
        A parameter value controlling the coarseness of the clustering. 
        Higher values lead to more clusters. The default is 1.
    phenograph_clustering_metric : string, optional
        Distance metric to define nearest neighbors. Note that performance will be slower for correlation and cosine. 
        Available methods- cityblock’, ‘cosine’, ‘euclidean’, ‘manhattan’, braycurtis’, ‘canberra’, ‘chebyshev’, 
        ‘correlation’, ‘dice’, ‘hamming’, ‘jaccard’, ‘kulsinski’, ‘mahalanobis’, ‘minkowski’, ‘rogerstanimoto’, 
        ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘yule’
        The default is 'euclidean'.
    nearest_neighbors : int, optional
        Number of nearest neighbors to use in first step of graph construction. 
        This parameter is used both in leiden and phenograph clustering.
        The default is 30.
    random_state : int, optional
        Change the initialization of the optimization. The default is 0.
    label : string, optional
        Key or optional column name for the returned data, stored in `adata.obs`. The default is adata.obs[spatial_method used].

    Returns
    -------
    adata : AnnData Object
        Returns an updated anndata object with a new column. check- adata.obs[spatial_method used]
        
    Example
    -------
    adata = sm.tl.spatial_cluster (adata, k= 10, method = 'kmeans') # results will be saved under adata.obs['spatial_kmeans']

    """
    
    # Make a copy of adata to modify
    adata_copy = adata.copy()
    
    # Error check
    try:
        adata_copy.uns[df_name]
    except KeyError:
        print (str('Supplied df_name not found, please run either `sm.tl.spatial_expression` or `sm.tl.spatial_count`'))
    
    # Crete a new anndata object with the user defined spatial information
    adata_new = ad.AnnData(adata_copy.uns[df_name].fillna(0),obs=adata_copy.obs)
    
    # Create a meaningful label name
    if label is None:
        label = 'spatial_' + str(method)
    
    # Run the clustering algorithm
    adata_new = cluster (adata = adata_new,
                         method = method,
                         k=k, 
                         n_pcs=n_pcs, 
                         resolution=resolution,
                         phenograph_clustering_metric=phenograph_clustering_metric,
                         nearest_neighbors=nearest_neighbors, 
                         sub_cluster_column=phenotype,
                         use_raw=False, 
                         random_state=random_state,
                         label=label)
    
    # Get the clusters and append that to original adata object
    result = adata_new.obs[label]
    result = result.reindex(adata.obs.index)
    adata.obs[label] = result
    
    
    # Return anndata object
    return adata
    
    
    
    