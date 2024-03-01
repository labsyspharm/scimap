#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Sat May 28 19:13:24 2022
# @author: Ajit Johnson Nirmal
# UMAP Function

"""
!!! abstract "Short Description"
    `sm.tl.umap`: This function enables dimensionality reduction on high-dimensional 
    datasets using UMAP, allowing for the visualization of complex data structures 
    in a lower-dimensional space. It supports customization through various parameters, 
    including data source selection, logarithmic transformation, and manifold 
    approximation settings, accommodating a wide range of analytical needs. Results 
    are stored in `adata.obsm`, ready for subsequent visualization or analysis.
    
## Function
"""

# libs
import umap as um
import numpy as np


# function
def umap (adata, 
          use_layer=None, 
          use_raw=False, 
          log=False,
          n_neighbors=15, 
          n_components=2, 
          metric='euclidean',
          min_dist=0.1, 
          random_state=0, 
          label='umap', **kwargs):
    """
Parameters:
        adata (anndata.AnnData):  
            Annotated data matrix or path to an AnnData object, containing spatial gene expression data.
            
        use_layer (str, optional):  
            Specifies a layer in `adata.layers` for UMAP. Defaults to using `adata.X`.
            
        use_raw (bool, optional):  
            Whether to use `adata.raw.X` for the analysis.
            
        log (bool, optional):  
            Applies natural log transformation to the data if `True`.
            
        n_neighbors (int, optional):  
            Number of neighboring points used in manifold approximation.
            
        n_components (int, optional):  
            Dimensionality of the target embedding space.
            
        metric (str, optional):  
            Metric used to compute distances in high-dimensional space.
            
        min_dist (float, optional):  
            Effective minimum distance between embedded points.
            
        random_state (int, optional):  
            Seed used by the random number generator for reproducibility.
            
        label (str, optional):  
            Key for storing UMAP results in `adata.obsm`.

Returns:
        adata (anndata.AnnData):  
            The input `adata` object, updated with UMAP embedding results in `adata.obsm[label]`.

Example:
        ```python
        
        # Basic UMAP reduction
        adata = sm.tl.umap(adata, n_neighbors=15, min_dist=0.1, label='umap_basic')
    
        # UMAP using specific layer and log transformation
        adata = sm.tl.umap(adata, use_layer='counts', use_raw=True, log=True, n_neighbors=30, min_dist=0.05, label='umap_layer_log')
    
        # UMAP with a different metric and higher dimensionality
        adata = sm.tl.umap(adata, metric='manhattan', n_components=3, n_neighbors=50, label='umap_manhattan_3d')
        
        # plot results
        sm.pl.umap(adata)
        
        ```
    """

    # adata_layer=None;use_raw=False;log=False;n_neighbors=15;n_components=2;metric='euclidean';min_dist=0.1;
    # random_state=0;
    # load data
    if use_layer is not None:
        data = adata.layers[use_layer]
    elif use_raw is True:
        data = adata.raw.X
    else:
        data = adata.X
    
    # log the data if user requests
    if log is True:
        data = np.log1p(data)
    
        
    # embedding
    embedding = um.UMAP(n_neighbors=n_neighbors,
                          n_components=n_components,
                          metric=metric,
                          min_dist=min_dist,
                          random_state=random_state).fit_transform(data)
    
    # plot
    #plt.scatter(embedding[:, 0], embedding[:, 1], s=5)
    
    # return data
    adata.obsm[label] = embedding
    return adata
