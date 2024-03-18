#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Created on Sun Mar 17 16:09:22 2024
#@author: Ajit Johnson Nirmal
#Create a neighbourhood graph

"""
!!! abstract "Short Description"
    `sm.pp.nGraph` constructs a k-neighbors graph from single-cell data contained within an AnnData object. It offers options for data preprocessing such as standard scaling and principal component analysis (PCA) before graph construction. The resulting graph is stored in the `.obsp['connectivities']` of the AnnData object. The function accommodates data from the raw layer, a specified layer, or the default data layer (`adata.X`), and allows for specifying the number of neighbors to define connectivity. 

## Function
"""


# lib
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import numpy as np
import igraph as ig


# function

def nGraph (adata, 
            layer='raw', 
            standardScale=False, 
            runPCA=False, 
            k_neighbors=15):
    
    """
Generates a k-neighbors graph from high-dimensional single-cell data, with options for preprocessing steps 
such as standard scaling and principal component analysis (PCA).

Parameters:
    adata (AnnData): 
        An AnnData object containing single-cell data. Must have `.X` for data matrix, `.raw.X` for raw data, 
        and `.layers` for additional data layers.
        
    layer (str, optional): 
        Specifies which layer of the `adata` to use for graph construction. The default is 'raw', indicating 
        that `adata.raw.X` will be used. If `None`, `adata.X` will be utilized. Otherwise, specifies the key 
        to use data from `adata.layers`.
        
    standardScale (bool, optional): 
        If `True`, applies standard scaling to the data, making the mean of each feature 0 and the variance 1. 
        
    runPCA (bool, optional): 
        If `True`, performs principal component analysis on the data and uses the principal components for 
        graph construction. This is often done to reduce dimensionality and noise. 
        
    k_neighbors (int, optional): 
        The number of neighbors to use for k-neighbors graph construction. This parameter determines the 
        connectivity of the graph. Defaults to 15.

Returns:
    adata (annData): 
        The input `adata` object is returned after adding the k-neighbors graph to `.obsp['connectivities']`.

Examples:
    ```python
    
    # Example 1: Basic usage with raw layer data and default settings
    adata = sm.pp.nGraph(adata)

    # Example 2: Using data from default layer, with standard scaling and PCA applied, specifying k_neighbors
    adata = sm.pp.nGraph(adata, layer=None, standardScale=True, runPCA=True, k_neighbors=20)
    ```
"""
    
    # prepare data
    if layer == 'raw':
        data = adata.raw.X.copy()
    elif layer is None:
        data = adata.X.copy()
    else:
        data = adata.layers[layer].copy()
    
    if standardScale:
        scaler = StandardScaler()
        data = scaler.fit_transform(data)
    
    if runPCA:
        # Initialize PCA object
        pca = PCA(n_components=None)  # 'None' to obtain all PCs
        # Fit PCA on the data
        pca.fit(data.T)
        # Transform the data
        X_pca = pca.transform(data.T)
        # X_pca now contains the principal components
        data = pca.components_.T
        
    
    # Generate a k-neighbors graph from the data
    graph = kneighbors_graph(X=data, n_neighbors=k_neighbors, mode='connectivity')
    adata.obsp['connectivities'] = graph
    
    # return graph
    return adata
