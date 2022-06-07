#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Sat May 28 19:13:24 2022
# @author: Ajit Johnson Nirmal
# UMAP Function

"""
!!! abstract "Short Description"
    `sm.tl.umap`: The function allows users to perform dimensionality reduction using UMAP 
## Function
"""

# libs
import umap as um
import numpy as np


# function
def umap (adata, use_layer=None, use_raw=False, log=False,
          n_neighbors=15, n_components=2, metric='euclidean',min_dist=0.1, random_state=0, 
          label='umap', **kwargs):
    """
Parameters:

    adata : AnnData Object  
    
    use_layer : string, optional  
        Pass name of any `Layer` in AnnData. The default is `None` and `adata.X` is used.
        
    use_raw : bool, optional  
        If set to `True`, values in `adata.raw.X` will be used. The default is False.
        
    log : bool, optional  
        If set to `True`, the data will natural log transformed using `np.log1p()`. The default is False.
        
    n_neighbors : int, optional  
        The size of local neighborhood (in terms of number of neighboring sample points) 
        used for manifold approximation. Larger values result in more global views of the manifold, 
        while smaller values result in more local data being preserved. In general 
        values should be in the range 2 to 100. The default is 15.
        
    n_components : int, optional  
        The dimension of the space to embed into. This defaults to 2 to provide easy visualization, 
        but can reasonably be set to any integer value in the range 2 to 100. The default is 2.
        
    metric : TYPE, optional  
        The metric to use to compute distances in high dimensional space. 
        Check `https://umap-learn.readthedocs.io/en/latest/api.html` for all available 
        options. The default is 'euclidean'.
        
    min_dist : float, optional  
        The effective minimum distance between embedded points. Smaller values will 
        result in a more clustered/clumped embedding where nearby points on the manifold 
        are drawn closer together, while larger values will result on a more even 
        dispersal of points. The value should be set relative to the spread value, 
        which determines the scale at which embedded points will be spread out. The default is 0.1.
        
    random_state : int, optional  
        If int, random_state is the seed used by the random number generator; 
        If RandomState instance, random_state is the random number generator; 
        If None, the random number generator is the RandomState instance used by np.random. The default is 0.
        
    label : string, optional  
        Key used to save the embeddings. 
        Check `adata.obsm[label]` for results. The default is 'umap'.
    
    **kwargs : Other `umap` parameters. Check `https://umap-learn.readthedocs.io/en/latest/api.html`  

Returns:

    adata : Modified AnnData object
        Embedding stored in `adata.obsm[label]`.
    
Example:
```python
# Run UMAP
adata = sm.tl.umap(adata)
    
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
