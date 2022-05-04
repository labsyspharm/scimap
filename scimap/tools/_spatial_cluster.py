#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# # Created on Mon Dec 14 17:45:14 2020
# @author: Ajit Johnson Nirmal 
"""
!!! abstract "Short Description"
    `sm.tl.spatial_cluster`: This function allows users to cluster the spatial neighbourhood matrix 
    genereated by either `sm.tl.spatial_expression` or `sm.tl.spatial_count` methods.

## Function
"""

# Import library
from scimap.tools._cluster import cluster
import anndata as ad
import argparse
import sys
import pathlib

def main(argv=sys.argv):
    parser = argparse.ArgumentParser(
        description='This function allows users to cluster the spatial neighbourhood matrix genereated by either `sm.tl.spatial_expression` or `sm.tl.spatial_count` methods.'
    )
    parser.add_argument(
        '--adata', required=True, 
        help='AnnData object loaded into memory or path to AnnData object.'
    )
    parser.add_argument(
        '--df_name', type=str, required=False, default='spatial_count',
        help='Label of the spatial analysis performed. By default if `sm.tl.spatial_count` was run the results will be saved under `spatial_count` and if `sm.tl.spatial_expression` was run, the results will be saved under `spatial_expression`.'
    )
    parser.add_argument(
        '--method', type=str, required=False, default='kmeans',
        help='Clustering method to be used- Implemented methods- kmeans, phenograph, leiden and parc.'
    )
    parser.add_argument(
        '--k', type=int, required=False, default=10,
        help='Number of clusters to return when using K-Means clustering.'
    )
    parser.add_argument(
        '--n_pcs', type=int, required=False, default=None,
        help='Number of PCs to be used in leiden clustering. By default it uses all PCs.'
    )
    parser.add_argument(
        '--resolution', type=float, required=False, default=1,
        help='A parameter value controlling the coarseness of the clustering. Higher values lead to more clusters.'
    )
    parser.add_argument(
        '--phenograph_clustering_metric', type=str, required=False, default='euclidean',
        help='Distance metric to define nearest neighbors. Note that performance will be slower for correlation and cosine. Available methods- cityblock’, ‘cosine’, ‘euclidean’, ‘manhattan’, braycurtis’, ‘canberra’, ‘chebyshev’,  ‘correlation’, ‘dice’, ‘hamming’, ‘jaccard’, ‘kulsinski’, ‘mahalanobis’, ‘minkowski’, ‘rogerstanimoto’, ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘yule’'
    )
    parser.add_argument(
        '--nearest_neighbors', type=int, required=False, default=30,
        help='Number of nearest neighbors to use in first step of graph construction. This parameter is used both in leiden and phenograph clustering.'
    )
    parser.add_argument(
        '--random_state', type=int, required=False, default=0,
        help='Change the initialization of the optimization.'
    )
    parser.add_argument(
        '--label', type=str, required=False, default=None,
        help='Key or optional column name for the returned data, stored in `adata.obs`. The default is adata.obs [method used].'
    )
    parser.add_argument(
        '--output_dir', type=str, required=False, default=None,
        help='Path to output directory.'
    )
    args = parser.parse_args(argv[1:])
    print(vars(args))
    spatial_cluster(**vars(args))


# Function
def spatial_cluster (adata, df_name='spatial_count', method = 'kmeans',k=10,
                     n_pcs=None, resolution=1, phenograph_clustering_metric='euclidean', 
                     nearest_neighbors=30, random_state=0,label=None, output_dir=None):
    """
    

Parameters:
    adata : AnnData object loaded into memory or path to AnnData object.
    
    df_name : string, required  
        Label of the spatial analysis performed.
        By default if `sm.tl.spatial_count` was run the results will be saved under `spatial_count` and
        if `sm.tl.spatial_expression` was run, the results will be saved under `spatial_expression`.

    method : string, optional  
        Clustering method to be used- Implemented methods- kmeans, phenograph and leiden.

    k : int, optional  
        Number of clusters to return when using K-Means clustering.

    phenotype : string, optional  
        The column name that contains the cluster/phenotype information.

    n_pcs : int, optional  
        Number of PC's to be used in leiden clustering. By default it uses all PC's.

    resolution : float, optional  
        A parameter value controlling the coarseness of the clustering. 
        Higher values lead to more clusters.

    phenograph_clustering_metric : string, optional  
        Distance metric to define nearest neighbors. Note that performance will be slower for correlation and cosine. 
        Available methods- cityblock’, ‘cosine’, ‘euclidean’, ‘manhattan’, braycurtis’, ‘canberra’, ‘chebyshev’, 
        ‘correlation’, ‘dice’, ‘hamming’, ‘jaccard’, ‘kulsinski’, ‘mahalanobis’, ‘minkowski’, ‘rogerstanimoto’, 
        ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘yule’

    nearest_neighbors : int, optional  
        Number of nearest neighbors to use in first step of graph construction. 
        This parameter is used both in leiden and phenograph clustering.

    random_state : int, optional  
        Change the initialization of the optimization.

    label : string, optional  
        Key or optional column name for the returned data, stored in `adata.obs`. The default is adata.obs [spatial_method used].
    
    output_dir : string, optional  
        Path to output directory.

Returns:
    adata : AnnData Object  
        Returns an updated anndata object with a new column. check- adata.obs [spatial_method used]
        
Example:
```python
    adata = sm.tl.spatial_cluster (adata, k= 10, method = 'kmeans') # results will be saved under adata.obs['spatial_kmeans']
```
    """

    # Load the andata object    
    if isinstance(adata, str):
        imid = str(adata.rsplit('/', 1)[-1])
        adata = ad.read(adata)
    else:
        adata = adata
    
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
                         use_raw=False, 
                         random_state=random_state,
                         label=label)
    
    # Get the clusters and append that to original adata object
    result = adata_new.obs[label]
    result = result.reindex(adata.obs.index)
    adata.obs[label] = result
    
    
    # Save data if requested
    if output_dir is not None:
        output_dir = pathlib.Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        adata.write(output_dir / imid)
    else:    
        # Return data
        return adata

if __name__ == '__main__':
    main()