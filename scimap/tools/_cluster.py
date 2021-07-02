#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Mon Oct 12 17:03:56 2020
# @author: Ajit Johnson Nirmal
"""
!!! abstract "Short Description"
    `sm.tl.cluster`: This function allows users to cluster the dataset. 
    The function supports four clustering algorithm (kmeans, phenograph, leiden and parc).

## Function
"""

# Import library
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from sklearn.cluster import KMeans
import argparse
import sys
import anndata
import pathlib
try:
    import parc
except:
    pass


def main(argv=sys.argv):
    parser = argparse.ArgumentParser(
        description='This function allows users to cluster the dataset. The function supports four clustering algorithm (kmeans, phenograph, leiden and parc).'
    )
    parser.add_argument(
        '--adata', required=True, 
        help='AnnData object loaded into memory or path to AnnData object.'
    )
    parser.add_argument(
        '--method', type=str, required=False, default='kmeans',
        help='Clustering method to be used- Implemented methods- kmeans, phenograph, leiden and parc.'
    )
    parser.add_argument(
        '--subset_genes', type=list, required=False, default=None,
        help='Pass a list of genes [`CD3D`, `CD20`, `KI67`] that should be included for the purpose of clustering. By default the algorithm uses all genes in the dataset.'
    )
    parser.add_argument(
        '--sub_cluster', type=bool, required=False, default=False,
        help='If the user has already performed clustering or phenotyping previously and would like to sub-cluster within a particular cluster/phenotype, this option can be used.'
    )
    parser.add_argument(
        '--sub_cluster_column', type=str, required=False, default='phenotype',
        help='The column name that contains the cluster/phenotype information to be sub-clustered. This is only required when sub_cluster is set to True.'
    )
    parser.add_argument(
        '--sub_cluster_group', type=list, required=False, default=None,
        help='By default the program will sub-cluster all groups within column passed through the argument sub_cluster_column. If user wants to sub cluster only a subset of phenotypes/clusters this option can be used. Pass them as list e.g. ["tumor", "b cells"].'
    )
    parser.add_argument(
        '--parc_small_pop', type=int, required=False, default=50,
        help='Smallest cluster population to be considered a community in PARC clustering.'
    )
    parser.add_argument(
        '--parc_too_big_factor', type=float, required=False, default=0.4,
        help='If a cluster exceeds this share of the entire cell population, then the PARC will be run on the large cluster. at 0.4 it does not come into play.'
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
        '--use_raw', type=bool, required=False, default=True,
        help='If True, raw data will be used for clustering. If False, normalized/scaled data within `adata.X` will be used.'
    )
    parser.add_argument(
        '--log', type=bool, required=False, default=True,
        help='If `True`, the log of raw data is used. Set use_raw = `True` for this to take effect. '
    )
    parser.add_argument(
        '--random_state', type=int, required=False, default=0,
        help='Change the initialization of the optimization.'
    )
    parser.add_argument(
        '--collapse_labels', type=bool, required=False, default=False,
        help='While sub clustering only a few phenotypes/clusters, this argument helps to group all the other phenotypes/clusters into a single category- Helps in visualisation.'
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
    cluster(**vars(args))



def cluster (adata, method='kmeans', subset_genes=None,
             sub_cluster=False, sub_cluster_column='phenotype', sub_cluster_group = None,
             parc_small_pop= 50, parc_too_big_factor=0.4, 
             k= 10, n_pcs=None, resolution=1, 
             phenograph_clustering_metric='euclidean', nearest_neighbors= 30, 
             use_raw=True, log=True, random_state=0, collapse_labels= False,
             label=None, output_dir=None):
    """
    
Parameters:

    adata : AnnData Object

    method : string, optional  
        Clustering method to be used- Implemented methods- kmeans, phenograph, leiden and parc.

    subset_genes : list, optional  
        Pass a list of genes ['CD3D', 'CD20', 'KI67'] that should be included for the purpose of clustering. 
        By default the algorithm uses all genes in the dataset.

    sub_cluster : Boolian, optional  
        If the user has already performed clustering or phenotyping previously and would like to
        sub-cluster within a particular cluster/phenotype, this option can be used.

    sub_cluster_column : string, optional  
        The column name that contains the cluster/phenotype information to be sub-clustered. 
        This is only required when sub_cluster is set to True.

    sub_cluster_group : list, optional  
        By default the program will sub-cluster all groups within column passed through the argument sub_cluster_column.
        If user wants to sub cluster only a subset of phenotypes/clusters this option can be used.
        Pass them as list e.g. ["tumor", "b cells"].     

    parc_small_pop : int, optional  
        Smallest cluster population to be considered a community in PARC clustering.

    parc_too_big_factor : float, optional  
        If a cluster exceeds this share of the entire cell population, then the PARC will be run on 
        the large cluster. at 0.4 it does not come into play.

    k : int, optional  
        Number of clusters to return when using K-Means clustering.

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

    use_raw : bool, optional  
        If True, raw data will be used for clustering. 
        If False, normalized/scaled data within `adata.X` will be used.

    log : boolian, optional  
        If `True`, the log of raw data is used. Set use_raw = `True` for this to take effect. 

    random_state : int, optional  
        Change the initialization of the optimization.

    collapse_labels : bool, optional  
        While sub clustering only a few phenotypes/clusters, this argument helps to 
        group all the other phenotypes/clusters into a single category- 
        Helps in visualisation.

    label : string, optional  
        Key or optional column name for the returned data, stored in `adata.obs`. The default is adata.obs [method used].

    output_dir : string, optional  
        Path to output directory.


Returns:

    adata : AnnData Object
        Returns an updated `anndata` object with a new column. check- adata.obs [method used]
        
Example:

```python
    adata = sm.tl.cluster (adata, k= 10, method = 'kmeans', 
    sub_cluster_column='phenotype', use_raw = True)
```

    """

    # Load the andata object    
    if isinstance(adata, str):
        imid = str(adata.rsplit('/', 1)[-1])
        adata = anndata.read(adata)
    else:
        adata = adata
        
    # dynamically adapt the number of neighbours
    if nearest_neighbors > adata.shape[0]:
        nearest_neighbors = adata.shape[0] - 3
    
    
    # Leiden clustering
    def leiden_clustering (pheno, adata, nearest_neighbors, n_pcs, resolution):
        
        # subset the data to be clustered
        if pheno is not None:
            cell_subset =  adata.obs[adata.obs[sub_cluster_column] == pheno].index
        else:
            cell_subset = adata.obs.index
        
        if use_raw == True:
            data_subset = adata[cell_subset]
            if log is True:
                data_subset.X = np.log1p(data_subset.raw.X)          
            else:
                data_subset.X = data_subset.raw.X
        else:
            data_subset = adata[cell_subset]
        
        # clustering
        if pheno is not None:
            print('Leiden clustering ' + str(pheno))
        else:
            print('Leiden clustering')
            
        sc.tl.pca(data_subset)
        if n_pcs is None:
            n_pcs = len(adata.var)
        sc.pp.neighbors(data_subset, n_neighbors=nearest_neighbors, n_pcs=n_pcs)
        sc.tl.leiden(data_subset,resolution=resolution, random_state=random_state)
        
        # Rename the labels
        cluster_labels = list(map(str,list(data_subset.obs['leiden'])))
        if pheno is not None:
            cluster_labels = list(map(lambda orig_string: pheno + '-' + orig_string, cluster_labels))

        # Make it into a dataframe
        cluster_labels = pd.DataFrame(cluster_labels, index = data_subset.obs.index)
        
        # return labels
        return cluster_labels

    # Kmeans clustering
    def k_clustering (pheno, adata, k, sub_cluster_column, use_raw, random_state):
        
        # subset the data to be clustered
        if pheno is not None:
            cell_subset =  adata.obs[adata.obs[sub_cluster_column] == pheno].index
        else:
            cell_subset = adata.obs.index
        
        # Usage of scaled or raw data
        if use_raw == True:
            if log is True:
                data_subset = pd.DataFrame(np.log1p(adata.raw[cell_subset].X), columns =adata[cell_subset].var.index, index = adata[cell_subset].obs.index)
            else:
                data_subset = pd.DataFrame(adata.raw[cell_subset].X, columns =adata[cell_subset].var.index, index = adata[cell_subset].obs.index)
        else:
            data_subset = pd.DataFrame(adata[cell_subset].X, columns =adata[cell_subset].var.index, index = adata[cell_subset].obs.index)
            
        # K-means clustering
        if pheno is not None:
            print('Kmeans clustering ' + str(pheno))
        else:
            print('Kmeans clustering')
        
        kmeans = KMeans(n_clusters=k, random_state=random_state).fit(data_subset)
        
        # Rename the labels
        cluster_labels = list(map(str,kmeans.labels_))
        if pheno is not None:
            cluster_labels = list(map(lambda orig_string: pheno + '-' + orig_string, cluster_labels))
                
        # Make it into a 
        cluster_labels = pd.DataFrame(cluster_labels, index = data_subset.index)
        
        # return labels
        return cluster_labels
        
    # Phenograph clustering
    def phenograph_clustering (pheno, adata, primary_metric, nearest_neighbors):
        
        # subset the data to be clustered
        if pheno is not None:
            cell_subset =  adata.obs[adata.obs[sub_cluster_column] == pheno].index
        else:
            cell_subset = adata.obs.index
        
        # Usage of scaled or raw data
        if use_raw == True:
            data_subset = adata[cell_subset]
            if log is True:
                data_subset.X = np.log1p(data_subset.raw.X)          
            else:
                data_subset.X = data_subset.raw.X
        else:
            data_subset = adata[cell_subset]
        
        # Phenograph clustering
        if pheno is not None:
            print('Phenograph clustering ' + str(pheno))
        else:
            print('Phenograph clustering')
        
        sc.tl.pca(data_subset)
        result = sce.tl.phenograph(data_subset.obsm['X_pca'], k = nearest_neighbors, primary_metric=phenograph_clustering_metric)
        
        # Rename the labels
        cluster_labels = list(map(str,result[0]))
        if pheno is not None:
            cluster_labels = list(map(lambda orig_string: pheno + '-' + orig_string, cluster_labels))
        
        # Make it into a dataframe
        cluster_labels = pd.DataFrame(cluster_labels, index = data_subset.obs.index)
        
        # return labels
        return cluster_labels
        
    
    # PARC clustering
    # https://github.com/ShobiStassen/PARC
    def parc_clustering (pheno, adata, random_state,resolution,parc_too_big_factor,parc_small_pop):
        
        # subset the data to be clustered
        if pheno is not None:
            cell_subset =  adata.obs[adata.obs[sub_cluster_column] == pheno].index
        else:
            cell_subset = adata.obs.index
        
        # Usage of scaled or raw data
        if use_raw == True:
            data_subset = adata[cell_subset]
            if log is True:
                data_subset.X = np.log1p(data_subset.raw.X)
            else:
                data_subset.X = data_subset.raw.X      
        else:
            data_subset = adata[cell_subset]
        
        # Phenograph clustering
        if pheno is not None:
            print('Parc clustering ' + str(pheno))
        else:
            print('Parc clustering')
        
        sc.tl.pca(data_subset)
        parc1 = parc.PARC(data_subset.obsm['X_pca'], random_seed=random_state, parc_small_pop = parc_small_pop, resolution_parameter=resolution,parc_too_big_factor=parc_too_big_factor)  
        parc1.run_PARC() # Run Parc
        

        # Rename the labels
        cluster_labels = list(map(str,parc1.labels))
        if pheno is not None:
            cluster_labels = list(map(lambda orig_string: pheno + '-' + orig_string, cluster_labels))
        
        # Make it into a dataframe
        cluster_labels = pd.DataFrame(cluster_labels, index = data_subset.obs.index)
        
        # return labels
        return cluster_labels
    
    # Use user defined genes for clustering
    if subset_genes is not None:
        bdata = adata[:,subset_genes]
        bdata.raw = bdata[:,subset_genes]
    else:
        bdata = adata.copy()
    
    # IF sub-cluster is True
    # What cells to run the clustering on?
    if sub_cluster is True:
        if sub_cluster_group is not None:
            if isinstance(sub_cluster_group, list):
                pheno = sub_cluster_group
            else:
                pheno = [sub_cluster_group]         
        else:
            # Make sure number of clusters is not greater than number of cells available
            if method == 'kmeans':
                pheno = (bdata.obs[sub_cluster_column].value_counts() > k+1).index[bdata.obs[sub_cluster_column].value_counts() > k+1]
            if method == 'phenograph':
                pheno = (bdata.obs[sub_cluster_column].value_counts() > nearest_neighbors+1).index[bdata.obs[sub_cluster_column].value_counts() > nearest_neighbors+1]
            if method == 'leiden':
                pheno = (bdata.obs[sub_cluster_column].value_counts() > 1).index[bdata.obs[sub_cluster_column].value_counts() > 1]
            if method == 'parc':
                pheno = (bdata.obs[sub_cluster_column].value_counts() > 1).index[bdata.obs[sub_cluster_column].value_counts() > 1]
                
        
    # Run the specified method
    if method == 'kmeans':
        if sub_cluster == True:  
            # Apply the Kmeans function
            r_k_clustering = lambda x: k_clustering(pheno=x, adata=bdata, k=k, sub_cluster_column=sub_cluster_column, use_raw=use_raw, random_state=random_state) # Create lamda function 
            all_cluster_labels = list(map(r_k_clustering, pheno)) # Apply function 
        else:
            all_cluster_labels = k_clustering(pheno=None, adata=bdata, k=k, sub_cluster_column=sub_cluster_column, use_raw=use_raw, random_state=random_state)
            
    if method == 'phenograph':
        if sub_cluster == True:
            r_phenograph_clustering = lambda x: phenograph_clustering(pheno=x, adata=bdata, primary_metric=phenograph_clustering_metric, nearest_neighbors=nearest_neighbors) # Create lamda function 
            all_cluster_labels = list(map(r_phenograph_clustering, pheno)) # Apply function      
        else:
            all_cluster_labels = phenograph_clustering(pheno=None, adata=bdata, primary_metric=phenograph_clustering_metric, nearest_neighbors=nearest_neighbors)
            
            
    if method == 'leiden':
        if sub_cluster == True:
            r_leiden_clustering = lambda x: leiden_clustering(pheno=x, adata=bdata, nearest_neighbors=nearest_neighbors, n_pcs=n_pcs, resolution=resolution) # Create lamda function 
            all_cluster_labels = list(map(r_leiden_clustering, pheno)) # Apply function 
        else:
            all_cluster_labels = leiden_clustering(pheno=None, adata=bdata, nearest_neighbors=nearest_neighbors, n_pcs=n_pcs, resolution=resolution)
            
            
    if method == 'parc':
        if sub_cluster == True:
            r_parc_clustering = lambda x: parc_clustering(pheno=x, adata=bdata, random_state=random_state,resolution=resolution,parc_too_big_factor=parc_too_big_factor,parc_small_pop=parc_small_pop) # Create lamda function 
            all_cluster_labels = list(map(r_parc_clustering, pheno)) # Apply function 
        else:
            all_cluster_labels = parc_clustering(pheno=None, adata=bdata, random_state=random_state,resolution=resolution,parc_too_big_factor=parc_too_big_factor,parc_small_pop=parc_small_pop)
       
    
    # Merge all the labels into one and add to adata
    if sub_cluster == True:
        sub_clusters = pd.concat(all_cluster_labels, axis=0, sort=False)
    else:
        sub_clusters = all_cluster_labels
        
    # Merge with all cells
    #sub_clusters = pd.DataFrame(bdata.obs[sub_cluster_column]).merge(sub_clusters, how='outer', left_index=True, right_index=True)
    sub_clusters = pd.DataFrame(bdata.obs).merge(sub_clusters, how='outer', left_index=True, right_index=True)


    # Transfer labels
    if collapse_labels is False and sub_cluster is True:
        sub_clusters = pd.DataFrame(sub_clusters[0].fillna(sub_clusters[sub_cluster_column]))
    
        
    # Get only the required column
    sub_clusters = sub_clusters[0]
    
    # re index the rows
    sub_clusters = sub_clusters.reindex(adata.obs.index)
    
    # Append to adata
    if label is None:
        adata.obs[method] = sub_clusters
    else:
        adata.obs[label] = sub_clusters
    
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