#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 17:03:56 2020
@author: Ajit Johnson Nirmal
Function to sub-cluster a cluster of interest. Particularly useful to check if 
there are sub-phenotypes after performing the gating based phenotyping.
"""

# Import library
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import parc
from sklearn.cluster import KMeans


def cluster (adata, cluster_group = None, small_pop= 50, too_big_factor=0.4, 
                          group_others= False, k= 10, n_pcs=None, method = 'kmeans', resolution=1, 
                          clustering_metric='euclidean', nearest_neighbors= 30, label='phenotype', 
                          use_raw_data = True, random_state=0, genes=None):
    """
    
    cluster any data
    cluster within groups
    cluster within specific group
    
    
    
    
    Parameters
    ----------
    adata : AnnData Object
    cluster_group : list, optional
    The phenotype/cell_types to sub cluster. Pass values as list e.g. ["tumor", "b cells"]. 
    By default it will run on all groups within column passed through the argument label. The default is None.
    group_others : bool, optional
    While sub clustering only a few phenotypes/ cell types, this argument helps to group all the other cell types into a single category- Helps in visualisation. The default is False.
    label : string, optional
    The column name that contains the phenotype information. The default is 'phenotype'.
    k : int, optional
        Number of clusters if KMeans clustering method is used. The default is 10.
    method : string, optional
    Clustering method to be used- Implemented methods- kmeans, phenograph, leiden, parc. The default is 'kmeans'.
    n_pcs : int, optional
        Number of PC's to use in leiden clustering. The default is 20.
    resolution : int, optional
        A parameter value controlling the coarseness of the clustering. Higher values lead to more clusters. The default is 1.
    clustering_metric : string, optional
        Distance metric to define nearest neighbors. Note that performance will be slower for correlation and cosine. 
        Available methods- cityblock’, ‘cosine’, ‘euclidean’, ‘manhattan’, braycurtis’, ‘canberra’, ‘chebyshev’, 
        ‘correlation’, ‘dice’, ‘hamming’, ‘jaccard’, ‘kulsinski’, ‘mahalanobis’, ‘minkowski’, ‘rogerstanimoto’, 
        ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘yule’
        The default is 'euclidean'.
    nearest_neighbors : int, optional
        Number of nearest neighbors to use in first step of graph construction. The default is 30.
    use_raw_data : bool, optional
        If True, log transformed raw data will be used for clustering. If False, normalized data will be used. The default is True.
    random_state : int, optional
        Change the initialization of the optimization. The default is 0.

    Returns
    -------
    adata : AnnData Object
        Returns an updated anndata object with a new column. check- adata.obs[method used]
        
    Example
    -------
    adata = phenotype_subcluster (adata, k= 10, method = 'kmeans', label='phenotype', use_raw_data = True)

    """
    
    # Leiden clustering
    def leiden_clustering (pheno, adata, nearest_neighbors, n_pcs, resolution):
        
        # subset the data to be clustered
        if pheno is not None:
            cell_subset =  adata.obs[adata.obs[label] == pheno].index
        else:
            cell_subset = adata.obs.index
        
        if use_raw_data == True:
            data_subset = adata[cell_subset]
            data_subset.X = np.log1p(data_subset.raw.X)          
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
    def k_clustering (pheno, adata, k, label, use_raw_data, random_state):
        
        # subset the data to be clustered
        if pheno is not None:
            cell_subset =  adata.obs[adata.obs[label] == pheno].index
        else:
            cell_subset = adata.obs.index
        
        # Usage of scaled or raw data
        if use_raw_data == True:
            data_subset = pd.DataFrame(np.log1p(adata.raw[cell_subset].X), columns =adata[cell_subset].var.index, index = adata[cell_subset].obs.index)
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
            cell_subset =  adata.obs[adata.obs[label] == pheno].index
        else:
            cell_subset = adata.obs.index
        
        # Usage of scaled or raw data
        if use_raw_data == True:
            data_subset = adata[cell_subset]
            data_subset.X = np.log1p(data_subset.raw.X)          
        else:
            data_subset = adata[cell_subset]
        
        # Phenograph clustering
        print('Phenograph clustering ' + str(pheno))
        sc.tl.pca(data_subset)
        result = sce.tl.phenograph(data_subset.obsm['X_pca'], k = nearest_neighbors, primary_metric=clustering_metric)
        
        # Rename the labels
        cluster_labels = list(map(str,result[0]))
        cluster_labels = list(map(lambda orig_string: pheno + '-' + orig_string, cluster_labels))
        
        # Make it into a dataframe
        cluster_labels = pd.DataFrame(cluster_labels, index = data_subset.obs.index)
        
        # return labels
        return cluster_labels
        
    
    # PARC clustering
    # https://github.com/ShobiStassen/PARC
    def parc_clustering (pheno, adata, random_state,resolution,too_big_factor,small_pop):
        
        # subset the data to be clustered
        if pheno is not None:
            cell_subset =  adata.obs[adata.obs[label] == pheno].index
        else:
            cell_subset = adata.obs.index
        
        # Usage of scaled or raw data
        if use_raw_data == True:
            data_subset = adata[cell_subset]
            data_subset.X = np.log1p(data_subset.raw.X)          
        else:
            data_subset = adata[cell_subset]
        
        # Phenograph clustering
        if pheno is not None:
            print('Parc clustering ' + str(pheno))
        else:
            print('Parc clustering')
        
        sc.tl.pca(data_subset)
        parc1 = parc.PARC(data_subset.obsm['X_pca'], random_seed=random_state, small_pop = small_pop, resolution_parameter=resolution,too_big_factor=too_big_factor)  
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
    if genes is not None:
        bdata = adata[:,genes]
        bdata.raw = bdata[:,genes]
    else:
        bdata = adata
        
    # What cells to run the clustering on?
    if cluster_group is not None:
        pheno = list(cluster_group)        
    else:
        pheno = bdata.obs[label].unique()
    
    # Run the specified method
    if method == 'kmeans':
        # Apply the Kmeans function
        r_k_clustering = lambda x: k_clustering(pheno=x, adata=bdata, k=k, label=label, use_raw_data=use_raw_data, random_state=random_state) # Create lamda function 
        all_cluster_labels = list(map(r_k_clustering, pheno)) # Apply function 
    if method == 'phenograph':
        r_phenograph_clustering = lambda x: phenograph_clustering(pheno=x, adata=bdata, primary_metric=clustering_metric, nearest_neighbors=nearest_neighbors) # Create lamda function 
        all_cluster_labels = list(map(r_phenograph_clustering, pheno)) # Apply function        
    if method == 'leiden':
        r_leiden_clustering = lambda x: leiden_clustering(pheno=x, adata=bdata, nearest_neighbors=nearest_neighbors, n_pcs=n_pcs, resolution=resolution) # Create lamda function 
        all_cluster_labels = list(map(r_leiden_clustering, pheno)) # Apply function 
    if method == 'parc':
        r_parc_clustering = lambda x: parc_clustering(pheno=x, adata=bdata, random_state=random_state,resolution=resolution,too_big_factor=too_big_factor,small_pop=small_pop) # Create lamda function 
        all_cluster_labels = list(map(r_parc_clustering, pheno)) # Apply function 
       
    # Merge all the labels into one and add to adata
    sub_clusters = pd.concat(all_cluster_labels, axis=0, sort=False)
        
    # Merge with all cells
    sub_clusters = pd.DataFrame(bdata.obs[label]).merge(sub_clusters, how='outer', left_index=True, right_index=True)
        
    
    #cluster_labels.columns = ['KMeans']
    #sub_clusters = pd.DataFrame(adata.obs).merge(cluster_labels, how='outer', left_index=True, right_index=True)

    
    # Transfer labels
    if group_others == False:
        sub_clusters = pd.DataFrame(sub_clusters[0].fillna(sub_clusters[label]))
        
    # Get only the required column
    sub_clusters = sub_clusters[0]
    
    # re index the rows
    sub_clusters = sub_clusters.reindex(adata.obs.index)
    
    # Append to adata
    adata.obs[method] = sub_clusters
    
    # Return adata
    return adata
