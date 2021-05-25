**scimap.tl.cluster**

!!! note "Function Call"
    `scimap.tl.cluster` (
      **adata,
      method='kmeans', 
      subset_genes=None,
      sub_cluster=False, 
      sub_cluster_column='phenotype', 
      sub_cluster_group=None,
      parc_small_pop=50, 
      parc_too_big_factor=0.4, 
      k=10, 
      n_pcs=None, 
      resolution=1, 
      phenograph_clustering_metric='euclidean', 
      nearest_neighbors=30, 
      use_raw=True, 
      random_state=0, 
      collapse_labels=False**)

**Short description**

The `cluster` function allows users to cluster single-cell data. <br>
The function currently supports four clustering algorithms- kmeans, phenograph, leiden and parc.

The function also allows users to sub-cluster existing clusters by setting `sub_cluster=True`. Check arguments
`sub_cluster_column` and `sub_cluster_group` for more information. <br>

Additionally, if the user wishes to use only a subset of genes for the purpose of clustering, it can be
acheived by passing the genes as a list to `subset_genes`.

The resultant clusters are saved under `adata.obs[method used]`.


**Parameters**

`adata` : AnnData Object  

`method` : string, optional *(The default is 'kmeans')*  
Clustering method to be used- Implemented methods- kmeans, phenograph, leiden and parc.

`subset_genes` : list, optional *(The default is None)*  
Pass a list of genes ['CD3D', 'CD20', 'KI67'] that should be included for the purpose of clustering. <br> 
By default the algorithm uses all genes in the dataset.

`sub_cluster` : bool, optional *(The default is False)*  
If the user has already performed clustering or phenotyping previously and would like to sub-cluster within a particular cluster/phenotype, this option can be used.

`sub_cluster_column` : string, optional *(The default is 'phenotype')*  
The column name that contains the cluster/phenotype information to be sub-clustered. This is only required when sub_cluster is set to True.

`sub_cluster_group` : list, optional *(The default is None)*  
By default the program will sub-cluster all groups within column passed through the argument sub_cluster_column. If user wants to sub cluster only a subset of phenotypes/clusters this option can be used. Pass them as list e.g. ["tumor", "b cells"].  

`parc_small_pop` : int, optional *(The default is 50)*  
Smallest cluster population to be considered a community in PARC clustering.

`parc_too_big_factor` : float, optional *(The default is 0.4)*  
If a cluster exceeds this share of the entire cell population, then the PARC will be run on the large cluster. at 0.4 it does not come into play.

`k` : int, optional *(The default is 10)*  
Number of clusters to return when using K-Means clustering.

`n_pcs` : int, optional *(The default is None)*  
 Number of PC's to be used in leiden clustering. By default it uses all PC's.

`resolution` : float, optional *(The default is 1)*  
A parameter value controlling the coarseness of the clustering. Higher values lead to more clusters.

`phenograph_clustering_metric` : string, optional *(The default is 'euclidean')*   
Distance metric to define nearest neighbors. Note that performance will be slower for correlation and cosine. <br>
Available methods- cityblock’, ‘cosine’, ‘euclidean’, ‘manhattan’, braycurtis’, ‘canberra’, ‘chebyshev’,  ‘correlation’, ‘dice’, ‘hamming’, ‘jaccard’, ‘kulsinski’, ‘mahalanobis’, ‘minkowski’, ‘rogerstanimoto’,  ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘yule’

`nearest_neighbors` : int, optional *(The default is 30)*  
Number of nearest neighbors to use in first step of graph construction. This parameter is used both in `leiden` and `phenograph` clustering.

`use_raw` : bool, optional *(The default is True)*  
If True, log transformed raw data will be used for clustering. If False, normalized/scaled data within `adata.X` will be used.

`random_state` : int, optional *(The default is 0)*  
Change the initialization of the optimization. 

`collapse_labels` : bool, optional *(The default is False)*  
While sub clustering only a few phenotypes/clusters, this argument helps to group all the other phenotypes/clusters into a single category -  Helps in visualisation.

`label` : string, optional *(The default is None)*  
Key or optional column name for the returned data, stored in `adata.obs`. The default is adata.obs[method used].


**Returns**
`AnnData` object with the results stored in `adata.obs[method used]`.

**Example**

```
# Running clustering on entire data using K-means method
adata = sm.tl.cluster (adata,  method = 'kmeans', k= 10, use_raw = True)

# Sub-cluster a already named cluster called `Tumor` using leiden clustering
adata = sm.tl.cluster (adata, method = 'leiden', resolution = 0.5, 
        nearest_neighbors = 20, use_raw = True,
        sub_cluster=True, sub_cluster_column='phenotype', sub_cluster_group='Tumor')
        
# Run phenograph clustering by only using a subset of genes
gene_subset = ['CD25', 'CD2', 'CD10', 'CD163', 'CD3D', 'CD5', 'CD30', 'ACTIN', 'CD45', 
                'CD206', 'CD68', 'PD1', 'KI67', 'CD11C', 'CD7', 'CD8A', 'FOXP3', 'CD20']
adata = sm.tl.cluster (adata,  subset_genes = gene_subset, method = 'phenograph', 
        nearest_neighbors = 10, use_raw = True)


```
