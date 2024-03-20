# 👑 Additional Helper Function to make Your Life Easier


```python
# import packages
import scimap as sm
import anndata as ad
```

    Running SCIMAP  1.3.14



```python
# Load the data that we saved in the last tutorial (with ROIs added)
adata = ad.read_h5ad('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/resources/exemplarData/scimapExampleData/scimapExampleData.h5ad')
```

### classify

The `sm.hl.classify` function allows users to annotate cells based on the presence or absence of certain markers, providing the option to apply classifications across the entire dataset or within specific subsets, such as groups of cells that have already been phenotyped or clustered. This functionality is especially useful for quickly determining the percentage of cells expressing a particular marker within a subset of interest. A prerequisite for using this function is that gating has been performed and the rescale function has been applied, as it relies on threshold-based classification.


```python
adata.var.index
```




    Index(['ELANE', 'CD57', 'CD45', 'CD11B', 'SMA', 'CD16', 'ECAD', 'FOXP3',
           'NCAM'],
          dtype='object')




```python
# I am going to find out how many cells are CD45 and FOXP3 positive and ECAD negative likely indicating Tregs
adata = sm.hl.classify(adata, pos=['CD45', 'FOXP3'], neg=['ECAD'], collapse_failed=False, label='T_cell_classification')
```


```python
# let's look at the results
adata.obs['T_cell_classification'].value_counts()
```


```python

```

### dropFeatures

The `sm.hl.dropFeatures` function simplifies the refinement of an adata object by allowing users to selectively exclude markers, cells, metadata columns, and particular cell groups.


```python
adata
```




    AnnData object with n_obs × n_vars = 11201 × 9
        obs: 'X_centroid', 'Y_centroid', 'Area', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', 'Solidity', 'Extent', 'Orientation', 'CellID', 'imageid', 'leiden', 'leiden_phenotype', 'ROI', 'phenotype', 'spatial_pscore', 'index_info', 'neigh_kmeans', 'RCNs', 'spatial_lda_kmeans', 'spatial_expression_kmeans', 'spatial_aggregate_radius', 'tumor_similarity_ROI1'
        uns: 'all_markers', 'foldchange_fc', 'foldchange_pval', 'gates', 'spatial_count', 'spatial_distance', 'spatial_expression', 'spatial_interaction_radius', 'spatial_interaction_radius_roi', 'spatial_lda', 'spatial_lda_probability', 'spatial_pscore', 'tumor_similarity'
        obsm: 'umap'
        layers: 'log'



The dataset now contains 11201 cella and 9 markers


```python
# Lets drop 2 markers
adata = sm.hl.dropFeatures(adata, drop_markers=['CD45', 'FOXP3'])
```


```python
adata
```




    AnnData object with n_obs × n_vars = 11201 × 7
        obs: 'X_centroid', 'Y_centroid', 'Area', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', 'Solidity', 'Extent', 'Orientation', 'CellID', 'imageid', 'leiden', 'leiden_phenotype', 'ROI', 'phenotype', 'spatial_pscore', 'index_info', 'neigh_kmeans', 'RCNs', 'spatial_lda_kmeans', 'spatial_expression_kmeans', 'spatial_aggregate_radius', 'tumor_similarity_ROI1'
        uns: 'all_markers', 'foldchange_fc', 'foldchange_pval', 'gates', 'spatial_count', 'spatial_distance', 'spatial_expression', 'spatial_interaction_radius', 'spatial_interaction_radius_roi', 'spatial_lda', 'spatial_lda_probability', 'spatial_pscore', 'tumor_similarity'
        obsm: 'umap'
        layers: 'log'



As you can see now the dataset contains only 7 markers


```python
# lets also drop some cells
adata = sm.hl.dropFeatures(adata, drop_groups='ROI3', groups_column='ROI')
```


```python
adata
```




    AnnData object with n_obs × n_vars = 9629 × 7
        obs: 'X_centroid', 'Y_centroid', 'Area', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', 'Solidity', 'Extent', 'Orientation', 'CellID', 'imageid', 'leiden', 'leiden_phenotype', 'ROI', 'phenotype', 'spatial_pscore', 'index_info', 'neigh_kmeans', 'RCNs', 'spatial_lda_kmeans', 'spatial_expression_kmeans', 'spatial_aggregate_radius', 'tumor_similarity_ROI1'
        uns: 'all_markers', 'foldchange_fc', 'foldchange_pval', 'gates', 'spatial_count', 'spatial_distance', 'spatial_expression', 'spatial_interaction_radius', 'spatial_interaction_radius_roi', 'spatial_lda', 'spatial_lda_probability', 'spatial_pscore', 'tumor_similarity'
        obsm: 'umap'
        layers: 'log'



As you can see now the dataset contains only 9629 cells now


```python

```
