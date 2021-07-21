**scimap.tl.spatial_lda**

!!! note "Function Call"
    `scimap.tl.spatial_lda` (
      **adata, 
      x_coordinate='X_centroid', 
      y_coordinate='Y_centroid', 
      phenotype='phenotype', 
      method='radius', radius=30, knn=10, 
      imageid='imageid', num_motifs=10, 
      random_state=0, subset=None, 
      label='spatial_lda', 
      **kwargs**)

**Short description**

The `spatial_lda` function applies Latent Dirichlet Allocation (LDA) method to identify spatial motifs. First local 
neighbourhoods are defined and then they are passed togeather through the LDA modelling to identify motifs of neighbourhoods.<br>

The function supports two methods to define a local neighbourhood <br>
**Radius method**: Can be used to identifies the neighbours within a user defined radius for every cell.
**KNN method**: Can be used to identifies the neighbours based on K nearest neigbours for every cell

The results of the method are stored in `adata.uns`<br>
The weights in latent space for each cell which is further used for clustering to identify the motifs are
saved under `adata.uns['spatial_lda']`<br>
The weights each cell-types contribution to the spatial motif is stored under `adata.uns['spatial_lda_probability']`<br>
The model itself is saved under `adata.uns['spatial_lda_model']`


**Parameters**

`adata` : AnnData Object  

`x_coordinate` : float, required *(The default is 'X_centroid')*  
Column name containing the x-coordinates values.  

`y_coordinate` : float, required *(The default is 'Y_centroid')*  
Column name containing the y-coordinates values.

`num_motifs` : int, optional *(The default is 10)*  
The number of requested latent motifs to be extracted from the training corpus. 

`method` : string, optional *(The default is 'radius')*  
Two options are available: a) 'radius', b) 'knn'.
a) radius - Identifies the neighbours within a given radius for every cell.
b) knn - Identifies the K nearest neigbours for every cell.

`radius` : int, optional *(The default is 30)*  
The radius used to define a local neighbhourhood.

`knn` : int, optional *(The default is 10)*  
Number of cells considered for defining the local neighbhourhood.

`imageid` : string, optional *(The default is 'imageid')*  
Column name of the column containing the image id.

`subset` : string, optional *(The default is None)*  
imageid of the image to be subsetted for analyis. 

`random_state` : int, optional *(The default is 0)*  
Either a randomState object or a seed to generate one. Useful for reproducibility.

`label` : string, optional *(The default is 'spatial_lda')*  
Key for the returned data, stored in `adata.uns`. 


**Returns**
`AnnData` object with the results stored in `adata.uns['spatial_lda']`.

**Example**

```
# Running the radius method
adata = spatial_lda (adata, num_motifs=10, radius=100)

```
