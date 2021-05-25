**scimap.tl.spatial_distance**

!!! note "Function Call"
    `scimap.tl.spatial_distance` (
      **adata,
      x_coordinate='X_centroid',
      y_coordinate='Y_centroid',
      phenotype='phenotype',subset=None,
      imageid='imageid',
      label='spatial_distance'**)
      

**Short description**

The `spatial_distance` function enables users to compute the shortest distance
between all cells and every phenotype defined in the dataset. This can be used
to understand the average distance between two cell-types of interest. For
example the average distance between Tumor cell and activated T-cells.
In the event of multiple images present within the dataset, one can compare and
contrast these distances to understand the distribution of cellular proximity in
various biological situations. <br>


The resultant distance matrix is saved with `adata.uns`.


**Parameters**

`adata` : AnnData Object  

`x_coordinate` : float, required *(The default is 'X_centroid')*  
Column name containing the x-coordinates values.  

`y_coordinate` : float, required *(The default is 'Y_centroid')*  
Column name containing the y-coordinates values.

`phenotype` : string, required *(The default is 'phenotype')*  
Column name of the column containing the phenotype information. It could also be any categorical assignment given to single cells. 

`imageid` : string, optional *(The default is 'imageid')*  
Column name of the column containing the image id.

`subset` : string, optional *(The default is None)*  
imageid of the image to be subsetted for analyis. 

`label` : string, optional *(The default is 'spatial_distance')*  
Key for the returned data, stored in `adata.uns`. 


**Returns**
`AnnData` object with the results stored in `adata.uns['spatial_distance']`.


**Example**

```
# Spatial distance
adata = sm.tl.spatial_distance (adata,x_coordinate='X_position',y_coordinate='Y_position',imageid='ImageId')
```
