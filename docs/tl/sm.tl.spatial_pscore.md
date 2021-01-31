**scimap.tl.spatial_pscore**

!!! note "Function Call"
    `scimap.tl.spatial_pscore` (
      **adata, 
      proximity, 
      score_by='imageid', 
      x_coordinate='X_centroid', 
      y_coordinate='Y_centroid', 
      phenotype='phenotype', 
      method='radius', 
      radius=20, knn=3, 
      imageid='imageid', 
      subset=None, 
      label='spatial_pscore'**)

**Short description**

The `spatial_pscore` function enables users to score interation between cell types.
The function generates two scores and saved at `adata.uns`: <br>
A) Proximity Density: Total number of interactions identified divided by the total number of 
cells of the cell-types that were used for interaction analysis. <br>s
B) Proximity Volume: Total number of interactions identified divided by the total number of all cells in the data.
The interaction sites are also recorded and saved in `adata.obs` <br>

The function supports two methods to define a local neighbourhood <br>
**Radius method**: Can be used to identifies the neighbours within a user defined radius for every cell.
**KNN method**: Can be used to identifies the neighbours based on K nearest neigbours for every cell

The resultant proportion matrix is saved with `adata.uns` and `adata.obs`.


**Parameters**

`adata` : AnnData Object  

`proximity` : list, required  
Pass a list of cell-types for which the proximity score needs to calculated. e.g. ['CellType-A', 'CellType-B'].
        
`score_by` : string, optional *(The default is 'imageid')* 
If the scores need to compared across region's of interest, the column name containing the ROI's 
should be passed. By default the score is calculated across the entire image. The default is 'imageid'.  

`x_coordinate` : float, required *(The default is 'X_centroid')*  
Column name containing the x-coordinates values.  

`y_coordinate` : float, required *(The default is 'Y_centroid')*  
Column name containing the y-coordinates values.

`phenotype` : string, required *(The default is 'phenotype')*  
Column name of the column containing the phenotype information. It could also be any categorical assignment given to single cells. 

`method` : string, optional *(The default is 'radius')*  
Two options are available: a) 'radius', b) 'knn'.
a) radius - Identifies the neighbours within a given radius for every cell.
b) knn - Identifies the K nearest neigbours for every cell.

`radius` : int, optional *(The default is 20)*  
The radius used to define a local neighbhourhood.

`knn` : int, optional *(The default is 3)*  
Number of cells considered for defining the local neighbhourhood.

`imageid` : string, optional *(The default is 'imageid')*  
Column name of the column containing the image id.

`subset` : string, optional *(The default is None)*  
imageid of the image to be subsetted for analyis. 

`label` : string, optional *(The default is 'spatial_pscore')*  
Key for the returned data, stored in `adata.uns` and `adata.obs`. 


**Returns**
`AnnData` object with the results stored in `adata.uns['spatial_pscore']` and `adata.obs['spatial_pscore']`.


**Example**

```
# Calculate the score for proximity between `Tumor CD30+` cells and `M2 Macrophages`
adata =  spatial_pscore (adata,proximity= ['Tumor CD30+', 'M2 Macrophages'], score_by = 'ImageId',
                        x_coordinate='X_position',y_coordinate='Y_position',
                        phenotype='phenotype',method='radius',radius=20,knn=3,
                        imageid='ImageId',subset=None, label='spatial_pscore')
```
