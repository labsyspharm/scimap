**scimap.tl.spatial_count**

!!! note "Function Call"
    `scimap.tl.spatial_count` (
      **adata,
      x_coordinate='X_centroid',
      y_coordinate='Y_centroid',
      phenotype='phenotype',
      method='radius',
      radius=30,
      knn=10,
      imageid='imageid',
      subset=None,
      label='spatial_count'**)

**Short description**

The `spatial_count` function allows users to compute the proportion of a categorical variable (e.g. cell-types) within the local neighbourhood of each cell. 

The function supports two methods to define a local neighbourhood <br>
**Radius method**: Can be used to identifies the neighbours within a user defined radius for every cell.
**KNN method**: Can be used to identifies the neighbours based on K nearest neigbours for every cell

The resultant proportion matrix is saved with `adata.uns`. This can be further clustered to identify similar neighbourhoods.


**Parameters**

`adata` : AnnData Object  

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

`radius` : int, optional *(The default is 30)*  
The radius used to define a local neighbhourhood.

`knn` : int, optional *(The default is 10)*  
Number of cells considered for defining the local neighbhourhood.

`imageid` : string, optional *(The default is 'imageid')*  
Column name of the column containing the image id.

`subset` : string, optional *(The default is None)*  
imageid of the image to be subsetted for analyis. 

`label` : string, optional *(The default is 'spatial_count')*  
Key for the returned data, stored in `adata.uns`. 


**Returns**
`AnnData` object with the results stored in `adata.uns['spatial_count']`.

**Example**

```
# Running the radius method
adata = sm.tl.spatial_count (adata,x_coordinate='X_centroid',y_coordinate='Y_centroid',
                       phenotype='phenotype',method='radius',radius=30,
                       imageid='imageid',subset=None,label='spatial_count_radius')
                           
# Running the knn method
adata = sm.tl.spatial_count (adata,x_coordinate='X_centroid',y_coordinate='Y_centroid',
                       phenotype='phenotype',method='knn',knn=10,
                       imageid='imageid',subset=None,label='spatial_count_knn')
```
