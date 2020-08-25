**scimap.tl.spatial_expression**

!!! note "Function Call"
    `scimap.tl.spatial_expression` (
      **adata,
      x_coordinate='X_centroid',
      y_coordinate='Y_centroid',
      method='radius', 
      radius=30, 
      knn=10, 
      imageid='imageid', 
      use_raw=True,
      subset=None,
      label='spatial_expression'**)

**Short description**

The `spatial_expression` function allows users to compute a proximity based weighted expression scoring. <br>
The function generates a neighbourhood for each cell and computes a score for all markers based on its proximity to cells within it's neighbourhood. 

The function supports two methods to define a local neighbourhood <br>
* **Radius method**: Can be used to identifies the neighbours within a user defined radius for every cell.
* **KNN method**: Can be used to identifies the neighbours based on K nearest neigbours for every cell

The resultant proportion matrix is saved with `adata.uns`. This can be further clustered to identify similar neighbourhoods.


**Parameters**

`adata` : AnnData Object  

`x_coordinate` : float, required *(The default is 'X_centroid')*  
Column name containing the x-coordinates values.  

`y_coordinate` : float, required *(The default is 'Y_centroid')*  
Column name containing the y-coordinates values.

`method` : string, optional *(The default is 'radius')*  
Two options are available: a) 'radius', b) 'knn'.
a) radius - Identifies the neighbours within a given radius for every cell.
b) knn - Identifies the K nearest neigbours for every cell.

`radius` : int, optional *(The default is 30)*  
The radius used to define a local neighbhourhood.

`knn` : int, optional *(The default is 10)*  
Number of cells considered for defining the local neighbhourhood.

`use_raw` : boolina, optional *(The default is True)*  
Argument to denote whether to use the raw data or scaled data (`sm.pp.rescale`).

`imageid` : string, optional *(The default is 'imageid')*  
Column name of the column containing the image id.

`subset` : string, optional *(The default is None)*  
imageid of the image to be subsetted for analyis. 

`label` : string, optional *(The default is 'spatial_expression')*  
Key for the returned data, stored in `adata.uns`. 


**Returns**
`AnnData` object with the results stored in `adata.uns['spatial_expression']`.

**Example**

```
# Running the radius method
adata = sm.tl.spatial_expression (adata, x_coordinate='X_centroid',y_coordinate='Y_centroid',
                                  method='radius', radius=30, imageid='imageid', 
                                  use_raw=True,subset=None,label='spatial_expression_radius')
# Running the knn method
adata = sm.tl.spatial_expression (adata, x_coordinate='X_centroid',y_coordinate='Y_centroid',
                                  method='knn', knn=10, imageid='imageid', 
                                  use_raw=True,subset=None,label='spatial_expression_knn')
```
