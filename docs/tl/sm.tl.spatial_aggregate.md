**scimap.tl.spatial_aggregate**

!!! note "Function Call"
    `scimap.tl.spatial_aggregate` (
      **adata, 
      x_coordinate='X_centroid',
      y_coordinate='Y_centroid',
      purity = 60, 
      phenotype='phenotype', 
      method='radius', 
      radius=30, 
      knn=10, 
      imageid='imageid',
      subset=None,
      label='spatial_aggregate'**)

**Short description**

The `spatial_aggregate` function allows users to find regions of aggregration of similar cells <br>
The `purity` parameter can be used to tune the granulatity of allowed cell-type heterogenity within local neighbourhood.

The function supports two methods to define a local neighbourhood <br>
**Radius method**: Can be used to identifies the neighbours within a user defined radius for every cell.
**KNN method**: Can be used to identifies the neighbours based on K nearest neigbours for every cell

The resultant proportion matrix is saved with `adata.obs`.


**Parameters**

`adata` : AnnData Object  

`x_coordinate` : float, required *(The default is 'X_centroid')*  
Column name containing the x-coordinates values.  

`y_coordinate` : float, required *(The default is 'Y_centroid')*  
Column name containing the y-coordinates values.

purity : int, required *(The default is 60)* 
Supply a value between 1 and 100. It is the percent purity of neighbouring cells. <br>
For e.g. if 60 is chosen, every neighbourhood is tested such that if a particular phenotype makes up greater than 60% of the total population it is annotated to be an aggregate of that particular phenotype. The default is 60.

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

`label` : string, optional *(The default is 'spatial_aggregate')*  
Key for the returned data, stored in `adata.obs`. 


**Returns**
`AnnData` object with the results stored in `adata.obs['spatial_aggregate']`.

**Example**

```
# Running the radius method
adata = sm.tl.spatial_aggregate (adata, x_coordinate='X_centroid',y_coordinate='Y_centroid',
                                 phenotype='phenotype', method='radius', radius=30,
                                 imageid='imageid',subset=None,label='spatial_aggregate_radius')
# Running the knn method
adata =  sm.tl.spatial_aggregate (adata, x_coordinate='X_centroid',y_coordinate='Y_centroid',
                                  phenotype='phenotype', method='knn', knn=10, 
                                  imageid='imageid',subset=None,label='spatial_aggregate_knn')
```
