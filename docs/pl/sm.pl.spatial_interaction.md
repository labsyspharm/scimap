**scimap.pl.spatial_interaction**

!!! note "Function Call"
    `scimap.pl.spatial_interaction` (
      **adata, 
      spatial_interaction='spatial_interaction',
      summarize_plot=True, p_val=0.05, nonsig_color='grey',
      row_cluster=False, col_cluster=False,
      cmap = 'vlag', **kwargs**)

**Short description**

The function allows users to generate a heatmap to visualize spatial interaction output from
`sm.tl.spatial_interaction`. The intensity represents number of interactions (scaled) 
and blank regions represent non-significant results.

**Parameters**

`adata` : AnnData object  

`spatial_interaction`: string *(The default is 'spatial_interaction')*  
In order to locate the spatial_interaction data within the AnnData object please provide the output 
label/columnname of `sm.tl.spatial_interaction` function. The default is 'spatial_interaction'.

`summarize_plot` : bool, optional *(The default is True)*   
In the event of analyzing multiple images, this argument allows users to 
plot the average cell-cell interaction across all images. The default is True.

`p_val` : float, optional *(The default is 0.05.)*  
P-value cut-off above which interactions are not considered significant. 

`nonsig_color` : string, optional *(The default is 'grey')*  
Color for non-significant interactions (Interactions above the P-value cut-off will use this color).

`row_cluster` : bool, optional *(The default is False)*  
Cluster Rows.

`col_cluster` : boo,, optional *(The default is False)*  
Cluster Columns. 

`cmap` : string, optional *(The default is `'vlag'`)*  
Color map to use for continous variables. Can be a name or a Colormap 
instance (e.g. `'magma'`, `'viridis'`). 

`**kwargs`: key:value pairs.  
Are passed to sns.clustermap. Pass other parameters that works with `sns.clustermap`. e.g. `linecolor='black'`
 

**Returns**

Heatmap

**Example**

```
# spatial_interaction heatmap for a single image
sm.pl.spatial_interaction(adata, summarize_plot=True, row_cluster=True, linewidths=0.75, linecolor='black')
    
# spatial_interaction heatmap for multiple images
sns.set(font_scale=0.6)
m.pl.spatial_interaction(adata, summarize_plot=False, row_cluster=True, col_cluster=True, yticklabels=True)
```
