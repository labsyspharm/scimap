**scimap.pl.spatial_interaction**

!!! note "Function Call"
    `scimap.pl.spatial_distance` (
      **adata, 
      spatial_distance='spatial_distance', 
      phenotype='phenotype', imageid='imageid', 
      method='heatmap', heatmap_summarize=True, 
      heatmap_na_color='grey',heatmap_cmap='vlag_r', 
      heatmap_row_cluster=False, heatmap_col_cluster=False, 
      heatmap_standard_scale=0, distance_from=None, 
      distance_to=None, x_axis = None, y_axis = None, 
      facet_by = None, plot_type = None, 
      **kwargs **)

**Short description**

The `spatial_distance` plotting function allows users to generate a heatmap to visualize spatial interaction output from
`sm.tl.spatial_interaction`. The intensity represents number of interactions (scaled) 
and blank regions represent non-significant results.

**Parameters**

`adata` : AnnData object  

`spatial_distance`: string *(The default is 'spatial_distance')*  
In order to locate the spatial_distance data within the AnnData object please provide the output 
label/columnname of `sm.tl.spatial_distance` function.

`phenotype` : string, required *(The default is 'phenotype')*  
Column name of the column containing the phenotype information. 
It could also be any categorical assignment given to single cells. 

`imageid` : string, optional *(The default is 'imageid')*  
Column name of the column containing the image id.

`method` : string, optional *(The default is 'heatmap')*  
Three options are available.  
1) `heatmap` - generates a heatmap of average shortest distance between all phenotypes.  
2) `numeric` - can be used to generate boxplot, violin plot etc between a given set of phenotypes.  
3) `distribution` - can be used to generate distribution plots between a given set of phenotypes.  

`heatmap_summarize` : bool, optional *(The default is True)*   
In the event multiple images are present in the dataset, True allows to calculate the 
average across all the images.

`heatmap_na_color` : string, optional *(The default is 'grey')*  
Color for NA values within the heatmap.

`heatmap_row_cluster` : bool, optional *(The default is False)*  
Cluster Rows.

`heatmap_col_cluster` : bool, optional *(The default is False)*  
Cluster Columns. 

`heatmap_cmap` : string, optional *(The default is `'vlag'`)*  
Color map to use for continous variables. Can be a name or a Colormap 
instance (e.g. `'magma'`, `'viridis'`). 

`heatmap_standard_scale` : int, optional *(The default is 0)*  
Either 0 (rows) or 1 (columns). Whether or not to standardize that dimension, 
meaning for each row or column, subtract the minimum and divide each by its maximum. The default is 0.
      
`distance_from` : string, optional *(The default is None)*  
In the event of using method = 'numeric' or 'distribution', this argument is required.  
Pass a phenotype of interest. If distance_from is provided and distance_to is not provided,
the function will plot the average distance from the phenotype of interest to all
phenotypes present within the dataset.

`distance_to` : string, optional *(The default is None)*  
In the event of using method = 'numeric' or 'distribution', this argument is required.  
Pass a phenotype of interest. The function will plot the average shortest between two phenotypes of
interest (distance_from and distance_to).

`x_axis` : string, optional *(The default is `distance`)*  
In the event of using method = 'numeric' or 'distribution', this argument is required.
This determines the elements present in the x-axis of the resultant plot.  
Allowed arguments are: 'group', 'distance', 'imageid'.

`y_axis` : string, optional *(The default is `group`)*  
In the event of using method = 'numeric' or 'distribution', this argument is required.
This determines the elements present in the y-axis of the numeric plot and if the user uses the distribution
plot this argument is used to overlaying multiple categories within the same distribution plot.  
Allowed arguments are: 'group', 'distance', 'imageid'.

`facet_by` : string, optional *(The default is None)*  
In the event of using method = 'numeric' or 'distribution', this argument can be used to
generate sub-plots. Allowed arguments are: 'group', 'imageid'.

`plot_type` : string, optional *(The default is None)*  
In the event of using method = 'numeric' or 'distribution', this argument is required.  
- For `numeric` plot, the following options are available: “strip”, “swarm”, “box”, “violin”, “boxen”, “point”, “bar”, or “count”.  
- For `distribution` plot, the following options are available: “hist”, “kde”, “ecdf”.  
1) The default for `numeric` plot is 'boxen'.  
2) The default for `distribution` plot is 'kde`.  

`**kwargs`: key:value pairs.  
Are passed to sns.clustermap. Pass other parameters that works with `sns.clustermap`, `sns.catplot` or `sns.displot`
e.g. `linecolor='black'`
 

**Returns**

Heatmap or Numeric Plot or Distribution Plot.

**Example**

```
# summary heatmap
sm.pl.spatial_distance (adata)
    
# Heatmap without summarizing the individual images
sm.pl.spatial_distance (adata, heatmap_summarize=False, imageid='ImageId')
    
# Numeric plot of shortest distance of phenotypes from tumor cells
sm.pl.spatial_distance (adata, method='numeric',distance_from='Tumor CD30+',imageid='ImageId')
    
# Distribution plot of shortest distance of phenotypes from tumor cells
sm.pl.spatial_distance (adata, method='distribution',distance_from='Tumor CD30+',imageid='ImageId', x_axis="distance", y_axis="imageid", plot_type="kde")
    
# Numeric plot of shortest distance of phenotypes from tumor cells to M2 Macrophages
sm.pl.spatial_distance (adata, method='numeric',distance_from='Tumor CD30+',distance_to = 'M2 Macrophages', imageid='ImageId')
    
# Distribution plot of shortest distance of phenotypes from tumor cells to M2 Macrophages
sm.pl.spatial_distance (adata, method='distribution',distance_from='Tumor CD30+',distance_to = 'M2 Macrophages',imageid='ImageId')
```
