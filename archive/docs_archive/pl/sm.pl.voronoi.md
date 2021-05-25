**scimap.pl.voronoi**

!!! note "Function Call"
    `scimap.pl.voronoi` (
      **adata, 
      color_by=None, colors=None, 
      x_coordinate='X_centroid', y_coordinate='Y_centroid',
      imageid='imageid',subset=None,
      voronoi_edge_color = 'facecolor',
      voronoi_line_width = 0.1, 
      voronoi_alpha = 0.5, 
      size_max=np.inf,
      overlay_points=None, 
      overlay_points_categories=None, 
      overlay_drop_categories=None,
      overlay_points_colors=None,
      overlay_point_size = 5, 
      overlay_point_alpha= 1, 
      overlay_point_shape=".", 
      plot_legend=True, 
      legend_size = 6,
      **kwargs **)

**Short description**

The `sm.pl.voronoi` function allows users to generate a voronoi diagram colored by user defined categories (e.g. cell types or neighbourhood).
The function also allows to overlay a scatter plot on top of the voronoi diagram (any categorical column). 

**Parameters**

`adata` : AnnData object   

`color_by` : string, optional *(The default is None)*  
        Color the voronoi diagram based on categorical variable (e.g. cell types or neighbourhoods).
        Pass the name of the column which contains the categorical variable.
        
`colors` : string or dict, optional *(The default is None)*  
        Custom coloring the voronoi diagram. The parameter accepts `sns color palettes` or a python dictionary
        mapping the categorical variable with the required color.
        
`x_coordinate` : float, optional *(The default is 'X_centroid')*  
        Column name containing the x-coordinates values. 
        
`y_coordinate` : float, optional *(The default is 'Y_centroid')*  
        Column name containing the y-coordinates values.
        
`imageid` : string, optional *(The default is 'imageid')*  
        Column name of the column containing the image id.

`subset` : list, optional *(The default is None)*  
        imageid of a single image to be subsetted for plotting. 

`voronoi_edge_color` : string, optional *(The default is 'black')*  
        A Matplotlib color for marking the edges of the voronoi. 
        If `facecolor` is passed, the edge color will always be the same as the face color.

`voronoi_line_width` : float, optional *(The default is 0.1)*  
        The linewidth of the marker edges. Note: The default edgecolors is 'face'. You may want to change this as well. 

`voronoi_alpha` : float, optional *(The default is 0.5)*  
        The alpha blending value, between 0 (transparent) and 1 (opaque).

`overlay_points` : string, optional *(The default is None)*  
        It is possible to overlay a scatter plot on top of the voronoi diagram.
        Pass the name of the column which contains categorical variable to be overlayed.

`overlay_points_categories` : list, optional *(The default is None)*  
       If the passed column in `overlay_points` contains multiple categories, however the user is only 
       interested in a subset of categories, those specific names can be passed as a list. By default all 
       categories will be overlayed on the voronoi diagram.

`overlay_drop_categories` : list, optional *(The default is None)*  
        Similar to `overlay_points_categories`. Here for ease of use, especially if large number of categories are present.
        The user can drop a set of categories.

`overlay_points_colors` : string or dict, optional *(The default is None)*  
        Similar to `colors`.  
        User can pass in a  
        a) solid color (like `black`)  
        b) sns palettes name (like `Set1`)  
        c) python dictionary mapping the categories with custom colors

`overlay_point_size` : float, optional *(The default is 5)*  
         Overlay scatter plot point size.
        
`overlay_point_alpha` : float, optional *(The default is 1)*  
        The alpha blending value for the overlay, between 0 (transparent) and 1 (opaque).
        
`overlay_point_shape` : string, optional *(The default is '.')*  
        The marker style. marker can be either an instance of the class or the text shorthand for a particular marker.
        
`plot_legend` : bool, optional *(The default is True)*  
        Define if the figure legend should be plotted.
        Please note the figure legend may be out of view and you may need to resize the image to see it, especially 
        the legend for the scatter plot which will be on the left side of the plot.
        
`legend_size` : float, optional *(The default is 6)*  
         Resize the legend if needed.

`**kwargs` : Additional keyword arguments that can be passed to `sns.scatterplot`  


**Example**

```
# Plot a voronoi diagram using the `phenotype` column 
sm.pl.voronoi(adata, color_by='phenotype', colors=None, x_coordinate='X_position', y_coordinate='Y_position',
             imageid='ImageId',subset=None,
             voronoi_edge_color = 'black',voronoi_line_width = 0.2, voronoi_alpha = 0.5, size_max=np.inf,
             overlay_points='phenotype', overlay_points_categories=None, overlay_drop_categories=None,
             overlay_point_size = 5, overlay_point_alpha= 1, overlay_point_shape=".", plot_legend=False, legend_size=6)
             
```
