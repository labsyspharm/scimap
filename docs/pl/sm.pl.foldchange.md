**scimap.pl.foldchange**

!!! note "Function Call"
    `scimap.pl.foldchange` (
      **adata, 
      label='foldchange', 
      p_val=0.05, 
      nonsig_color='grey', 
      subset_xaxis=None, 
      subset_yaxis=None, 
      cmap = 'vlag', 
      log=True, 
      center=0, 
      method='heatmap', 
      invert_axis=None, 
      parallel_coordinates_color=None, 
      matplotlib_bbox_to_anchor=(1.04,1), 
      matplotlib_legend_loc='upper left', 
      xticks_rotation=90, 
      return_data = False 
      **kwargs **)

**Short description**

The `sm.pl.foldchange` function allows users to visualize foldchanges in celltypes between samples/ROI's. 
Run `sm.tl.foldchange` first to compute the foldchange. <br>
<br>
The function incorportates two methods.<br>
a) A Heatmap plot - Use `method = 'heatmap'` <br>
b) Parallel Coordinates plot : Use `method = 'parallel_coordinates'`

**Parameters**

`adata` : AnnData object   

`label` : string, optional *(The default is `foldchange`)*  
        label used when running `sm.tl.foldchange`.
        
`p_val` : float, optional *(The default is `0.05`)*  
        p_val cut-off above which is considered not-significant. The cells containing 
        non-significant changes will be highlighted in the heatmap.
        
`method` : string, optional *(The default is `heatmap`)*  
        Two methods are available for plotting the foldchanges  
        a) Heatmap: Use `heatmap`  
        b) parallel coordinates plot : Use `parallel_coordinates`  
        
`nonsig_color` : string, optional *(The default is `grey`)*  
        Color used to highlight non-significant fold changes in the heatmap.
        
`subset_xaxis` : list, optional *(The default is None)*  
        Subset x-axis before plotting. Pass in a list of categories. eg- `subset_xaxis = ['CelltypeA', 'CellTypeB']`.
        
`subset_yaxis` : list, optional *(The default is None)*  
        Subset y-axis before plotting. Pass in a list of categories. eg- `subset_yaxis = ['ROI_1', 'ROI_5']`. 
        
`cmap` : string, optional *(The default is `vlag`)*  
        Color map. Can be a name or a Colormap instance (e.g. 'magma', 'viridis'). 
        
`log` : bool, optional *(The default is `True`)*  
        Convert foldchange to log2 scale. 
        
`center` : float, optional *(The default is `0`)*  
        The center value to be used in heatmap.

`invert_axis` : bool, optional *(The default is `None`)*  
        Flip the axis of the plot.
        
`parallel_coordinates_color` : list, optional *(The default is None)*  
        Custom colors for each category.
        
`xticks_rotation` : int, optional *(The default is `90`)*  
        Angle the x-axis ticks. 
        
`matplotlib_bbox_to_anchor` : tuple, optional *(The default is `(1.04,1)`)*  
        Bounding box argument used along with matplotlib_legend_loc to control
        the legend location when using the matplotlib method.
        
`matplotlib_legend_loc` : int/string, optional *(The default is `upper left`)*  
        Location of legend used along with matplotlib_bbox_to_anchor to control
        the legend location when using the matplotlib method.
        
`return_data` : bool, optional *(The default is False)*  
        When True, return the data used for plotting.
        
`**kwargs` : Additional keyword arguments passed to:  
        1) Pandas [DataFrame.parallel_coordinates())](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.plotting.parallel_coordinates.html) when using the `parallel_coordinates` method.  
        2) [sns.clustermap()](https://seaborn.pydata.org/generated/seaborn.clustermap.html) when using the `heatmap` method.


**Returns**

If return_data is set to `True`, the function returns a dataframe of the data used for the plotting.

**Example**

```
# Heatmap of foldchnage
sm.pl.foldchange (adata, label='foldchange', method='heatmap',
                  p_val=0.05, nonsig_color='grey',
                  cmap = 'vlag', log=True, center=0, linecolor='black',linewidths=0.7,
                  vmin=-5, vmax=5, row_cluster=False)
    
# Parallel_coordinates plot of the foldchanges
foldchange (adata, label='foldchange', 
            log=True, method='parallel_coordinates', invert_axis=True,
            parallel_coordinates_color=['black','blue','green','red','#000000'],
            matplotlib_bbox_to_anchor=(1.04,1),
            matplotlib_legend_loc='upper left',
            xticks_rotation=90,
            return_data = False
```
