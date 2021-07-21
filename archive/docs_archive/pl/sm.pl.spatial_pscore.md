**scimap.pl.spatial_pscore**

!!! note "Function Call"
    `scimap.pl.spatial_pscore` (
      **adata, 
      label='spatial_pscore', 
      plot_score='both', 
      order_xaxis = None,
      color='grey', 
      **kwargs **)

**Short description**

The `sm.pl.spatial_pscore` function allows users to generate a barplot of the spatial pscores. Using
the function the users can plot either the `Proximity Density` score, `Proximity Volume` score or both togeather.

**Parameters**

`adata` : AnnData object   

`label` : string, required *(The default is 'spatial_pscore')*  
        The label under which the data is saved. This is the same `label` parameter 
        passed when running the `sm.tl.spatial_pscore` function.  
        
`plot_score` : string, required *(The default is 'both')*  
        Three option are available.  
        A) Plot only the *Proximity Density* by passing in `Proximity Density`  
        B) Plot only the *Proximity Volume* by passing in `Proximity Volume`  
        C) Plot both side by side by passing `both`  
        
`order_xaxis` : list, optional *(The default is None)*  
        If the user wants to re-order the x-axis, pass all the names in the 
        x-axis in the desired order as a list. e.g. ['ROI2', 'ROI4', "ROI1"] 
        
`color` : list, optional *(The default is 'grey')*  
        Color of the bars.
        
`**kwargs` : Additional keyword arguments that can be passed to `sns.barplot`  


**Example**

```
sm.pl.spatial_pscore (adata, color='Black', plot_score='Proximity Volume', order_xaxis=['ROI2', 'ROI4', "ROI1"])
```
