**scimap.pl.pie**

!!! note "Function Call"
    `scimap.pl.pie` (
      **adata, 
      phenotype='phenotype', 
      group_by='imageid', 
      ncols=None,
      subset_phenotype=None, 
      subset_groupby=None,
      label='auto', 
      title='auto', 
      colors=None, 
      autopct='%1.1f%%',
      wedgeprops = {'linewidth': 0}, 
      return_data=False 
      **kwargs **)

**Short description**

The `sm.pl.pie` function allows users to visualize proportions of cell-types or any
categorical data in the form of a pie plot.


**Parameters**

`adata` : AnnData object   

`phenotype` : string, optional *(The default is `phenotype`)*  
        Column contaning the cell-type inforamtion or any categorical data to be displayed
        in the form of a pie plot.
        
`group_by` : string, optional *(The default is `imageid`)*  
        Column that contains inforamtion on data groupings which leads to generation of
        pie plot for each group (e.g. image-id). If `None` is passed,
        the entire data is considered as a single group.
        
`ncols` : int, optional *(The default is `None`)*  
        In case group_by is used, a grid of plots are returned. This paramenter
        controls the number of columns in that grid. 
        
`subset_phenotype` : list, optional *(The default is `None`)*  
        User can subset a list of categories within `phenotype` before plotting.
        
`subset_groupby` : list, optional *(The default is `None`)*  
        User can subset a list of categories within `group_by` before plotting.
        
`label` : list, optional *(The default is `auto`)*  
        A list of strings providing the labels for each wedge.
        
`title` : string, optional *(The default is `auto`)*  
        If `None`, the title of the pieplot is not plotted.
        
`colors` : list, optional *(The default is `None`)*  
        A sequence of colors through which the pie chart will cycle. If None, will use the 
        colors in the currently active cycle. 
        
`autopct` : None or str or callable, optional *(The default is `'%1.1f%%'`)*  
        If not None, is a string or function used to label the wedges with their numeric value. 
        The label will be placed inside the wedge. If it is a format string, 
        the label will be fmt % pct. If it is a function, it will be called.

`wedgeprops` : dict, optional *(The default is `{'linewidth': 0}`)*  
        Dict of arguments passed to the wedge objects making the pie. For example, you can pass in 
        wedgeprops = {'linewidth': 3} to set the width of the wedge border lines equal to 3. 
        For more details, look at the doc/arguments of the wedge object.

`return_data` : bool, optional *(The default is False)*  
        When True, return the data used for plotting.
        
`**kwargs` : Keyword arguments to pass on to `matplotlib.pyplot.pie`


**Returns**

If return_data is set to `True`, the function returns a dataframe of the data used for the plotting.

**Example**

```
# pie plot showing stromal tumor content among the different samples
sm.pl.pie (adata, phenotype='Tumor_Stroma', group_by='imageid', 
           autopct='%1.1f%%',
           textprops={'fontsize': 8, 'color': '#1d3557', 'fontweight': 'bold'},
           ncols=5, label=None, title=None, 
           colors=['#a8dadc','#e63946'], 
           wedgeprops = {'linewidth': 0.8})
```
