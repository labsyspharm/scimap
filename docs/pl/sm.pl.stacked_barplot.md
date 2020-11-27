**scimap.pl.stacked_barplot**

!!! note "Function Call"
    `scimap.pl.stacked_barplot` (
      **adata,
      x_axis='imageid',y_axis='phenotype',
      subset_xaxis=None,subset_yaxis=None,
      method='percent',
      plot_tool='matplotlib',matplotlib_cmap=None,
      matplotlib_bbox_to_anchor=(1,1.02), matplotlib_legend_loc=2, 
      return_data=False,
      **kwargs **)

**Short description**

The `sm.pl.stacked_barplot` function allows users to generate a stacked barplot of any two columns
within the `anndata.obs` object. It offers two methods- i.e either to plot the absolute number of cells 
or plot the number of cells as a proportion of total cells (normalized to 100%). 
The function also incorporates two plotting tools- `matplotlib` and `plotly`. When `plotly` is used,
the bar plot is opened in a local browser. This is especially useful when there are a large number of
categories and is difficult to match the colors with the keys. `plotly` allows users to hover over the
data and identify the group it belongs to. 

**Parameters**

`adata` : AnnData object   

`x_axis` : string, required *(The default is 'imageid')*  
        Column name of the data that need to be plotted in the x-axis.
        
`y_axis` : string, required *(The default is 'phenotype')*  
        Column name of the data that need to be plotted in the y-axis.
        
`subset_xaxis` : list, optional *(The default is None)*  
        Subset x-axis before plotting. Pass in a list of categories. `eg- subset_xaxis = ['ROI_1', 'ROI_5']`
        
`subset_yaxis` : list, optional *(The default is None)*  
        Subset y-axis before plotting. Pass in a list of categories. `eg- subset_yaxis = ['Celltype_A', 'Celltype_B']`
        
`method` : string, optional *(The default is 'percent')*  
        Available options: 'percent' and 'absolute'.  
        1) Use Percent to plot the percent proportion.  
        2) Use 'absolute' to plot the plot the absolute number.  
        
`plot_tool` : string, optional *(The default is 'matplotlib')*  
        Available options: 'matplotlib' and 'plotly'.  
        1) matplotlib uses the standard python plotting method  
        2) plotly opens the plot in a local browser. Advantage is to be able 
        to hover over the plot and retreive data for plots with large number of categories.
        
`matplotlib_cmap` : string, optional *(The default is None)*  
        Colormap to select colors from. If string, load colormap with that name from matplotlib. 
        
`matplotlib_bbox_to_anchor` : tuple, optional *(The default is (1,1.02))*  
        Bounding box argument used along with matplotlib_legend_loc to control
        the legend location when using the matplotlib method.
        
`matplotlib_legend_loc` : int, optional *(The default is 2)*  
        Location of legend used along with matplotlib_bbox_to_anchor to control
        the legend location when using the matplotlib method.
        
`return_data` : bool, optional *(The default is False)*  
        When True, return the data used for plotting.
        
`**kwargs` : Additional keyword arguments passed to:  
        1) Pandas [DataFrame.plot()](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.plot.html#pandas.DataFrame.plot) when using the `matplotlib` method.  
        2) [Plotly.bar()](https://plotly.com/python-api-reference/generated/plotly.express.bar.html) when using the `plotly` method.


[Duck Duck Go](https://duckduckgo.com)

**Returns**

Stacked bar plot. If return_data is set to `True` also returns a dataframe of the data used for the plot.

**Example**

```
# Plot the absolute number of phenotypes using the matplotlib tool across differnt ROI's
# ROI column is `epidermis_roi` and phenotypes are stored under `phenotype`

sm.pl.stacked_barplot (adata,x_axis='epidermis_roi',y_axis='phenotype',
                     method='absolute',plot_tool='matplotlib',figsize=(10, 10))
    
# Plot the number of cells normalized to 100% using the plotly tool across differnt ROI's
    
sm.pl.stacked_barplot (adata,x_axis='epidermis_roi',y_axis='phenotype',
                     method='percent',plot_tool='plotly',color_discrete_sequence=px.colors.qualitative.Alphabet)
    
# Same as above but visualize only a subset of ROI's and a subset of phenotypes
subset_xaxis = ['epidermis_1', 'epidermis_5', 'epidermis_6']
subset_yaxis = ['APCs', 'Keratinocytes', 'Macrophages']
    
sm.pl.stacked_barplot (adata,x_axis='epidermis_roi',y_axis='phenotype',
                     subset_xaxis=subset_xaxis,subset_yaxis=subset_yaxis,
                     method='percent',plot_tool='plotly')
    
# Visualize absolute number of phenotypes and return the data into a dataframe `absolute_number`
absolute_number = sm.pl.stacked_barplot (adata,x_axis='epidermis_roi',y_axis='phenotype',
                     method='absolute',plot_tool='matplotlib', return_data=True)
```
