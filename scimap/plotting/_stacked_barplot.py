#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Fri Nov 27 09:34:22 2020
# @author: Ajit Johnson Nirmal

"""
!!! abstract "Short Description"
    `sm.pl.stacked_barplot`: The function allows users to generate a stacked bar plot of a categorical column. 
    The function can generate the plots using matplotlib and Plotly libraries. Plotly is browser based and so 
    it can be used for interactive data exploration.

## Function
"""

# Required library
import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns; sns.set(color_codes=True)
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default = 'browser'
sns.set(style="white")


# Function
def stacked_barplot (adata, x_axis='imageid', y_axis='phenotype', subset_xaxis=None, subset_yaxis=None, 
                     order_xaxis=None, order_yaxis=None,
                     method='percent', plot_tool='matplotlib', matplotlib_cmap=None, 
                     matplotlib_bbox_to_anchor=(1,1.02), matplotlib_legend_loc=2, 
                     return_data=False, **kwargs):
    """
Parameters:
    adata : AnnData Object

    x_axis : string, required  
        Column name of the data that need to be plotted in the x-axis.

    y_axis : string, required  
        Column name of the data that need to be plotted in the y-axis.

    subset_xaxis : list, optional  
        Subset x-axis before plotting. Pass in a list of categories. eg- subset_xaxis = ['ROI_1', 'ROI_5']

    subset_yaxis : list, optional  
        Subset y-axis before plotting. Pass in a list of categories. eg- subset_yaxis = ['Celltype_A', 'Celltype_B']

    order_xaxis : list, optional  
        Order the x-axis of the plot as needed. Pass in a list of categories. eg- order_xaxis = ['ROI_5', 'ROI_1']
        The default is None and will be plotted based on alphabetic order. Please note that if you change the order, pass all categories, failure to do so
        will generate NaN's.

    order_yaxis : list, optional  
        Order the y-axis of the plot as needed. Pass in a list of categories. eg- order_yaxis = ['Celltype_B', 'Celltype_A']
        The default is None and will be plotted based on alphabetic order. Please note that if you change the order, pass all categories, failure to do so
        will generate NaN's.

    method : string, optional  
        Available options: 'percent' and 'absolute'. 
        1) Use Percent to plot the percent proportion.  
        2) Use 'absolute' to plot the plot the absolute number.  

    plot_tool : string, optional  
        Available options: 'matplotlib' and 'plotly'.  
        1) matplotlib uses the standard python plotting method  
        2) plotly opens the plot in a local browser. Advantage is to be able   
        to hover over the plot and retreive data for plots with large number of categories.

    matplotlib_cmap : string, optional  
        Colormap to select colors from. If string, load colormap with that name from matplotlib. 

    matplotlib_bbox_to_anchor : tuple, optional  
        Bounding box argument used along with matplotlib_legend_loc to control
        the legend location when using the matplotlib method.

    matplotlib_legend_loc : int, optional  
        Location of legend used along with matplotlib_bbox_to_anchor to control
        the legend location when using the matplotlib method.

    return_data : bool, optional  
        When True, return the data used for plotting.

    **kwargs : Additional keyword arguments passed to:  
        1) Pandas DataFrame.plot() when using the `matplotlib` method (https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.plot.html#pandas.DataFrame.plot))  
        2) Plotly.bar() when using the `plotly` method (https://plotly.com/python-api-reference/generated/plotly.express.bar.html))  

Returns:  
    Stacked bar plot. If return_data is set to `True` also returns a dataframe of the data used for the plot.

    
Example:
```python
    # Plot the absolute number of phenotypes using the matplotlib 
    # tool across differnt ROI's
    # ROI column is `epidermis_roi` and phenotypes are stored under `phenotype`

    sm.pl.stacked_barplot (adata,x_axis='epidermis_roi',y_axis='phenotype',
                     method='absolute',plot_tool='matplotlib',
                     figsize=(10, 10))
    
    # Plot the number of cells normalized to 100% using the plotly 
    # tool across differnt ROI's
    
    sm.pl.stacked_barplot (adata,x_axis='epidermis_roi',y_axis='phenotype',
                     method='percent',plot_tool='plotly',
                     color_discrete_sequence=px.colors.qualitative.Alphabet)
    
    # Same as above but visualize only a subset of ROI's and a subset of 
    # phenotypes
    subset_xaxis = ['epidermis_1', 'epidermis_5', 'epidermis_6']
    subset_yaxis = ['APCs', 'Keratinocytes', 'Macrophages']
    
    sm.pl.stacked_barplot (adata,x_axis='epidermis_roi',y_axis='phenotype',
                            subset_xaxis=subset_xaxis,subset_yaxis=subset_yaxis,
                            method='percent',plot_tool='plotly')
    
    # Visualize absolute number of phenotypes and return the data into a 
    # dataframe `absolute_number`
    absolute_number = sm.pl.stacked_barplot (adata,x_axis='epidermis_roi',
                      y_axis='phenotype', method='absolute',
                      plot_tool='matplotlib', return_data=True)

```
    """
    
    
    # create the dataframe with details
    data = pd.DataFrame(adata.obs)[[x_axis,y_axis]].astype(str)
    
    # subset the data if needed
    #if subset_data is not None:data = data[data[list(subset_data.keys())[0]].isin(list(subset_data.values())[0])]
    
    if subset_xaxis is not None:
        if isinstance(subset_xaxis, str):
            subset_xaxis = [subset_xaxis]
        data = data[data[x_axis].isin(subset_xaxis)]
    if subset_yaxis is not None:
        if isinstance(subset_yaxis, str):
            subset_yaxis = [subset_yaxis]
        data = data[data[y_axis].isin(subset_yaxis)]
        
        
    # Method: Absolute or Percentile
    if method == 'percent':
        total = data.groupby([x_axis,y_axis]).size().unstack().fillna(0).sum(axis=1)
        rg = pd.DataFrame(data.groupby([x_axis,y_axis]).size().unstack().fillna(0).div(total, axis=0).stack())
    elif method == 'absolute':
        rg = pd.DataFrame(data.groupby([x_axis,y_axis]).size().unstack().fillna(0).stack())
    else:
        raise ValueError('method should be either percent or absolute')
    
    # change column name
    rg.columns = ['count']
    
    # Add the index as columns in the data frame    
    rg.reset_index(inplace=True)  
    
    # re-order the x oy y axis if requested by user
    if order_xaxis is not None:
        rg[x_axis] = rg[x_axis].astype('category')
        rg[x_axis] = rg[x_axis].cat.reorder_categories(order_xaxis)
        rg = rg.sort_values(x_axis)
    if order_yaxis is not None:
        rg[y_axis] = rg[y_axis].astype('category')
        rg[y_axis] = rg[y_axis].cat.reorder_categories(order_yaxis)
        rg = rg.sort_values(y_axis)
    if order_xaxis and order_yaxis is not None:
        rg = rg.sort_values([x_axis, y_axis])
        
    pivot_df = rg.pivot(index=x_axis, columns=y_axis, values='count')
        
    # Plotting tool
    if plot_tool == 'matplotlib':
        
        if matplotlib_cmap is None:
            if len(rg[y_axis].unique()) <= 9:
                matplotlib_cmap = "Set1"        
            elif len(rg[y_axis].unique()) > 9 and len(rg[y_axis].unique()) <=20:
                matplotlib_cmap = plt.cm.tab20      #tab20  
            else:
                matplotlib_cmap = plt.cm.gist_ncar
        
        # Plotting
        # add width if not passed via parameters
        try:
            width
        except NameError:
            width=0.9
        # actual plotting   
        p = pivot_df.plot.bar(stacked=True, cmap=matplotlib_cmap, width=width, **kwargs)
        p.legend(bbox_to_anchor=matplotlib_bbox_to_anchor, loc=matplotlib_legend_loc)
    
    elif plot_tool == 'plotly':
        
        fig = px.bar(rg, x=x_axis, y="count", color=y_axis, **kwargs)
        fig.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)',
                           'paper_bgcolor': 'rgba(0, 0, 0, 0)'},
                          xaxis = dict(tickmode='linear') #type = 'category'
                          )
        fig.show()
        
        
    else:
        
        raise ValueError('plot_tool should be either matplotlib or plotly')
    
    # Return data
    if return_data is True:
        return pivot_df
        

    