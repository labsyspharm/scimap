#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 22:05:14 2021
@author: Ajit Johnson Nirmal
Pie plot for categories of interest.
A pie plot is a proportional representation of the numerical data in a column. 
"""

# Lib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Function
def pie (adata, phenotype='phenotype', group_by='imageid', ncols=None,
         subset_phenotype=None, subset_groupby=None,
         label='auto', title='auto', colors=None, autopct='%1.1f%%',
         wedgeprops = {'linewidth': 0}, return_data=False, **kwargs):
    """
    

    Parameters
    ----------
    adata : AnnData object

    phenotype : string, optional
        Column contaning the cell-type inforamtion or any categorical data to be displayed
        in the form of a pie plot. The default is 'phenotype'.
    group_by : string, optional
        Column that contains inforamtion on data groupings which leads to generation of
        pie plot for each group (e.g. image-id). If `None` is passed,
        the entire data is considered as a single group. The default is 'imageid'.
    ncols : int, optional
        In case group_by is used, a grid of plots are returned. This paramenter
        controls the number of columns in that grid. The default is None.
    subset_phenotype : list, optional
        User can subset a list of categories within `phenotype` before plotting. The default is None.
    subset_groupby : TYPE, optional
        User can subset a list of categories within `group_by` before plotting. The default is None.
    label : list, optional
        A list of strings providing the labels for each wedge. The default is 'auto'.
    title : string, optional
        If `None`, the title of the pieplot is not plotted. The default is 'auto'.
    colors : list, optional
        A sequence of colors through which the pie chart will cycle. If None, will use the 
        colors in the currently active cycle. The default is None.
    autopct : None or str or callable, optional
        If not None, is a string or function used to label the wedges with their numeric value. 
        The label will be placed inside the wedge. If it is a format string, 
        the label will be fmt % pct. If it is a function, it will be called. The default is '%1.1f%%'.
    wedgeprops : dict, optional
        Dict of arguments passed to the wedge objects making the pie. For example, you can pass in 
        wedgeprops = {'linewidth': 3} to set the width of the wedge border lines equal to 3. 
        For more details, look at the doc/arguments of the wedge object. By default clip_on=False.
        The default is {'linewidth': 0}.
    return_data : bool, optional
        Returns the data used for plotting. The default is `False`.
    **kwargs : 
        Keyword arguments to pass on to `matplotlib.pyplot.pie`.

    Returns
    -------
    Returns data used for plotting if `return_data = True`
    
    Example
    -------
    # pie plot showing stromal tumor content among the different samples
    sm.pl.pie (adata, phenotype='Tumor_Stroma', group_by='imageid', 
               autopct='%1.1f%%',
               textprops={'fontsize': 8, 'color': '#1d3557', 'fontweight': 'bold'},
               ncols=5, label=None, title=None, 
               colors=['#a8dadc','#e63946'], 
               wedgeprops = {'linewidth': 0.8})

    """
    
    
    # convert subset into list
    if subset_phenotype is not None:
        if isinstance (subset_phenotype, str): 
            subset_phenotype = [subset_phenotype]
    if subset_groupby is not None:
        if isinstance (subset_groupby, str):
            subset_groupby = [subset_groupby]
    
    # create copy of the required data
    if group_by is not None:
        data = adata.obs[[phenotype,group_by]]
    else:
        data = adata.obs[[phenotype]]

    # subset data if needed
    if subset_groupby is not None:
        data = data[data[group_by].isin(subset_groupby)]
        data[group_by] = data[group_by].astype('str').astype('category')
        data[group_by] = data[group_by].cat.reorder_categories(subset_groupby)
        data = data.sort_values(group_by)
    if subset_phenotype is not None:
        data = data[data[phenotype].isin(subset_phenotype)]
        data[phenotype] = data[phenotype].astype('str').astype('category')      
        data[phenotype] = data[phenotype].cat.reorder_categories(subset_phenotype)
        data = data.sort_values(phenotype)
    if group_by and phenotype is not None:
        data = data.sort_values([phenotype, group_by])
        
        
    # calculate the proportion
    if group_by is None:
        prop = data[phenotype].value_counts().reset_index(inplace=False)
        prop.columns = [phenotype, 'value']
        prop['group_by'] = phenotype
        labels = np.unique(prop[phenotype])

    else:
        # if group_by is provided
        prop = pd.DataFrame(data.groupby([group_by,phenotype]).size()).reset_index(inplace=False)
        prop.columns = ['group_by',phenotype,'value']
        labels = np.unique(prop[phenotype])
        #
        if ncols is not None:
            g = prop.groupby('group_by')
            rows = int(np.ceil(len(g)/ncols))
        else:
            g = prop.groupby('group_by')
            rows = 1
            ncols = len(g)
    
    # remove label if requested 
    if label == 'auto':
        label = labels
    elif label is None:
        label = None
    else:
        label = label
    
    
    # plot
    if group_by is None:
        fig, ax = plt.subplots()
        #ax.pie(prop.value, labels=label,colors=colors, wedgeprops = wedgeprops)
        ax.pie(prop.value, labels=label,colors=colors, wedgeprops = wedgeprops, **kwargs)
        if title is None:
            pass
        else:
            ax.set_title(phenotype)
    else:
        # plot the figure
        fig, axes = plt.subplots(ncols=ncols, nrows=rows)
        for (c, grp), ax in zip(g, axes.flat):
            #ax.pie(grp.value, labels=label, colors=colors, wedgeprops =wedgeprops)
            ax.pie(grp.value, labels=label, colors=colors, wedgeprops = wedgeprops, **kwargs)
            if title is None:
                pass
            else:
                ax.set_title(c)               
    plt.show()
    

    # return data
    if return_data is True:
        return prop
    
    