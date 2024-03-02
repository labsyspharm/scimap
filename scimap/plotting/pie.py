#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Tue Apr  6 22:05:14 2021
# @author: Ajit Johnson Nirmal
"""
!!! abstract "Short Description"
    `sm.pl.pie`: This function facilitates the creation of pie charts to visually 
    represent the proportions of categories within any selected categorical column in 
    an AnnData object. It provides an intuitive and straightforward way to assess 
    the distribution of cell types, clusters, or any other categorical annotations, 
    offering insights into the composition of the dataset.

## Function
"""

# Lib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

# Function
def pie (adata, 
         phenotype='phenotype', 
         group_by='imageid', 
         ncols=None,
         subset_phenotype=None, 
         subset_groupby=None,
         label='auto', 
         title='auto', 
         colors=None, 
         autopct='%1.1f%%',
         legend=False,
         legend_loc='upper right',
         wedgeprops = {'linewidth': 0}, 
         return_data=False, **kwargs):
    """
Parameters:
        adata (anndata.AnnData):  
            The annotated data matrix.

        phenotype (str, optional):  
            Column in `adata.obs` containing the categorical data for pie chart visualization. 

        group_by (str, optional):  
            Column in `adata.obs` for defining groups. Each group will have its own pie chart. 
            Default is 'imageid'. Pass None to treat all data as a single group.

        ncols (int, optional):  
            Number of columns in the grid layout when displaying multiple pie charts. 
            Only applicable if `group_by` is used.

        subset_phenotype (list, optional):  
            List of categories within `phenotype` to include in the visualization.

        subset_groupby (list, optional):  
            List of groups within `group_by` to include in the visualization.

        label (list, optional):  
            Labels for each wedge in the pie charts. If 'auto', labels are automatically 
            derived from `phenotype` categories.

        title (str, optional):  
            Title for the pie chart(s). If 'auto', titles are derived from `group_by` categories.

        colors (list, optional):  
            Custom color sequence for the pie chart wedges.

        autopct (str or callable, optional):  
            String or function used to label wedges with their numeric value. Default is '%1.1f%%'.

        legend (bool, optional):  
            Whether to display a legend for the pie chart. Default is False.

        legend_loc (str, optional):  
            Location of the legend. Default is 'upper right'.

        wedgeprops (dict, optional):  
            Properties passed to the wedge objects, such as `{'linewidth': 3}`.

        return_data (bool, optional):  
            If True, returns the data used for plotting instead of the pie chart(s).

        **kwargs:  
            Additional keyword arguments passed to `matplotlib.pyplot.pie`.

Returns:
    plot and dataFrame (matplotlib, pandas DF):
        If `return_data` is True, returns a pandas DataFrame used for plotting.

Example:
    ```python
    
    # Basic pie chart visualization of cell phenotypes
    sm.pl.pie(adata, phenotype='cell_type', group_by='sample_id', ncols=3)

    # Advanced visualization with custom colors and pie chart properties
    sm.pl.pie(adata, phenotype='cell_type', group_by='condition', ncols=4, colors=['#ff9999','#66b3ff','#99ff99'],
              wedgeprops={'edgecolor': 'black', 'linewidth': 2}, autopct='%1.1f%%', legend=True, legend_loc='best')

    # Subsetted visualization focusing on specific phenotypes and groups
    sm.pl.pie(adata, phenotype='cell_type', group_by='treatment', subset_phenotype=['T cells', 'B cells'],
              subset_groupby=['Control', 'Treated'], ncols=2, legend=True, legend_loc='lower left')
    
    ```
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
        ax.pie(prop.value, labels=label,colors=colors, wedgeprops = wedgeprops)
        #ax.pie(prop.value, labels=label,colors=colors, wedgeprops = wedgeprops, **kwargs)
        if title is None:
            pass
        else:
            ax.set_title(phenotype)
    else:
        # plot the figure
        # Ground work for removing unwanted axes
        total_axes = list(range(ncols * rows))
        required_axes = list(range(len(np.unique(prop['group_by']))))
        final_axes = list(set(total_axes) ^ set(required_axes))      
        # Plot
        fig, axes = plt.subplots(ncols=ncols, nrows=rows)
        for (c, grp), ax in zip(g, axes.flat):
            ax.pie(grp.value, labels=label, colors=colors, wedgeprops =wedgeprops)
            #ax.pie(grp.value, labels=label, colors=colors, wedgeprops = wedgeprops, **kwargs)
            if title is None:
                pass
            else:
                ax.set_title(c)        
        # removing unwanted axis
        for i in final_axes:
            fig.delaxes(axes.flatten()[i])
        
        if legend is True:
            plt.legend(labels, loc=legend_loc, framealpha=1)
            
    plt.show()
    

    # return data
    if return_data is True:
        return prop
    
    