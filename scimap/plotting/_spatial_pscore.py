#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 21:38:21 2021
@author: Ajit Johnson Nirmal
Function to plot the proximity scores
"""


# Library
import seaborn as sns; sns.set_theme(style='white', color_codes=True)
import matplotlib.pyplot as plt

# Function
def spatial_pscore (adata, label='spatial_pscore', plot_score='both', 
                    order_xaxis = None,
                    color='grey', **kwargs):
    """
    

    Parameters
    ----------
    adata : AnnData object
    
    label : string, optional
        The label under which the data is saved. This is the same `label` parameter 
        passed when running the `sm.tl.spatial_pscore` function.
        The default is 'spatial_pscore'.
    plot_score : string, optional
        Three option are available. 
        A) Plot only the *Proximity Density* by passing in `Proximity Density`
        B) Plot only the *Proximity Volume* by passing in `Proximity Volume`
        C) Plot both side by side by passing `both`
        The default is 'both'.
    order_xaxis : list, optional
        If the user wants to re-order the x-axis, pass all the names in the x-axis
        in the desired order as a list. e.g. ['ROI2', 'ROI4', "ROI1"] 
        The default is None.
    color : string, optional
        Color of the bars. The default is 'grey'.
    **kwargs : string
        Other arguments that can be passed into `sns.barplot`

    Example
    -------
    sm.pl.spatial_pscore (adata, color='Black', plot_score='Proximity Volume', order_xaxis=['ROI2', 'ROI4', "ROI1"])

    """
    
    
    # Isolate the data from anndata object
    data = adata.uns[label]
    
    # Order the x-axis if needed
    if order_xaxis is not None:
        data = data.reindex(order_xaxis)
        
     
    # Generate the x and y axis
    x  = data.index
    y_pd = data['Proximity Density'].values
    y_pv = data['Proximity Volume'].values
    
    # Plot what user requests
    if plot_score == 'Proximity Density':
        ax = sns.barplot(x, y_pd, color=color, **kwargs).set_title('Proximity Density')
        ax = plt.xticks(rotation=90)
        plt.tight_layout()
    if plot_score == 'Proximity Volume':
        ax = sns.barplot(x, y_pv, color=color, **kwargs).set_title('Proximity Volume')
        ax = plt.xticks(rotation=90)
        plt.tight_layout()
    if plot_score == 'both':
        fig, ax = plt.subplots(1,2)
        sns.barplot(x, y_pd, color=color, ax=ax[0], **kwargs).set_title('Proximity Density')
        ax[0].tick_params(axis='x', rotation=90)
        sns.barplot(x, y_pv, color=color, ax=ax[1], **kwargs).set_title('Proximity Volume')
        ax[1].tick_params(axis='x', rotation=90)
        plt.tight_layout()
        fig.show()