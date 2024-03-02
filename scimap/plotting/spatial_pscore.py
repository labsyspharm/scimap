#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Sat Jan 30 21:38:21 2021
# @author: Ajit Johnson Nirmal

"""
!!! abstract "Short Description"
    `sm.pl.spatial_pscore`: This function offers a visual representation of 
    proximity volume and density scores, essential for understanding the spatial 
    relationships and interactions among cell types or phenotypes within tissue samples. 
    To ensure accurate and meaningful visualizations, it is crucial to compute 
    these scores beforehand using `sm.tl.spatial_pscore`. Through customizable bar 
    plots, users can delve into the intricacies of spatial co-occurrence patterns, 
    facilitating deeper insights into the cellular microenvironment.

## Function
"""

# Library
import seaborn as sns; sns.set_theme(style='white', color_codes=True)
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

# Function
def spatial_pscore (adata, 
                    label='spatial_pscore', 
                    plot_score='both', 
                    order_xaxis = None,
                    color='grey', **kwargs):
    """
Parameters:
        adata (anndata.AnnData):  
            The annotated data matrix with spatial proximity scores.

        label (str, optional):  
            The label/key used to access the spatial proximity scores stored in `adata.uns`. 
            This should match the label used during the `sm.tl.spatial_pscore` computation. Default is 'spatial_pscore'.

        plot_score (str, optional):  
            Determines which score(s) to plot. Options are:
            - 'Proximity Density' for plotting only the Proximity Density scores,
            - 'Proximity Volume' for plotting only the Proximity Volume scores,
            - 'both' for plotting both scores side by side. Default is 'both'.

        order_xaxis (list, optional):  
            Custom order for the x-axis categories. Pass a list of category names in the desired order.
            This can be useful for comparing specific regions or samples in a specific sequence.

        color (str, optional):  
            Color to use for the bar plots. This can enhance plot readability or align with publication themes.

        **kwargs:  
            Additional keyword arguments passed directly to seaborn's barplot function, allowing for further customization of the plots.

Returns:
        Plot (matplotlib): 
            Displays the generated bar plots.

Example:
    ```python
    
    # Basic visualization of both Proximity Density and Proximity Volume
    sm.pl.spatial_pscore(adata, label='spatial_pscore', plot_score='both', color='skyblue')

    # Customized plot for Proximity Density with ordered x-axis and specific color
    sm.pl.spatial_pscore(adata, plot_score='Proximity Density', order_xaxis=['Sample1', 'Sample2', 'Sample3'], color='salmon', edgecolor: 'black'})

    # Focused plot on Proximity Volume with seaborn customization through kwargs
    sm.pl.spatial_pscore(adata, plot_score='Proximity Volume', color='lightgreen', saturation: 0.8, alpha: 0.7})
    
    ```
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
        ax = sns.barplot(x=x, y=y_pd, color=color, **kwargs).set_title('Proximity Density')
        ax = plt.xticks(rotation=90)
        plt.tight_layout()
    if plot_score == 'Proximity Volume':
        ax = sns.barplot(x=x, y=y_pv, color=color, **kwargs).set_title('Proximity Volume')
        ax = plt.xticks(rotation=90)
        plt.tight_layout()
    if plot_score == 'both':
        fig, ax = plt.subplots(1,2)
        sns.barplot(x=x, y=y_pd, color=color, ax=ax[0], **kwargs).set_title('Proximity Density')
        ax[0].tick_params(axis='x', rotation=90)
        sns.barplot(x=x, y=y_pv, color=color, ax=ax[1], **kwargs).set_title('Proximity Volume')
        ax[1].tick_params(axis='x', rotation=90)
        plt.tight_layout()
        fig.show()