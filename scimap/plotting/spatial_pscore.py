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
import seaborn as sns

sns.set_theme(style='white', color_codes=True)
import matplotlib.pyplot as plt
import os

# plt.rcParams['figure.dpi'] = 100
# plt.rcParams['savefig.dpi']=300
# plt.rcParams['font.family']='sans serif'
# plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['pdf.fonttype'] = 42


# Function
def spatial_pscore(
    adata,
    label='spatial_pscore',
    plot_score='both',
    order_xaxis=None,
    color='grey',
    figsize=None,
    fileName='spatial_pscore.pdf',
    saveDir=None,
    **kwargs,
):
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

            fileName (str, optional):
                Name of the file to save the plot. Relevant only if `saveDir` is not None.

            saveDir (str, optional):
                Directory to save the generated plot. If None, the plot is not saved.

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
    x = data.index
    y_pd = data['Proximity Density'].values
    y_pv = data['Proximity Volume'].values

    if figsize is None:
        # Dynamically calculate figsize based on the data size
        figsize_width_scale = (
            1.0  # Adjust based on your preference and the expected data length
        )
        figsize_height_scale = 0.5  # Adjust this to change how tall the plots are
        # For 'both', we might want a wider figure
        figsize_width = max(12, len(x) * figsize_width_scale)
        figsize_height = max(8, len(x) * figsize_height_scale)
        figsize = (figsize_width, figsize_height)

    # Plot what user requests
    if plot_score == 'Proximity Density':
        fig, ax = plt.subplots(figsize=figsize)
        sns.barplot(x=x, y=y_pd, color=color, **kwargs).set_title('Proximity Density')
        plt.xticks(rotation=90)
        plt.tight_layout()
    elif plot_score == 'Proximity Volume':
        fig, ax = plt.subplots(figsize=figsize)
        sns.barplot(x=x, y=y_pv, color=color, **kwargs).set_title('Proximity Volume')
        plt.xticks(rotation=90)
        plt.tight_layout()
    elif plot_score == 'both':
        fig, axs = plt.subplots(1, 2, figsize=figsize)
        sns.barplot(x=x, y=y_pd, color=color, ax=axs[0], **kwargs).set_title(
            'Proximity Density'
        )
        axs[0].tick_params(axis='x', rotation=90)
        sns.barplot(x=x, y=y_pv, color=color, ax=axs[1], **kwargs).set_title(
            'Proximity Volume'
        )
        axs[1].tick_params(axis='x', rotation=90)
        plt.tight_layout()

    # Saving the figure if saveDir and fileName are provided
    if saveDir:
        if not os.path.exists(saveDir):
            os.makedirs(saveDir)
        full_path = os.path.join(saveDir, fileName)
        fig.savefig(full_path, dpi=300)
        plt.close(fig)
        print(f"Saved plot to {full_path}")
    else:
        plt.show()
