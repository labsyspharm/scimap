# -*- coding: utf-8 -*-
#Created on Fri Mar 17 10:00:18 2023
#@author: Ajit Johnson Nirmal
#Plotting function to visualize postive cells

"""
!!! abstract "Short Description"
    The scatterPlot function offers a convenient way to generate scatter plots 
    for visualizing single-cell spatial data. By utilizing this function, 
    users can effectively visualize the spatial distribution of cells while 
    overlaying expression levels or categorical columns onto the plot. 
    This functionality allows for a comprehensive understanding of the 
    relationship between cell location and specific features of interest 
    within the dataset.

## Function
"""

# Libs
import anndata as ad
import pathlib
import matplotlib.pyplot as plt
import argparse
import pandas as pd
import math
import numpy as np
import matplotlib.patches as mpatches
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42


# Function
def spatial_scatterPlot (adata, 
                         colorBy, 
                         x_coordinate='X_centroid',
                         y_coordinate='Y_centroid',
                         imageid='imageid',
                         layer=None,
                         subset=None,
                         s=None,
                         ncols=None,
                         alpha=1,
                         dpi=200,
                         fontsize=None,
                         plotLegend=True,
                         cmap='RdBu_r',
                         catCmap='tab20',
                         vmin=None,
                         vmax=None,
                         customColors=None,
                         figsize=(5, 5),
                         invert_yaxis=True,
                         outputDir=None,
                         outputFileName='scimapScatterPlot.png',
                         **kwargs):
    """
Parameters:
    adata (anndata.AnnData):   
        Pass the `adata` loaded into memory or a path to the `adata` 
        file (.h5ad).
    
    colorBy (str):  
    The column name that will be used for color-coding the points. This can be 
    either markers (data stored in `adata.var`) or observations (data stored in `adata.obs`).

    x_coordinate (str, optional):  
        The column name in `spatial feature table` that records the
        X coordinates for each cell.
    
    y_coordinate (str, optional):  
        The column name in `single-cell spatial table` that records the
        Y coordinates for each cell. 
        
    imageid (str, optional):  
        The column name in `spatial feature table` that contains the image ID 
        for each cell. 
        
    layer (str or None, optional):  
        The layer in `adata.layers` that contains the expression data to use. 
        If `None`, `adata.X` is used. use `raw` to use the data stored in `adata.raw.X`.
        
    subset (list or None, optional):  
        `imageid` of a single or multiple images to be subsetted for plotting purposes.
        
    s (float, optional):  
        The size of the markers. 
        
    ncols (int, optional):  
        The number of columns in the final plot when multiple variables are plotted.
        
    alpha (float, optional):   
        The alpha value of the points (controls opacity).
        
    dpi (int, optional):   
        The DPI of the figure.
    
    fontsize (int, optional):    
        The size of the fonts in plot.
        
    plotLegend (bool, optional):   
        Whether to include a legend. 
        
    cmap (str, optional):   
        The colormap to use for continuous data.
        
    catCmap (str, optional):   
        The colormap to use for categorical data.
        
    vmin (float or None, optional):   
        The minimum value of the color scale. 
        
    vmax (float or None, optional):   
        The maximum value of the color scale. 
        
    customColors (dict or None, optional):   
        A dictionary mapping color categories to colors. 
        
    figsize (tuple, optional):   
        The size of the figure. Default is (5, 5).

    invert_yaxis (bool, optional):  
        Invert the Y-axis of the plot. 
        
    outputDir (str or None, optional):   
        The directory to save the output plot. If None, the plot will not be saved. 
        
    outputFileName (str, optional):   
        The name of the output file. Use desired file format as 
        suffix (e.g. `.png` or `.pdf`). Default is 'scimapScatterPlot.png'.
        
    **kwargs:  
        Additional keyword arguments to be passed to the matplotlib scatter function.


Returns:
    Plot (image):
        If `outputDir` is provided the plot will saved within the
        provided outputDir.

Example:
        ```python
        
        customColors = { 'Unknown' : '#e5e5e5',
                        'CD8+ T' : '#ffd166',
                        'Non T CD4+ cells' : '#06d6a0',
                        'CD4+ T' : '#118ab2',
                        'ECAD+' : '#ef476f',
                        'Immune' : '#073b4c',
                        'KI67+ ECAD+' : '#000000'               
            }
        
        sm.pl.spatial_scatterPlot (adata=core6, 
                         colorBy = ['ECAD', 'phenotype_gator'], 
                         subset = 'unmicst-6_cellMask',
                         figsize=(4,4),
                         s=0.5,
                         plotLegend=True,
                         fontsize=3,
                         dpi=300,
                         vmin=0,
                         vmax=1,
                         customColors=customColors,
                         outputFileName='scimapScatterPlot.svg',
                         outputDir='/Users/aj/Downloads')


        ```
        
    """
    
    # Load the andata object
    if isinstance(adata, str):
        adata = ad.read(adata)
    else:
        adata = adata.copy()

    # subset data if neede
    if subset is not None:
        if isinstance (subset, str):
            subset = [subset]
        if layer == 'raw':
            bdata=adata.copy()
            bdata.X = adata.raw.X
            bdata = bdata[bdata.obs[imageid].isin(subset)]
        else:
            bdata=adata.copy()
            bdata = bdata[bdata.obs[imageid].isin(subset)]
    else:
        bdata=adata.copy()

    # isolate the data
    if layer is None:
        data = pd.DataFrame(bdata.X, index=bdata.obs.index, columns=bdata.var.index)
    elif layer == 'raw':
        data = pd.DataFrame(bdata.raw.X, index=bdata.obs.index, columns=bdata.var.index)
    else:
        data = pd.DataFrame(bdata.layers[layer], index=bdata.obs.index, columns=bdata.var.index)
    
    # isolate the meta data
    meta = bdata.obs
    
    # identify the things to color
    if isinstance (colorBy, str):
        colorBy = [colorBy]   
    # extract columns from data and meta
    data_cols = [col for col in data.columns if col in colorBy]
    meta_cols = [col for col in meta.columns if col in colorBy]
    # combine extracted columns from data and meta
    colorColumns = pd.concat([data[data_cols], meta[meta_cols]], axis=1)
        
    # identify the x and y coordinates
    x = meta[x_coordinate]
    y = meta[y_coordinate]

    
    # auto identify rows and columns in the grid plot
    def calculate_grid_dimensions(num_items, num_columns=None):
        """
        Calculates the number of rows and columns for a square grid
        based on the number of items.
        """
        if num_columns is None:
            num_rows_columns = int(math.ceil(math.sqrt(num_items)))
            return num_rows_columns, num_rows_columns
        else:
            num_rows = int(math.ceil(num_items / num_columns))
            return num_rows, num_columns

    # calculate the number of rows and columns
    nrows, ncols = calculate_grid_dimensions(len(colorColumns.columns), num_columns = ncols)
    
    
    # resolve figsize
    #figsize = (figsize[0]*ncols, figsize[1]*nrows)
    
    # Estimate point size
    if s is None:
        s = (10000 / bdata.shape[0]) / len(colorColumns.columns)
    
    # Define the categorical colormap (optional)
    cmap_cat = plt.get_cmap(catCmap)
        
    # FIIGURE
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize, dpi=dpi)
    
    # Flatten the axs array for easier indexing
    if nrows == 1 and ncols == 1:
        axs = [axs]  # wrap single subplot in a list
    else:
        axs = axs.flatten()

    
    # Loop over the columns of the DataFrame
    for i, col in enumerate(colorColumns):
        # Select the current axis
        ax = axs[i]
        
        # invert y-axis
        if invert_yaxis is True:
            ax.invert_yaxis()
        
        # Scatter plot for continuous data
        if colorColumns[col].dtype.kind in 'iufc':
            scatter = ax.scatter(x=x, y=y, 
                                 c=colorColumns[col], 
                                 cmap=cmap, 
                                 s=s,
                                 vmin=vmin,
                                 vmax=vmax,
                                 linewidths=0,
                                 alpha=alpha, **kwargs)
            if plotLegend is True:
                cbar = plt.colorbar(scatter, ax=ax, pad=0)
                cbar.ax.tick_params(labelsize=fontsize)
        
        # Scatter plot for categorical data
        else:
            # Get the unique categories in the column
            categories = colorColumns[col].unique()
            
            # Map the categories to colors using either the custom colors or the categorical colormap
            if customColors:
                colors = [customColors.get(cat, cmap_cat(i)) for i, cat in enumerate(categories)]
            else:
                colors = [cmap_cat(i) for i in np.linspace(0, 1, len(categories))]
            
            # Map the categories to numeric codes for plotting
            codes = [np.where(categories == cat)[0][0] for cat in colorColumns[col]]

            
            # Plot the scatter plot with categorical colors
            c = [colors[code] for code in codes]
            scatter = ax.scatter(x=x, y=y, c=c, s=s, linewidths=0, alpha=alpha, **kwargs)
            if plotLegend is True:
                # Create the categorical legend outside the plot
                handles = [mpatches.Patch(color=colors[i], label=cat) for i, cat in enumerate(categories)]
                ax.legend(handles=handles, bbox_to_anchor=(1.0, 1.0), loc='upper left', bbox_transform=ax.transAxes, fontsize=fontsize)
        
        ax.set_title(col) # fontsize=fontsize
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
    
    # Remove any empty subplots
    num_plots = len(colorColumns.columns)
    for i in range(num_plots, nrows * ncols):
        ax = axs[i]
        fig.delaxes(ax)
                
    # Adjust the layout of the subplots grid
    plt.tick_params(axis='both', labelsize=fontsize)
    plt.tight_layout()

    # save figure
    if outputDir is not None:
        plt.savefig(pathlib.Path(outputDir) / outputFileName, dpi=dpi)
    
    #plt.show()
    



