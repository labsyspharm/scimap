# -*- coding: utf-8 -*-
#Created on Mon Apr 24 10:42:32 2023
#@author: Ajit Johnson Nirmal
#Create distribution plots
"""
!!! abstract "Short Description"
    The `sm.pl.distPlot` function is used to create distribution plots of 
    marker intensity data.

## Function
"""

import pandas as pd
import math
import matplotlib.pyplot as plt
import itertools
import pathlib
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def distPlot(adata, 
             layer=None, 
             markers=None, 
             subset=None, 
             imageid='imageid',
             vline=None,
             plotGrid=True, 
             ncols=None, 
             color=None, 
             xticks=None, 
             figsize=(5, 5), 
             fontsize=None, 
             dpi=200, 
             outputDir=None, 
             outputFileName='scimapDistPlot.png'):

    """
Parameters:
    adata (AnnData): 
        Annotated data object.
    
    layer (str, optional): 
        Layer of data to plot.
    
    markers (list, optional): 
        List of marker genes to plot.
    
    subset (list or None, optional):  
        `imageid` of a single or multiple images to be subsetted for plotting purposes.
    
    imageid (str, optional):  
        The column name in `spatial feature table` that contains the image ID 
        for each cell. 
    
    vline (float or 'auto', optional):  
        The x-coordinate of the vertical line to plot. If set to `None`, a vertical line is not plotted.
        Use 'auto' to draw a vline at the center point. 
    
    plotGrid (bool, optional):  
        Whether to plot each marker in it's own sub plot. If `False` and multiple markers 
        are passed in via `markers`, all distributions will be plotted within a single plot.
    
    ncols (int, optional):  
        The number of columns in the final plot when multiple variables are plotted.
    
    color (str, optional):   
        Color of the distribution plot. 
    
    xticks (list of float, optional):  
        Custom x-axis tick values.
    
    figsize (tuple, optional):   
        Figure size. Defaults to (5, 5).
    
    fontsize (int, optional):  
        The size of the font of the axis labels.
    
    dpi (int, optional):  
        The DPI of the figure. Use this to control the point size. Lower the dpi, larger the point size.
    
    outputDir (str, optional):  
        The directory to save the output plot.

    outputFileName (str, optional):  
        The name of the output file. Use desired file format as suffix (e.g. `.png` or `.pdf`).

Returns:
    Plot (image):
        If `outputDir` is provided the plot will saved within the
        provided outputDir.

Example:

        ```python
        
        sm.pl.distPlot(adata, 
                     layer=None, 
                     markers=['CD45','CD3D','CD20'], 
                     plotGrid=True, 
                     ncols=5)
    
    """
    
    # testing
    # layers=None; markers=None; plotGrid=True; ncols=None; color=None; figsize=(10, 10); fontsize=None; subset=None; imageid='imageid'; xticks=None; dpi=200; outputDir=None; 
    # outputFileName='distPlot.png'
    # color = {'markerA': '#000000', 'markerB': '#FF0000'}
    # outputDir = r"C:\Users\aj\Downloads"
    
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
    
    # keep only columns that are required
    if markers is not None:
        if isinstance(markers, str):
            markers = [markers]
        # subset the list
        data = data[markers]

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
        
    
    if plotGrid is False:
        # Create a figure and axis object
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)           
        # Loop through each column in the DataFrame and plot a KDE with the
        # user-defined color or the default color (grey)
        if color is None:
            for column in data.columns:
                data[column].plot.kde(ax=ax, label=column)
        else:
            for column in data.columns:
                c = color.get(column, 'grey')
                data[column].plot.kde(ax=ax, label=column, color=c)
        ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=fontsize)
        ax.tick_params(axis='both', which='major', width=1, labelsize=fontsize)
        plt.tight_layout()
        if xticks is not None:
            ax.set_xticks(xticks)
            ax.set_xticklabels([str(x) for x in xticks])       
            

        if vline == 'auto':
            ax.axvline((data[column].max() + data[column].min()) / 2, color='black')
        elif vline is None:
            pass
        else:
            ax.axvline(vline, color='black')
                
                
        # save figure
        if outputDir is not None:
            plt.savefig(pathlib.Path(outputDir) / outputFileName)
            
    else:
        # calculate the number of rows and columns
        num_rows, num_cols = calculate_grid_dimensions(len(data.columns), num_columns = ncols)
        
        # set colors
        if color is None:
            # Define a color cycle of 10 colors
            color_cycle = itertools.cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
            # Assign a different color to each column
            color = {col: next(color_cycle) for col in data.columns}
                
        # Set the size of the figure
        fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=figsize, dpi=dpi)
        # Set the spacing between subplots
        #fig.subplots_adjust(bottom=0.1, hspace=0.1)

        # Loop through each column in the DataFrame and plot a KDE with the
        # user-defined color or the default color (grey) in the corresponding subplot
        for i, column in enumerate(data.columns):
            c = color.get(column, 'grey')
            row_idx = i // num_cols
            col_idx = i % num_cols
            data[column].plot.kde(ax=axes[row_idx, col_idx], label=column, color=c)
            axes[row_idx, col_idx].set_title(column)
            axes[row_idx, col_idx].tick_params(axis='both', which='major', width=1, labelsize=fontsize)
            axes[row_idx, col_idx].set_ylabel('')
            
            if vline == 'auto':
                axes[row_idx, col_idx].axvline((data[column].max() + data[column].min()) / 2, color='black')
            elif vline is None:
                pass
            else:
                axes[row_idx, col_idx].axvline(vline, color='black')
            
            if xticks is not None:
                axes[row_idx, col_idx].set_xticks(xticks)
                axes[row_idx, col_idx].set_xticklabels([str(x) for x in xticks])
        
        # Remove any empty subplots
        num_plots = len(data.columns)
        for i in range(num_plots, num_rows * num_cols):
            row_idx = i // num_cols
            col_idx = i % num_cols
            fig.delaxes(axes[row_idx, col_idx])
        
        # Set font size for tick labels on both axes
        plt.tick_params(axis='both', labelsize=fontsize)
        plt.tight_layout()
    
        # Save the figure to a file
        # save figure
        if outputDir is not None:
            plt.savefig(pathlib.Path(outputDir) / outputFileName)
        

