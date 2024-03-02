# -*- coding: utf-8 -*-
#Created on Wed May  3 16:09:24 2023
#@author: Ajit Johnson Nirmal
#Plot a 2D density plot of the expression of two markers.

"""
!!! abstract "Short Description"
    The `sm.pl.densityPlot2D` function crafts 2D density plots to visualize expression 
    levels of one or two specified markers. When a single marker is provided, it depicts 
    its expression distribution across the dataset. With two markers, it contrasts the expression 
    of the first against the second, offering insights into their co-expression or distribution patterns. 

## Function
"""

import pandas as pd
import math
import matplotlib.pyplot as plt
import mpl_scatter_density
import pathlib
import matplotlib as mpl
import copy
from matplotlib import cm
from matplotlib.colors import LogNorm
mpl.rcParams['pdf.fonttype'] = 42


def densityPlot2D (adata, 
                   markerA,  markerB=None, 
                   layer=None, 
                   subset=None, 
                   imageid='imageid', 
                   ncols=None, 
                   cmap='jet', 
                   figsize=(3, 3), 
                   hline = 'auto', vline = 'auto',
                   fontsize=None, 
                   dpi=100, 
                   xticks=None,
                   yticks=None,
                   outputDir=None, 
                   outputFileName='densityPlot2D.pdf'):

    """
Parameters:
    adata (anndata.AnnData): 
        Annotated data matrix containing single-cell gene expression data.

    markerA (str):  
        The name of the first marker whose expression will be plotted.

    markerB (list, optional):  
        The name of the second marker or a list of second markers whose expression will be plotted. 
        If not provided, a 2D density plot of `markerA` against all markers in the dataset will be plotted.

    layer (str or list of str, optional):  
        The layer in adata.layers that contains the expression data to use. 
        If None, adata.X is used. use `raw` to use the data stored in `adata.raw.X`

    subset (list, optional):  
        `imageid` of a single or multiple images to be subsetted for plotting purposes.

    imageid (str, optional):  
        Column name of the column containing the image id. Use in conjunction with `subset`.

    ncols (int, optional):  
        The number of columns in the grid of density plots.

    cmap (str, optional):  
        The name of the colormap to use. Defaults to 'jet'.

    figsize (tuple, optional):  
        The size of the figure in inches.

    hline (float or 'auto', optional):  
        The y-coordinate of the horizontal line to plot. If set to `None`, a horizontal line is not plotted. 
        Use 'auto' to draw a vline at the center point. 

    vline (float or 'auto', optional):  
        The x-coordinate of the vertical line to plot. If set to `None`, a vertical line is not plotted. 
        Use 'auto' to draw a vline at the center point. 

    fontsize (int, optional):  
        The size of the font of the axis labels.

    dpi (int, optional):  
        The DPI of the figure. Use this to control the point size. Lower the dpi, larger the point size.

    xticks (list of float, optional):  
        Custom x-axis tick values.

    yticks (list of float, optional):  
        Custom y-axis tick values.

    outputDir (str, optional):  
        The directory to save the output plot.

    outputFileName (str, optional):  
        The name of the output file. Use desired file format as suffix (e.g. `.png` or `.pdf`).

Returns:
    Plot (image):  
        If `outputDir` is not provided, the plot is displayed on the screen. 
        Otherwise, the plot is saved in the provided `outputDir` directory.

Example:
    ```python
    
    # create a 2D density plot of the expression of 'CD3D' against 'CD8A' in the dataset 'adata'
    sm.pl.densityPlot2D(adata, markerA='CD3D', markerB='CD8A')

    # create a 2D density plot of the expression of 'CD3D' against all markers in the dataset 'adata'
    sm.pl.densityPlot2D(adata, markerA='CD3D')
    ```

    """
    # testing
    # import anndata as ad
    # adata = ad.read(r"C:\Users\aj\Dropbox (Partners HealthCare)\nirmal lab\softwares\scimap\scimap\tests\_data\example_data.h5ad")
    # adata = ad.read('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/softwares/scimap/scimap/tests/_data/example_data.h5ad')
    # markerA ='CD3E'; layers=None; markerB='CD163'; plotGrid=True; ncols=None; color=None; figsize=(10, 10); fontsize=None; subset=None; imageid='imageid'; xticks=None; dpi=200; outputDir=None; 
    # hline = 'auto'; vline = 'auto'
    # outputFileName='densityPlot2D.png'
    # color = {'markerA': '#000000', 'markerB': '#FF0000'}
    # outputDir = r"C:\Users\aj\Downloads"
    
    #densityPlot2D (adata, markerA='CD3D', markerB=['CD2', 'CD10', 'CD163'], dpi=50, outputDir=r"C:\Users\aj\Downloads")
    
    
    # set color
    cp = copy.copy(cm.get_cmap(cmap))
    cp.set_under(alpha=0)
    
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
    x = data[markerA]
    
    if markerB is None:
        y = data.drop(markerA, axis=1)
    else:
        if isinstance(markerB, str):
            markerB = [markerB]
        y = data[markerB]
        
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
    num_rows, num_cols = calculate_grid_dimensions(len(y.columns), num_columns = ncols)
    
    
    
    fig, axs = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(num_cols*figsize[0],num_rows*figsize[0]), subplot_kw={'projection': 'scatter_density'})
    if num_rows == 1 and num_cols == 1:
        axs = [axs]  # wrap single subplot in a list
    else:
        axs = axs.flatten()
    for i, col in enumerate(y.columns):
        ax = axs[i]
        ax.scatter_density(x, y[col], dpi=dpi, cmap=cp, norm=LogNorm(vmin=0.5, vmax=x.size))
        ax.set_xlabel(markerA, size = fontsize)
        ax.set_ylabel(col, size = fontsize)
        
        if hline == 'auto':
            ax.axhline((y[col].max() + y[col].min()) / 2, color='grey')
        elif hline is None:
            pass
        else:
            ax.axhline(hline, color='grey')

        if vline == 'auto':
            ax.axvline((x.max() + x.min()) / 2, color='grey')
        elif vline is None:
            pass
        else:
            ax.axvline(vline, color='grey')
        
        # control and x and y ticks
        if xticks is not None:
            ax.set_xticks(xticks)
            ax.set_xticklabels([str(x) for x in xticks])
        
        if yticks is not None:
            ax.set_yticks(yticks)
            ax.set_yticklabels([str(x) for x in yticks])
        
        
    # Remove any empty subplots
    num_plots = len(y.columns)
    for i in range(num_plots, num_rows * num_cols):
        ax = axs[i]
        fig.delaxes(ax)
    
    plt.tick_params(axis='both', labelsize=fontsize)
    plt.tight_layout()

    # Save the figure to a file
    # save figure
    if outputDir is not None:
        plt.savefig(pathlib.Path(outputDir) / outputFileName)
    



        
        
        
        