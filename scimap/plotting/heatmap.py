# -*- coding: utf-8 -*-
#Created on Thu Nov 10 09:53:28 2022
#@author: Ajit Johnson Nirmal
#Function to incorporate the feature matrix and probability matrix into an adata object

"""
!!! abstract "Short Description"
    The `sm.pl.heatmap` function generates a comprehensive visualization of marker expression or other relevant features across various groups or clusters identified within spatial datasets. Through customizable clustering, normalization, and annotation features, it supports detailed exploratory data analysis and comparison across different conditions or phenotypes. This function effectively consolidates complex datasets into intuitive visual representations, enhancing the interpretability of high-dimensional data.

## Function
"""


# import packages
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.preprocessing import StandardScaler
import os
from matplotlib.colors import Normalize
import pandas as pd
import anndata as ad
import argparse


plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi']=300
plt.rcParams['font.family']='sans serif'
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['pdf.fonttype']=42


# function
def heatmap (adata,
             groupBy,
             layer=None,
             subsetMarkers=None,
             subsetGroups=None,
             clusterRows=True,
             clusterColumns=True,
             standardScale=None,
             orderRow=None, 
             orderColumn=None, 
             showPrevalence=False,
             cmap='vlag',
             figsize=None,
             saveDir=None, 
             fileName='heatmap.pdf',
             verbose=True,
             **kwargs
             ):
    """

Parameters:
    adata (AnnData): 
        An AnnData object or `path` to an Anndata object containing the dataset to be visualized. It should have features as variables and observations as rows. 
        
    groupBy (str): 
        The key in `adata.obs` on which to group observations. Typically, this will be a clustering pr phenotype label like 'leiden' or 'phenotype'.
        
    layer (str, optional): 
        Specifies the layer of `adata` to use for the heatmap. If None, the `.X` attribute is used. If you want to plot the raw data use `raw`
        
    subsetMarkers (list of str, optional): 
        A list of marker genes or features to include in the heatmap. If None, all markers are used.
        
    subsetGroups (list of str, optional): 
        A list of group labels to include in the heatmap. Useful for focusing on specific clusters or conditions.
        
    clusterRows (bool): 
        Whether to cluster rows (observations). 
        
    clusterColumns (bool): 
        Whether to cluster columns (features). 
        
    standardScale (str, optional): 
        Determines if and how to normalize the data across rows or columns. Acceptable values are 'row', 'column', or None.
        
    orderRow (list of str, optional): 
        Specifies a custom order for the rows based on group labels. 
        
    orderColumn (list of str, optional): 
        Specifies a custom order for the columns based on feature names.
        
    showPrevalence (bool): 
        If True, adds a bar showing the prevalence of the feature across the groups. 
        
    cmap (str): 
        The colormap for the heatmap. 
        
    figsize (tuple of float, optional): 
        The size of the figure to create. If None, the size is inferred.
        
    saveDir (str, optional): 
        Directory to save the generated heatmap. If None, the heatmap is not saved. 
        
    fileName (str, optional): 
        Name of the file to save the heatmap. Relevant only if `saveDir` is not None. 
        
    verbose (bool): 
        If True, print additional information during execution.
        
    **kwargs: 
        Additional keyword arguments are passed to the underlying matplotlib plotting function.

Returns:
    plot (matplotlib): 
        Returns a plot, if `saveDir` and `fileName` are provided, the plot is saved in the given directory.

Example:
        ```python
    
        # Example 1: Basic usage with clustering and standard scale by column.
    
        sm.pl.heatmap(adata, groupBy='leiden', standardScale='column')
    
        # Example 2: Advanced usage with specified subset markers, custom grouping, and file saving.
    
        subsetMarkers = ['ELANE', 'CD57', 'CD45', 'CD11B', 'SMA', 'CD16', 'ECAD']
        subsetGroups = ['0', '1', '3', '6']
        orderRow = ['6', '3', '0', '1']
        orderColumn = ['SMA', 'CD16', 'ECAD', 'ELANE', 'CD57', 'CD45', 'CD11B']
        saveDir = '/path/to/save'
        fileName = 'custom_heatmap.pdf'
        
        sm.pl.heatmap(adata, groupBy='leiden', subsetMarkers=subsetMarkers, subsetGroups=subsetGroups, clusterRows=False, clusterColumns=False, standardScale='column', orderRow=orderRow, orderColumn=orderColumn, showPrevalence=True, figsize=(10, 5), saveDir=saveDir, fileName=fileName, vmin=0, vmax=1)
        
        ```
"""
    
    # load adata
    if isinstance(adata, str):
        adata = ad.read_h5ad(adata)
    
    
    # check if the location is provided if the user wishes to save the image
    if (saveDir is None and fileName is not None) or (saveDir is not None and fileName is None):
        raise ValueError("Both 'saveDir' and 'fileName' must be provided together or not at all.")


    # subset data if user requests
    subsetadata = None # intialize subsetted data
    if subsetGroups: 
        subsetGroups = [subsetGroups] if isinstance(subsetGroups, str) else subsetGroups # convert to list 
        subsetadata = adata[adata.obs[groupBy].isin(subsetGroups)]
        # also identify the categories to be plotted
        categories = subsetadata.obs[groupBy].values
    else:
        # also identify the categories to be plotted
        categories = adata.obs[groupBy].values
    
    # subset the markers if user requests
    if subsetMarkers:
        subsetMarkers = [subsetMarkers] if isinstance(subsetMarkers, str) else subsetMarkers # convert to list 
        if subsetadata:
            # isolate the data
            if layer == 'raw':
                data = subsetadata[:, subsetMarkers].raw.X
            elif layer is None:
                data = subsetadata[:, subsetMarkers].X
            else:
                data = subsetadata[:, subsetMarkers].layers[layer]
        else:
            # isolate the data
            if layer == 'raw':
                data = adata[:, subsetMarkers].raw.X
            elif layer is None:
                data = adata[:, subsetMarkers].X
            else:
                data = adata[:, subsetMarkers].layers[layer]
    else:
        # take the whole data if the user does not subset anything
        if layer == 'raw':
            data = adata.raw.X
        elif layer is None:
            data = adata.X
        else:
            data = adata.layers[layer]

    
    # intialize the markers to be plotted
    if subsetMarkers is None:
        subsetMarkers = adata.var.index.tolist()
    
    # The actual plotting function
    def plot_category_heatmap_vectorized(data,
                                         marker_names, 
                                         categories, 
                                         clusterRows, 
                                         clusterColumns, 
                                         standardScale, 
                                         orderRow, 
                                         orderColumn, 
                                         showPrevalence,
                                         cmap,
                                         figsize,
                                         saveDir, 
                                         fileName,
                                         **kwargs):
        # Validate clustering and ordering options
        if (clusterRows or clusterColumns) and (orderRow is not None or orderColumn is not None):
            raise ValueError("Cannot use clustering and manual ordering together. Please choose one or the other.")
        
        if standardScale not in [None, 'row', 'column']:
            raise ValueError("standardScale must be 'row', 'column', or None.")
    
        # Convert marker_names to list if it's a pandas Index
        #if isinstance(marker_names, pd.Index):
        #    marker_names = marker_names.tolist()
        
        # Data preprocessing
        sorted_indices = np.argsort(categories)
        data = data[sorted_indices, :]
        categories = categories[sorted_indices]
        unique_categories, category_counts = np.unique(categories, return_counts=True)
        
        # Compute mean values for each category
        mean_data = np.array([np.mean(data[categories == category, :], axis=0) for category in unique_categories])
        
        # Apply standard scaling if specified
        if standardScale == 'row':
            scaler = StandardScaler()
            mean_data = scaler.fit_transform(mean_data)
        elif standardScale == 'column':
            scaler = StandardScaler()
            mean_data = scaler.fit_transform(mean_data.T).T
        
        # Apply manual ordering if specified
        if orderRow:
            # Ensure orderRow is a list
            if isinstance(orderRow, pd.Index):
                orderRow = orderRow.tolist()
            row_order = [unique_categories.tolist().index(r) for r in orderRow]
            mean_data = mean_data[row_order, :]
            unique_categories = [unique_categories[i] for i in row_order]
            category_counts = [category_counts[i] for i in row_order]
    
        if orderColumn:
            # Ensure orderColumn is a list
            if isinstance(orderColumn, pd.Index):
                orderColumn = orderColumn.tolist()
            col_order = [marker_names.index(c) for c in orderColumn]
            mean_data = mean_data[:, col_order]
            marker_names = [marker_names[i] for i in col_order]
        
            # Clustering
        if clusterRows:
            # Perform hierarchical clustering
            row_linkage = linkage(pdist(mean_data), method='average')
            # Reorder data according to the clustering
            row_order = dendrogram(row_linkage, no_plot=True)['leaves']
            mean_data = mean_data[row_order, :]
            unique_categories = unique_categories[row_order]
            category_counts = category_counts[row_order]
        
        if clusterColumns:
            # Perform hierarchical clustering
            col_linkage = linkage(pdist(mean_data.T), method='average')
            # Reorder data according to the clustering
            col_order = dendrogram(col_linkage, no_plot=True)['leaves']
            mean_data = mean_data[:, col_order]
            marker_names = [marker_names[i] for i in col_order]
    
        # Plotting
        # Dynamic figsize calculation
        if figsize is None:
            base_size = 0.5  # Base size for each cell in inches
            figsize_width = max(10, len(marker_names) * base_size)
            figsize_height = max(8, len(unique_categories) * base_size)
            figsize=(figsize_width, figsize_height)
        
        fig, ax = plt.subplots(figsize=figsize)
        
        
        # Heatmap
        # Extract vmin and vmax from kwargs if present, else default to min and max of mean_data
        vmin = kwargs.pop('vmin', np.min(mean_data))
        vmax = kwargs.pop('vmax', np.max(mean_data))
    
        # Create the Normalize instance with vmin and vmax
        norm = Normalize(vmin=vmin, vmax=vmax)
    
        c = ax.imshow(mean_data, aspect='auto', cmap=cmap, norm=norm, **kwargs)
        
        # Prevalence text
        if showPrevalence:
            # Calculate text offset from the last column of the heatmap
            text_offset = mean_data.shape[1] * 0.001  # Small offset from the right edge of the heatmap
            
            for index, count in enumerate(category_counts):
                # Position text immediately to the right of the heatmap
                ax.text(mean_data.shape[1] + text_offset, index, f"n={count}", va='center', ha='left')
        
    
        # Setting the tick labels
        ax.set_xticks(np.arange(mean_data.shape[1]))
        ax.set_xticklabels(marker_names, rotation=90, ha="right")
        ax.set_yticks(np.arange(mean_data.shape[0]))
        ax.set_yticklabels(unique_categories)
        
        # Move the colorbar to the top left corner
        cbar_ax = fig.add_axes([0.125, 0.92, 0.2, 0.02]) # x, y, width, height
        cbar = plt.colorbar(c, cax=cbar_ax, orientation='horizontal')
        cbar_ax.xaxis.set_ticks_position('top')
        cbar_ax.xaxis.set_label_position('top')
        cbar.set_label('Mean expression in group')
        
        ax.set_xlabel('Markers')
        ax.set_ylabel('Categories')
        
        plt.tight_layout(rect=[0, 0, 0.9, 0.9]) # Adjust the layout
        
        # Saving the figure if saveDir and fileName are provided
        if saveDir:
            if not os.path.exists(saveDir):
                os.makedirs(saveDir)
            full_path = os.path.join(saveDir, fileName)
            plt.savefig(full_path, dpi=300)
            plt.close(fig)
            print(f"Saved heatmap to {full_path}")
        else:
            plt.show()
        
        
    # call the plotting function
    plot_category_heatmap_vectorized(data=data,
                                         marker_names=subsetMarkers, 
                                         categories=categories, 
                                         clusterRows=clusterRows, 
                                         clusterColumns=clusterColumns, 
                                         standardScale=standardScale, 
                                         orderRow=orderRow, 
                                         orderColumn=orderColumn, 
                                         showPrevalence=showPrevalence, 
                                         cmap=cmap,
                                         figsize=figsize,
                                         saveDir=saveDir, 
                                         fileName=fileName,
                                         **kwargs)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate a heatmap visualization of features across different groups.')
    
    parser.add_argument('--adata', type=str, required=True, help='Path to an AnnData object file or AnnData object containing the dataset to be visualized.')
    parser.add_argument('--groupBy', type=str, required=True, help="The key in `adata.obs` on which to group observations.")
    parser.add_argument('--layer', type=str, default=None, help="Specifies the layer of `adata` to use for the heatmap. Use 'raw' to plot raw data.")
    parser.add_argument('--subsetMarkers', type=str, nargs='+', default=None, help="A list of marker genes or features to include in the heatmap.")
    parser.add_argument('--subsetGroups', type=str, nargs='+', default=None, help="A list of group labels to include in the heatmap.")
    parser.add_argument('--clusterRows', type=bool, default=True, help="Whether to cluster rows (observations).")
    parser.add_argument('--clusterColumns', type=bool, default=True, help="Whether to cluster columns (features).")
    parser.add_argument('--standardScale', type=str, choices=['row', 'column', None], default=None, help="Normalize the data across rows or columns.")
    parser.add_argument('--orderRow', type=str, nargs='+', default=None, help="Custom order for the rows based on group labels.")
    parser.add_argument('--orderColumn', type=str, nargs='+', default=None, help="Custom order for the columns based on feature names.")
    parser.add_argument('--showPrevalence', type=bool, default=False, help="Adds a bar showing the prevalence of the feature across the groups.")
    parser.add_argument('--cmap', type=str, default='vlag', help="The colormap for the heatmap.")
    parser.add_argument('--figsize', type=float, nargs=2, default=None, help="The size of the figure to create.")
    parser.add_argument('--saveDir', type=str, default=None, help="Directory to save the generated heatmap.")
    parser.add_argument('--fileName', type=str, default=None, help="Name of the file to save the heatmap.")
    parser.add_argument('--verbose', type=bool, default=True, help="If True, print additional information during execution.")
    
    args = parser.parse_args()
    
    # Load AnnData if path is provided instead of AnnData object
    if args.adata.endswith('.h5ad'):  # Simple check to assume a path is provided
        adata = ad.read_h5ad(args.adata)
    else:
        raise ValueError("Please provide a valid path to an AnnData object file.")
    
    # Convert bool strings to bool
    args.clusterRows = args.clusterRows in ['True', 'true', True]
    args.clusterColumns = args.clusterColumns in ['True', 'true', True]
    args.showPrevalence = args.showPrevalence in ['True', 'true', True]
    args.verbose = args.verbose in ['True', 'true', True]

    # Call heatmap function
    heatmap(adata,
            args.groupBy,
            layer=args.layer,
            subsetMarkers=args.subsetMarkers,
            subsetGroups=args.subsetGroups,
            clusterRows=args.clusterRows,
            clusterColumns=args.clusterColumns,
            standardScale=args.standardScale,
            orderRow=args.orderRow,
            orderColumn=args.orderColumn,
            showPrevalence=args.showPrevalence,
            cmap=args.cmap,
            figsize=args.figsize,
            saveDir=args.saveDir,
            fileName=args.fileName,
            verbose=args.verbose)