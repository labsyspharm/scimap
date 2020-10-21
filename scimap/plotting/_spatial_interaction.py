#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 11:00:57 2020
@author: Ajit Johnson Nirmal
Heatmap plot to visualize spatial interaction output. The intensity represents 
number of interactions (scaled) observed and blank regions represent non-significant results.
"""

# Library
import seaborn as sns; sns.set(color_codes=True)
sns.set_style("white")

# Function
def spatial_interaction (adata, spatial_interaction='spatial_interaction',
                         summarize_plot=True, 
                         row_cluster=False, col_cluster=False,
                         cmap = 'vlag', **kwargs):
    """

    Parameters
    ----------
    adata : AnnData object

    spatial_interaction : string, optional
        In order to locate the spatial_interaction data within the AnnData object please provide the output 
        label/columnname of `sm.tl.spatial_interaction` function. The default is 'spatial_interaction'.
    summarize_plot : bool, optional
        In the event of analyzing multiple images, this argument allows users to
        plot the average cell-cell interaction across all images. The default is True.
    row_cluster : bool, optional
        Cluster Rows. The default is False.
    col_cluster : bool, optional
        Cluster Columns. The default is False.
    cmap : string, optional
        Color map to use for continous variables. 
        Can be a name or a Colormap instance (e.g. 'magma', 'viridis'). The default is 'vlag'.
    **kwargs : key:value pairs
        Are passed to sns.clustermap. Pass other parameters that works with `sns.clustermap`. e.g. `linecolor='black'`

    Example
    -------
    # spatial_interaction heatmap for a single image
    sm.pl.spatial_interaction(adata, summarize_plot=True, row_cluster=True, linewidths=0.75, linecolor='black')
    
    # spatial_interaction heatmap for multiple images
    sns.set(font_scale=0.6)
    sm.pl.spatial_interaction(adata, summarize_plot=False, row_cluster=True, col_cluster=True, yticklabels=True)

    """
    

    # Copy the interaction results from anndata object
    try:
        interaction_map = adata.uns[spatial_interaction].copy()
    except KeyError:
        raise ValueError('spatial_interaction not found- Please run sm.tl.spatial_interaction first')

    
    if summarize_plot == True:
        # convert first two columns to multi-index column
        interaction_map = interaction_map.set_index(['phenotype','neighbour_phenotype'])
        # If multiple images are present, take the average of interactions
        interaction_map['mean'] = interaction_map.mean(axis=1).values
        interaction_map = interaction_map[['mean']] # keep only the mean column
        interaction_map = interaction_map['mean'].unstack()
        # Plotting heatmap
        mask = interaction_map.isnull() # identify the NAN's for masking 
        interaction_map = interaction_map.fillna(0) # replace nan's with 0 so that clustering will work
        # heatmap
        sns.clustermap(interaction_map, cmap=cmap, row_cluster=row_cluster, col_cluster=col_cluster, mask=mask, **kwargs)
        
    else:
        if len(interaction_map.columns) <= 3:
            raise ValueError('Data for only a single image is available please set summarize_plot=True and try again')
        # convert first two columns to multi-index column
        interaction_map = interaction_map.set_index(['phenotype','neighbour_phenotype'])
        # Plotting heatmap
        mask = interaction_map.isnull() # identify the NAN's for masking 
        interaction_map = interaction_map.fillna(0) # replace nan's with 0 so that clustering will work
        # Plot
        sns.clustermap(interaction_map, cmap=cmap, row_cluster=row_cluster, col_cluster=col_cluster, mask=mask, **kwargs)