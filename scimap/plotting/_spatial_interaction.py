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
import numpy as np
import pandas as pd
import matplotlib
sns.set_style("white")

# Function
def spatial_interaction (adata, spatial_interaction='spatial_interaction',
                         summarize_plot=True, p_val=0.05,
                         row_cluster=False, col_cluster=False,
                         cmap = 'vlag', nonsig_color='grey', 
                         binary_view=False, **kwargs):
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
    p_val : float, optional
        P-value cut-off above which interactions are not considered significant. The default is 0.05.
    row_cluster : bool, optional
        Cluster Rows. The default is False.
    col_cluster : bool, optional
        Cluster Columns. The default is False.
    cmap : string, optional
        Color map to use for continous variables. 
        Can be a name or a Colormap instance (e.g. 'magma', 'viridis'). The default is 'vlag'.
    nonsig_color : string, optional
        Color for non-significant interactions (Interactions above the P-value cut-off will use this color).
        The default is 'grey'.
    binary_view : bool, optional
        Removes the intensity of intreaction and plots significant interactions and avoidance in a binary format.
        The default is 'False'.

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
    
    # set color for heatmap
    #cmap_updated = copy.copy(matplotlib.cm.get_cmap(cmap))
    cmap_updated = matplotlib.cm.get_cmap(cmap)
    cmap_updated.set_bad(color=nonsig_color)
    

    # Copy the interaction results from anndata object
    try:
        interaction_map = adata.uns[spatial_interaction].copy()
    except KeyError:
        raise ValueError('spatial_interaction not found- Please run sm.tl.spatial_interaction first')
        
    # Seperate Interaction intensity from P-value
    p_value = interaction_map.filter(regex='pvalue_')    
    p_val_df = pd.concat([interaction_map[['phenotype','neighbour_phenotype']], p_value], axis=1, join='outer')
    p_val_df = p_val_df.set_index(['phenotype','neighbour_phenotype'])
    interaction_map = interaction_map[interaction_map.columns.difference(p_value.columns)]
    interaction_map = interaction_map.set_index(['phenotype','neighbour_phenotype'])
    
    # Binarize the values if user requests
    if binary_view == True:
        interaction_map[interaction_map > 0] = 1
        interaction_map[interaction_map <= 0] = -1
        

    if summarize_plot == True:
        # convert first two columns to multi-index column
        #interaction_map = interaction_map.set_index(['phenotype','neighbour_phenotype'])
        #p_val_df = p_val_df.set_index(['phenotype','neighbour_phenotype'])
        
        # If multiple images are present, take the average of interactions
        interaction_map['mean'] = interaction_map.mean(axis=1).values
        interaction_map = interaction_map[['mean']] # keep only the mean column
        interaction_map = interaction_map['mean'].unstack()
        # Do the same for P-values
        p_val_df['mean'] = p_val_df.mean(axis=1).values
        p_val_df = p_val_df[['mean']] # keep only the mean column
        # set the P-value threshold
        p_val_df.loc[p_val_df[p_val_df['mean'] > p_val].index,'mean'] = np.NaN       
        p_val_df = p_val_df['mean'].unstack()
        
        # Plotting heatmap
        mask = p_val_df.isnull() # identify the NAN's for masking 
        interaction_map = interaction_map.fillna(0) # replace nan's with 0 so that clustering will work
        # heatmap
        sns.clustermap(interaction_map, cmap=cmap, row_cluster=row_cluster, col_cluster=col_cluster,  mask=mask, **kwargs)
        
    else:
        if len(interaction_map.columns) <= 3:
            raise ValueError('Data for only a single image is available please set summarize_plot=True and try again')
        # convert first two columns to multi-index column
        #interaction_map = interaction_map.set_index(['phenotype','neighbour_phenotype'])
        #p_val_df = p_val_df.set_index(['phenotype','neighbour_phenotype'])
        
        # P value threshold
        p_val_df = p_val_df.apply(lambda x: np.where(x > p_val,np.nan,x))
        
        # Remove rows that are all nan
        idx = p_val_df.index[p_val_df.isnull().all(1)] # Find all nan rows
        interaction_map = interaction_map.loc[interaction_map.index.difference(idx)] # clean intensity data
        p_val_df = p_val_df.loc[p_val_df.index.difference(idx)] # clean p-value data
        
        # Plotting heatmap
        mask = p_val_df.isnull() # identify the NAN's for masking 
        interaction_map = interaction_map.fillna(0) # replace nan's with 0 so that clustering will work
        mask.columns = interaction_map.columns
        # Plot
        sns.clustermap(interaction_map, cmap=cmap, row_cluster=row_cluster, col_cluster=col_cluster, mask=mask, **kwargs)
        

