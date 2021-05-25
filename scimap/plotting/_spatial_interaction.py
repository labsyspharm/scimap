#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Wed Oct 21 11:00:57 2020
# @author: Ajit Johnson Nirmal
"""
!!! abstract "Short Description"
    `sm.pl.spatial_interaction`: The function allows users to plot a heatmap to visualize spatial interaction output. 
    The intensity represents abundance of co-occurrence (scaled) observed and blank regions represent non-significant results.

## Function
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
                         subset_phenotype=None, subset_neighbour_phenotype=None,
                         binary_view=False, return_data=False, **kwargs):
    """
Parameters:
    adata : AnnData object

    spatial_interaction : string, optional  
        In order to locate the spatial_interaction data within the AnnData object please provide the output 
        label/columnname of `sm.tl.spatial_interaction` function.

    summarize_plot : bool, optional  
        In the event of analyzing multiple images, this argument allows users to
        plot the average cell-cell interaction across all images.

    p_val : float, optional  
        P-value cut-off above which interactions are not considered significant.

    row_cluster : bool, optional  
        Cluster Rows.

    col_cluster : bool, optional  
        Cluster Columns.

    subset_phenotype : list, optional  
        If user requires to visualize a subset of phenotypes, it can be passed here. 
        e.g.  `subset_phenotype = ['celltype_A', 'celltype_B']`.

    subset_neighbour_phenotype : list, optional  
        If user requires to visualize a subset of interacting phenotypes, it can be passed here. 
        e.g.  `subset_neighbour_phenotype = ['celltype_C', 'celltype_D']`.

    cmap : string, optional  
        Color map to use for continous variables. 
        Can be a name or a Colormap instance (e.g. 'magma', 'viridis').

    nonsig_color : string, optional  
        Color for non-significant interactions (Interactions above the P-value cut-off will use this color).

    binary_view : bool, optional  
        Removes the intensity of intreaction and plots significant interactions and avoidance in a binary format.

    return_data : bool, optional  
        When True, return the data used for plotting.

    **kwargs : key:value pairs  
        Pass other parameters that works with `sns.clustermap`. e.g. `linecolor='black'`

Example:
```python
    # spatial_interaction heatmap for a single image
    sm.pl.spatial_interaction(adata, summarize_plot=True, 
    row_cluster=True, linewidths=0.75, linecolor='black')
    
    # spatial_interaction heatmap for multiple images
    sns.set(font_scale=0.6)
    sm.pl.spatial_interaction(adata, summarize_plot=False, 
    row_cluster=True, col_cluster=True, yticklabels=True)
```
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
        
    # subset the data if user requests
    if subset_phenotype is not None:
        if isinstance(subset_phenotype, str):
            subset_phenotype = [subset_phenotype]
        # subset the phenotype
        interaction_map = interaction_map[interaction_map['phenotype'].isin(subset_phenotype)]
    
    if subset_neighbour_phenotype is not None:
        if isinstance(subset_neighbour_phenotype, str):
            subset_neighbour_phenotype = [subset_neighbour_phenotype]
        # subset the phenotype
        interaction_map = interaction_map[interaction_map['neighbour_phenotype'].isin(subset_neighbour_phenotype)]
        
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
        
        # change to the order passed in subset
        if subset_phenotype is not None:
            interaction_map = interaction_map.reindex(subset_phenotype)
            p_val_df = p_val_df.reindex(subset_phenotype)
        if subset_neighbour_phenotype is not None:
            interaction_map = interaction_map.reindex(columns=subset_neighbour_phenotype)
            p_val_df = p_val_df.reindex(columns=subset_neighbour_phenotype)
        
        # Plotting heatmap
        mask = p_val_df.isnull() # identify the NAN's for masking 
        im = interaction_map.fillna(0) # replace nan's with 0 so that clustering will work
        # heatmap
        sns.clustermap(im, cmap=cmap, row_cluster=row_cluster, col_cluster=col_cluster,  mask=mask, **kwargs)
        
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
        
        # order the plot as needed
        if subset_phenotype or subset_neighbour_phenotype is not None:
            interaction_map.reset_index(inplace=True)
            p_val_df.reset_index(inplace=True)
            if subset_phenotype is not None:
                interaction_map['phenotype'] = interaction_map['phenotype'].astype('category')
                interaction_map['phenotype'] = interaction_map['phenotype'].cat.reorder_categories(subset_phenotype)
                interaction_map = interaction_map.sort_values('phenotype')
                # Do same for Pval
                p_val_df['phenotype'] = p_val_df['phenotype'].astype('category')
                p_val_df['phenotype'] = p_val_df['phenotype'].cat.reorder_categories(subset_phenotype)
                p_val_df = p_val_df.sort_values('phenotype')
            if subset_neighbour_phenotype is not None:
                interaction_map['neighbour_phenotype'] = interaction_map['neighbour_phenotype'].astype('category')
                interaction_map['neighbour_phenotype'] = interaction_map['neighbour_phenotype'].cat.reorder_categories(subset_neighbour_phenotype)
                interaction_map = interaction_map.sort_values('neighbour_phenotype')
                # Do same for Pval
                p_val_df['neighbour_phenotype'] = p_val_df['neighbour_phenotype'].astype('category')
                p_val_df['neighbour_phenotype'] = p_val_df['neighbour_phenotype'].cat.reorder_categories(subset_neighbour_phenotype)
                p_val_df = p_val_df.sort_values('neighbour_phenotype')
            if subset_phenotype and subset_neighbour_phenotype is not None:
                 interaction_map = interaction_map.sort_values(['phenotype', 'neighbour_phenotype'])
                 p_val_df = p_val_df.sort_values(['phenotype', 'neighbour_phenotype'])
            
            # convert the data back into multi-index
            interaction_map = interaction_map.set_index(['phenotype', 'neighbour_phenotype'])
            p_val_df = p_val_df.set_index(['phenotype', 'neighbour_phenotype'])
                 
        # Plotting heatmap
        mask = p_val_df.isnull() # identify the NAN's for masking 
        im = interaction_map.fillna(0) # replace nan's with 0 so that clustering will work
        mask.columns = im.columns
                
        # covert the first two columns into index
        # Plot
        sns.clustermap(im, cmap=cmap, row_cluster=row_cluster, col_cluster=col_cluster, mask=mask, **kwargs)
    
    if return_data is True:
        # perpare data for export
        map_data = interaction_map.copy()
        p_val_data = mask.copy()
        map_data.reset_index(inplace=True)
        p_val_data.reset_index(inplace=True)
        # remove the first two colums
        map_data = map_data.drop(['phenotype','neighbour_phenotype'],axis=1)
        p_val_data = p_val_data.drop(['phenotype','neighbour_phenotype'],axis=1)
        p_val_data.columns = map_data.columns
        # remove the mased values
        final_Data = map_data.where(~p_val_data, other=np.nan)
        final_Data.index = interaction_map.index
        return final_Data
        

