#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Wed Oct 21 11:00:57 2020
# @author: Ajit Johnson Nirmal
"""
!!! abstract "Short Description"
    `sm.pl.spatial_interaction`: This function provides a sophisticated approach to 
    visualizing spatial interactions between cell types or phenotypes through heatmaps. 
    It adeptly highlights the frequency and significance of co-occurrence patterns, 
    where the intensity of colors reflects the scaled abundance of interactions. 
    Non-significant results are clearly marked with blank regions, offering a clear 
    demarcation of areas where interactions do not meet the specified threshold of significance. 
    This visualization tool is invaluable for uncovering intricate spatial relationships 
    and potential signaling networks within complex tissue environments.

## Function
"""

# Library
import seaborn as sns

sns.set(color_codes=True)
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os


sns.set_style("white")
plt.rcParams['pdf.fonttype'] = 42


# Function
def spatial_interaction(
    adata,
    spatial_interaction='spatial_interaction',
    summarize_plot=True,
    p_val=0.05,
    row_cluster=False,
    col_cluster=False,
    cmap='vlag',
    nonsig_color='grey',
    subset_phenotype=None,
    subset_neighbour_phenotype=None,
    binary_view=False,
    return_data=False,
    fileName='spatial_interaction.pdf',
    saveDir=None,
    **kwargs,
):
    """
    Parameters:
            adata (anndata.AnnData):
                The annotated data matrix with spatial interaction calculations.

            spatial_interaction (str, optional):
                Key in `adata.uns` where spatial interaction data is stored, typically the output of `sm.tl.spatial_interaction`.

            summarize_plot (bool, optional):
                If True, summarizes cell-cell interactions across all images or samples to provide an aggregated view.

            p_val (float, optional):
                Threshold for significance of interactions. Interactions with a P-value above this threshold are considered non-significant.

            row_cluster, col_cluster (bool, optional):
                If True, performs hierarchical clustering on rows or columns in the heatmap to group similar patterns of interaction.

            cmap (str, optional):
                Colormap for the heatmap visualization. Default is 'vlag'.

            nonsig_color (str, optional):
                Color used to represent non-significant interactions in the heatmap.

            subset_phenotype, subset_neighbour_phenotype (list, optional):
                Subsets of phenotypes or neighboring phenotypes to include in the analysis and visualization.

            binary_view (bool, optional):
                If True, visualizes interactions in a binary manner, highlighting presence or absence of significant interactions without intensity gradation.

            return_data (bool, optional):
                If True, returns the DataFrame used for plotting instead of the plot itself.

            fileName (str, optional):
                Name of the file to save the plot. Relevant only if `saveDir` is not None.

            saveDir (str, optional):
                Directory to save the generated plot. If None, the plot is not saved.

            **kwargs:
                Additional keyword arguments for seaborn's clustermap function, such as `linecolor` and `linewidths`.

    Returns:
            pandas.DataFrame (dataframe):
                Only if `return_data` is True. The DataFrame containing the data used for plotting.

    Example:
        ```python

        # Basic visualization of spatial interactions with default settings
        sm.pl.spatial_interaction(adata)

        # Detailed heatmap of spatial interactions, excluding non-significant interactions
        sm.pl.spatial_interaction(adata, summarize_plot=False, p_val=0.01, cmap='coolwarm', nonsig_color='lightgrey',
                            binary_view=True, row_cluster=True, col_cluster=True)

        # Visualizing specific phenotypes interactions, with custom colormap and binary view
        sm.pl.spatial_interaction(adata, subset_phenotype=['T cells', 'B cells'], subset_neighbour_phenotype=['Macrophages'],
                            cmap='seismic', binary_view=True, row_cluster=True, col_cluster=False,
                            figsize=(10, 8), dendrogram_ratio=(.1, .2), cbar_pos=(0, .2, .03, .4))
        ```
    """

    # set color for heatmap
    # cmap_updated = copy.copy(matplotlib.cm.get_cmap(cmap))
    # cmap_updated = matplotlib.cm.get_cmap(cmap)
    cmap_updated = matplotlib.colormaps[cmap]
    cmap_updated.set_bad(color=nonsig_color)

    # Copy the interaction results from anndata object
    try:
        interaction_map = adata.uns[spatial_interaction].copy()
    except KeyError:
        raise ValueError(
            'spatial_interaction not found- Please run sm.tl.spatial_interaction first'
        )

    # subset the data if user requests
    if subset_phenotype is not None:
        if isinstance(subset_phenotype, str):
            subset_phenotype = [subset_phenotype]
        # subset the phenotype
        interaction_map = interaction_map[
            interaction_map['phenotype'].isin(subset_phenotype)
        ]

    if subset_neighbour_phenotype is not None:
        if isinstance(subset_neighbour_phenotype, str):
            subset_neighbour_phenotype = [subset_neighbour_phenotype]
        # subset the phenotype
        interaction_map = interaction_map[
            interaction_map['neighbour_phenotype'].isin(subset_neighbour_phenotype)
        ]

    # Seperate Interaction intensity from P-value
    p_value = interaction_map.filter(regex='pvalue_')
    p_val_df = pd.concat(
        [interaction_map[['phenotype', 'neighbour_phenotype']], p_value],
        axis=1,
        join='outer',
    )
    p_val_df = p_val_df.set_index(['phenotype', 'neighbour_phenotype'])
    interaction_map = interaction_map[
        interaction_map.columns.difference(p_value.columns)
    ]
    interaction_map = interaction_map.set_index(['phenotype', 'neighbour_phenotype'])

    # Binarize the values if user requests
    if binary_view == True:
        interaction_map[interaction_map > 0] = 1
        interaction_map[interaction_map <= 0] = -1

    if summarize_plot == True:
        # convert first two columns to multi-index column
        # interaction_map = interaction_map.set_index(['phenotype','neighbour_phenotype'])
        # p_val_df = p_val_df.set_index(['phenotype','neighbour_phenotype'])

        # If multiple images are present, take the average of interactions
        interaction_map['mean'] = interaction_map.mean(axis=1).values
        interaction_map = interaction_map[['mean']]  # keep only the mean column
        interaction_map = interaction_map['mean'].unstack()
        # Do the same for P-values
        p_val_df['mean'] = p_val_df.mean(axis=1).values
        p_val_df = p_val_df[['mean']]  # keep only the mean column
        # set the P-value threshold
        p_val_df.loc[p_val_df[p_val_df['mean'] > p_val].index, 'mean'] = np.NaN
        p_val_df = p_val_df['mean'].unstack()

        # change to the order passed in subset
        if subset_phenotype is not None:
            interaction_map = interaction_map.reindex(subset_phenotype)
            p_val_df = p_val_df.reindex(subset_phenotype)
        if subset_neighbour_phenotype is not None:
            interaction_map = interaction_map.reindex(
                columns=subset_neighbour_phenotype
            )
            p_val_df = p_val_df.reindex(columns=subset_neighbour_phenotype)

        # Plotting heatmap
        mask = p_val_df.isnull()  # identify the NAN's for masking
        im = interaction_map.fillna(
            0
        )  # replace nan's with 0 so that clustering will work
        # heatmap
        plot = sns.clustermap(
            im,
            cmap=cmap_updated,
            row_cluster=row_cluster,
            col_cluster=col_cluster,
            mask=mask,
            **kwargs,
        )

    else:
        if len(interaction_map.columns) < 2:
            raise ValueError(
                'Data for only a single image is available please set summarize_plot=True and try again'
            )
        # convert first two columns to multi-index column
        # interaction_map = interaction_map.set_index(['phenotype','neighbour_phenotype'])
        # p_val_df = p_val_df.set_index(['phenotype','neighbour_phenotype'])

        # P value threshold
        p_val_df = p_val_df.apply(lambda x: np.where(x > p_val, np.nan, x))

        # Remove rows that are all nan
        idx = p_val_df.index[p_val_df.isnull().all(1)]  # Find all nan rows
        interaction_map = interaction_map.loc[
            interaction_map.index.difference(idx)
        ]  # clean intensity data
        p_val_df = p_val_df.loc[p_val_df.index.difference(idx)]  # clean p-value data

        # order the plot as needed
        if subset_phenotype or subset_neighbour_phenotype is not None:
            interaction_map.reset_index(inplace=True)
            p_val_df.reset_index(inplace=True)
            if subset_phenotype is not None:
                interaction_map['phenotype'] = (
                    interaction_map['phenotype'].astype('str').astype('category')
                )
                interaction_map['phenotype'] = interaction_map[
                    'phenotype'
                ].cat.reorder_categories(subset_phenotype)
                interaction_map = interaction_map.sort_values('phenotype')
                # Do same for Pval
                p_val_df['phenotype'] = (
                    p_val_df['phenotype'].astype('str').astype('category')
                )
                p_val_df['phenotype'] = p_val_df['phenotype'].cat.reorder_categories(
                    subset_phenotype
                )
                p_val_df = p_val_df.sort_values('phenotype')
            if subset_neighbour_phenotype is not None:
                interaction_map['neighbour_phenotype'] = (
                    interaction_map['neighbour_phenotype']
                    .astype('str')
                    .astype('category')
                )
                interaction_map['neighbour_phenotype'] = interaction_map[
                    'neighbour_phenotype'
                ].cat.reorder_categories(subset_neighbour_phenotype)
                interaction_map = interaction_map.sort_values('neighbour_phenotype')
                # Do same for Pval
                p_val_df['neighbour_phenotype'] = (
                    p_val_df['neighbour_phenotype'].astype('str').astype('category')
                )
                p_val_df['neighbour_phenotype'] = p_val_df[
                    'neighbour_phenotype'
                ].cat.reorder_categories(subset_neighbour_phenotype)
                p_val_df = p_val_df.sort_values('neighbour_phenotype')
            if subset_phenotype and subset_neighbour_phenotype is not None:
                interaction_map = interaction_map.sort_values(
                    ['phenotype', 'neighbour_phenotype']
                )
                p_val_df = p_val_df.sort_values(['phenotype', 'neighbour_phenotype'])

            # convert the data back into multi-index
            interaction_map = interaction_map.set_index(
                ['phenotype', 'neighbour_phenotype']
            )
            p_val_df = p_val_df.set_index(['phenotype', 'neighbour_phenotype'])

        # Plotting heatmap
        mask = p_val_df.isnull()  # identify the NAN's for masking
        im = interaction_map.fillna(
            0
        )  # replace nan's with 0 so that clustering will work
        mask.columns = im.columns

        # covert the first two columns into index
        # Plot
        plot = sns.clustermap(
            im,
            cmap=cmap_updated,
            row_cluster=row_cluster,
            col_cluster=col_cluster,
            mask=mask,
            **kwargs,
        )

    # Saving the figure if saveDir and fileName are provided
    if saveDir:
        if not os.path.exists(saveDir):
            os.makedirs(saveDir)
        full_path = os.path.join(saveDir, fileName)
        plot.savefig(full_path, dpi=300)
        plt.close()
        print(f"Saved plot to {full_path}")
    else:
        plt.show()

    if return_data is True:
        # perpare data for export
        map_data = interaction_map.copy()
        p_val_data = mask.copy()
        map_data.reset_index(inplace=True)
        p_val_data.reset_index(inplace=True)
        # remove the first two colums
        map_data = map_data.drop(['phenotype', 'neighbour_phenotype'], axis=1)
        p_val_data = p_val_data.drop(['phenotype', 'neighbour_phenotype'], axis=1)
        p_val_data.columns = map_data.columns
        # remove the mased values
        final_Data = map_data.where(~p_val_data, other=np.nan)
        final_Data.index = interaction_map.index
        return final_Data
