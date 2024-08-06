#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Wed Apr  7 17:46:29 2021
# @author: Ajit Johnson Nirmal
"""
!!! abstract "Short Description"
    `sm.pl.foldchange`: This function facilitates the visualization of fold changes in cell type 
    abundance across samples or Regions of Interest (ROIs), offering insights into differential 
    expression or abundance patterns. It is designed to work with data processed by `sm.tl.foldchange`, 
    which should be executed beforehand to calculate the fold changes. Through heatmap or parallel 
    coordinates presentations, users can effectively interpret variations, highlighting significant 
    shifts and guiding further analyses.

## Function
"""

# lib
import seaborn as sns

sns.set(color_codes=True)
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from pandas.plotting import parallel_coordinates
import os

sns.set_style("white")
plt.rcParams['pdf.fonttype'] = 42


# Function
def foldchange(
    adata,
    label='foldchange',
    p_val=0.05,
    nonsig_color='grey',
    subset_xaxis=None,
    subset_yaxis=None,
    cmap='vlag',
    log=True,
    center=0,
    method='heatmap',
    invert_axis=None,
    parallel_coordinates_color=None,
    matplotlib_bbox_to_anchor=(1.04, 1),
    matplotlib_legend_loc='upper left',
    xticks_rotation=90,
    return_data=False,
    fileName='foldchange.pdf',
    saveDir=None,
    **kwargs,
):
    """
    Parameters:
            adata (anndata.AnnData):
                The annotated data matrix with fold change calculations.

            label (str):
                Label key from `adata.uns` indicating the fold change data to visualize.

            p_val (float):
                P-value threshold for significance; values above this threshold are considered non-significant.

            nonsig_color (str):
                Color for non-significant changes in the visualization.

            subset_xaxis (list, optional):
                Categories to include from the x-axis for plotting.

            subset_yaxis (list, optional):
                Categories to include from the y-axis for plotting.

            cmap (str):
                Colormap for the heatmap visualization.

            log (bool):
                If True, fold changes are converted to log2 scale for visualization.

            center (float):
                The value at which the colormap is centered.

            method (str):
                Plotting method, either 'heatmap' for a heatmap or 'parallel_coordinates' for a parallel coordinates plot.

            invert_axis (bool, optional):
                If True, inverts the x and y axes in the plot. Default is False.

            parallel_coordinates_color (list, optional):
                Specifies custom colors for parallel coordinates plot.

            matplotlib_bbox_to_anchor (tuple):
                Adjusts the legend position in parallel coordinates plot.

            matplotlib_legend_loc (str):
                Specifies the location of the legend.

            xticks_rotation (int):
                Rotation angle for x-axis tick labels.

            return_data (bool):
                If True, returns the data frame used for plotting instead of the plot.

    Returns:
            Dataframe; plot (pandas, matplotlib):
                If `return_data` is True, returns a pandas DataFrame used for plotting. Otherwise, displays the plot.

    Example:
            ```python

            # Generate a heatmap of fold changes with custom settings
            sm.pl.foldchange(adata, label='foldchange', method='heatmap', cmap='coolwarm', log=True,
                             p_val=0.05, nonsig_color='lightgrey', xticks_rotation=45)

            # Create a parallel coordinates plot to visualize fold changes across groups
            sm.pl.foldchange(adata, label='foldchange', method='parallel_coordinates', log=True,
                             parallel_coordinates_color=['red', 'blue', 'green'], invert_axis=True)

            # Return the data frame used for fold change visualization
            df_foldchange = sm.pl.foldchange(adata, label='foldchange', return_data=True)

            ```
    """

    # set color for heatmap
    # cmap_updated = copy.copy(matplotlib.cm.get_cmap(cmap))
    cmap_updated = matplotlib.cm.get_cmap(cmap)
    cmap_updated.set_bad(color=nonsig_color)

    # get the data
    fc = adata.uns[str(label) + '_fc']
    p = adata.uns[str(label) + '_pval']

    # fold
    fold = fc.copy()
    p_mask = p.copy()

    # reference image
    ref = fold.index.name

    # log
    if log is True:
        fold = np.log2(fold)

    # create a mask for non-sig values
    p_mask[p_mask > p_val] = np.nan

    # subset x axis data
    if subset_xaxis is not None:
        if isinstance(subset_xaxis, str):
            subset_xaxis = [subset_xaxis]
        fold = fold[subset_xaxis]
        p_mask = p_mask[subset_xaxis]
        # reorder

    # subset y axis data
    if subset_yaxis is not None:
        if isinstance(subset_yaxis, str):
            subset_yaxis = [subset_yaxis]
        fold = fold.loc[subset_yaxis]
        p_mask = p_mask.loc[subset_yaxis]
        # reorder

    # invert axis if user requests
    if invert_axis is True:
        fold = fold.T
        p_mask = p_mask.T

    # mask
    mask = p_mask.isnull()  # identify the NAN's for masking

    if method == 'heatmap':
        # heatmap of the foldchange
        # g= sns.clustermap(fold, cmap=cmap, mask=mask, center=center, col_cluster=False, row_cluster=False)
        fig = sns.clustermap(fold, cmap=cmap, mask=mask, center=center, **kwargs)
        plt.suptitle('reference: ' + str(ref))
        plt.setp(fig.ax_heatmap.get_xticklabels(), rotation=xticks_rotation)
        plt.tight_layout()

    if method == 'parallel_coordinates':
        fold['sample'] = fold.index
        # plotting
        fig, axes = plt.subplots()
        if parallel_coordinates_color is not None:
            parallel_coordinates(
                fold, 'sample', color=parallel_coordinates_color, **kwargs
            )
        else:
            # parallel_coordinates(fold, 'sample', colormap=cmap_updated)
            parallel_coordinates(fold, 'sample', colormap=cmap_updated, **kwargs)
        axes.grid(False)
        plt.legend(bbox_to_anchor=matplotlib_bbox_to_anchor, loc=matplotlib_legend_loc)
        plt.axhline(y=0, color='black', linestyle='-')
        plt.xticks(rotation=xticks_rotation)
        plt.suptitle('reference: ' + str(ref))
        fig.tight_layout()

    # save
    if saveDir:
        if not os.path.exists(saveDir):
            os.makedirs(saveDir)
        full_path = os.path.join(saveDir, fileName)
        plt.savefig(full_path, dpi=300)
        plt.close()
        print(f"Saved plot to {full_path}")
    else:
        plt.show()

    # return data
    if return_data is True:
        return fold
