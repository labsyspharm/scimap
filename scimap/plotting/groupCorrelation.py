#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Mon Mar 18 08:48:29 2024
# @author: Ajit Johnson Nirmal

"""
!!! abstract "Short Description"
    The `sm.pl.groupCorrelation` function calculates and visualizes the correlation between group abundances across various conditions within an `AnnData` object. Customizable features such as normalization, hierarchical clustering, and manual ordering are available.

## Function
"""

# lib
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
import os
import anndata as ad
import warnings
from scipy.stats import zscore
import argparse


# plt.rcParams['figure.dpi'] = 100
# plt.rcParams['savefig.dpi']=300
# plt.rcParams['font.family']='sans serif'
# plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['pdf.fonttype'] = 42


# function
def groupCorrelation(
    adata,
    groupBy,
    condition,
    normalize=False,
    subsetGroups=None,
    orderRow=None,
    orderColumn=None,
    clusterRows=True,
    clusterColumns=True,
    cmap='vlag',
    figsize=None,
    overlayValues=False,
    fontSize=10,
    fontColor='black',
    fileName='groupCorrelation.pdf',
    saveDir=None,
    **kwargs,
):
    """
    Parameters:
        adata (AnnData or str):
            An AnnData object containing the dataset, or a string path to an AnnData file to be loaded.

        groupBy (str):
            The column in `adata.obs` used for defining groups.

        condition (str):
            The column in `adata.obs` that distinguishes different conditions or samples.

        normalize (bool, optional):
            If True, apply z-score normalization to the group counts across conditions.

        subsetGroups (list of str, optional):
            A list specifying a subset of groups to include in the analysis. If None, all groups are included.

        orderRow (list of str, optional):
            Custom order for the rows in the heatmap. If None, the order is determined by clustering or the original group order.

        orderColumn (list of str, optional):
            Custom order for the columns in the heatmap.

        clusterRows (bool, optional):
            Whether to apply hierarchical clustering to rows.

        clusterColumns (bool, optional):
            Whether to apply hierarchical clustering to columns.

        cmap (str, optional):
            The colormap for the heatmap.

        figsize (tuple of float, optional):
            The size of the figure to create (width, height). If None, the size is inferred.

        overlayValues (bool, optional):
            If True, overlays the correlation coefficient values on the heatmap.

        fontSize (int, optional):
            Font size for overlay values.

        fontColor (str, optional):
            Color of the font used for overlay values.

        fileName (str, optional):
            Name of the file to save the heatmap. Relevant only if `saveDir` is specified.

        saveDir (str, optional):
            Directory to save the generated heatmap. If None, the heatmap is not saved.

    Returns:
            plot (matplotlib):
                Displays or saves a heatmap visualizing the correlation between specified groups.

    Example:
        ```python

        # Basic usage with auto-detected conditions and groups
        sm.pl.groupCorrelation(adata, groupBy='cell_type', condition='patient_id')

        # Normalized group counts with specific groups and custom clustering disabled
        sm.pl.groupCorrelation(adata, groupBy='cell_type', condition='patient_id', normalize=True,
                         subsetGroups=['B cells', 'T cells'], clusterRows=False, clusterColumns=False)

        # Using custom ordering and overlaying values with specified font size and color
        sm.pl.groupCorrelation(adata, groupBy='cell_type', condition='patient_id', overlayValues=True,
                         orderRow=['T cells', 'B cells'], fontSize=12, fontColor='blue',
                         saveDir='/path/to/results', fileName='customGroupCorrelation.pdf')
        ```

    """

    # Load adata if a path is provided
    if isinstance(adata, str):
        adata = ad.read_h5ad(adata)

    # Calculate group counts
    group_counts = (
        adata.obs.groupby([condition, groupBy], observed=False)
        .size()
        .unstack(fill_value=0)
    )

    # Subset groups if needed
    if subsetGroups:
        group_counts = group_counts[subsetGroups]

    # Normalize if requested
    if normalize:
        group_counts = group_counts.apply(zscore, axis=0)

    # Calculate correlation
    corr_matrix = group_counts.corr()

    # var_names for axis labels, directly from group_counts columns
    var_names = group_counts.columns.tolist()

    # Manual ordering takes precedence over clustering
    if orderRow and clusterRows:
        warnings.warn(
            "Both orderRow and clusterRows were provided. Proceeding with orderRow and ignoring clusterRows."
        )
        clusterRows = False
    if orderColumn and clusterColumns:
        warnings.warn(
            "Both orderColumn and clusterColumns were provided. Proceeding with orderColumn and ignoring clusterColumns."
        )
        clusterColumns = False

    # Apply manual ordering or clustering
    if orderRow:
        row_order = [var_names.index(name) for name in orderRow]
    else:
        row_order = range(len(var_names))  # Default order if no manual ordering
        if clusterRows:
            linkage_row = linkage(pdist(corr_matrix, 'euclidean'), method='average')
            row_order = dendrogram(linkage_row, no_plot=True)['leaves']

    if orderColumn:
        col_order = [var_names.index(name) for name in orderColumn]
    else:
        col_order = range(len(var_names))  # Default order if no manual ordering
        if clusterColumns:
            linkage_col = linkage(pdist(corr_matrix.T, 'euclidean'), method='average')
            col_order = dendrogram(linkage_col, no_plot=True)['leaves']

    # Reorder the matrix based on row_order and col_order
    corr_matrix = corr_matrix.iloc[row_order, col_order]

    # Plotting
    if figsize is None:
        figsize_width = max(10, len(corr_matrix.columns) * 0.5)
        figsize_height = max(8, len(corr_matrix.index) * 0.5)
        figsize = (figsize_width, figsize_height)

    plt.figure(figsize=figsize)
    im = plt.imshow(corr_matrix, cmap=cmap, aspect='auto', **kwargs)
    plt.colorbar(im)

    if overlayValues:
        for i in range(len(row_order)):
            for j in range(len(col_order)):
                plt.text(
                    j,
                    i,
                    f"{corr_matrix.iloc[i, j]:.2f}",
                    ha="center",
                    va="center",
                    color=fontColor,
                    fontsize=fontSize,
                )

    # Set tick labels
    plt.xticks(
        ticks=np.arange(len(col_order)),
        labels=[var_names[i] for i in col_order],
        rotation=90,
    )
    plt.yticks(
        ticks=np.arange(len(row_order)), labels=[var_names[i] for i in row_order]
    )

    plt.tight_layout()

    # Save or show the figure
    if saveDir:
        if not os.path.exists(saveDir):
            os.makedirs(saveDir)
        full_path = os.path.join(saveDir, fileName)
        plt.savefig(full_path, dpi=300)
        plt.close()
        print(f"Saved heatmap to {full_path}")
    else:
        plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Calculates and visualizes the correlation between group abundances across various conditions within an AnnData object.'
    )

    parser.add_argument(
        '--adata',
        type=str,
        required=True,
        help='Path to an AnnData object file containing the dataset to be visualized.',
    )
    parser.add_argument(
        '--groupBy',
        type=str,
        required=True,
        help="The column in `adata.obs` used for defining groups.",
    )
    parser.add_argument(
        '--condition',
        type=str,
        required=True,
        help="The column in `adata.obs` that distinguishes different conditions or samples.",
    )
    parser.add_argument(
        '--normalize',
        action='store_true',
        help="Apply z-score normalization to the group counts across conditions. Defaults to False.",
    )
    parser.add_argument(
        '--subsetGroups',
        type=str,
        nargs='+',
        default=None,
        help="A list specifying a subset of groups to include in the analysis.",
    )
    parser.add_argument(
        '--orderRow',
        type=str,
        nargs='+',
        default=None,
        help="Custom order for the rows in the heatmap.",
    )
    parser.add_argument(
        '--orderColumn',
        type=str,
        nargs='+',
        default=None,
        help="Custom order for the columns in the heatmap. Equivalent to `orderRow` due to the square nature of the correlation matrix but kept for consistency.",
    )
    parser.add_argument(
        '--clusterRows',
        action='store_true',
        help="Whether to apply hierarchical clustering to rows. Defaults to True unless --no-clusterRows is specified.",
    )
    parser.add_argument(
        '--no-clusterRows',
        action='store_false',
        dest='clusterRows',
        help="Do not cluster rows.",
    )
    parser.add_argument(
        '--clusterColumns',
        action='store_true',
        help="Whether to apply hierarchical clustering to columns. Defaults to True unless --no-clusterColumns is specified.",
    )
    parser.add_argument(
        '--no-clusterColumns',
        action='store_false',
        dest='clusterColumns',
        help="Do not cluster columns.",
    )
    parser.add_argument(
        '--cmap',
        type=str,
        default='vlag',
        help="The colormap for the heatmap. Defaults to 'vlag'.",
    )
    parser.add_argument(
        '--figsize',
        type=float,
        nargs=2,
        default=None,
        help="The size of the figure to create (width, height). If None, the size is inferred.",
    )
    parser.add_argument(
        '--overlayValues',
        action='store_true',
        help="Overlay the actual correlation values on the heatmap.",
    )
    parser.add_argument(
        '--fontSize',
        type=int,
        default=10,
        help="Font size for overlay values. Defaults to 10.",
    )
    parser.add_argument(
        '--fontColor',
        type=str,
        default='black',
        help="Color of the font used for overlay values. Defaults to 'black'.",
    )
    parser.add_argument(
        '--fileName',
        type=str,
        default='groupCorrelation.pdf',
        help="Name of the file to save the heatmap. Defaults to 'groupCorrelation.pdf'.",
    )
    parser.add_argument(
        '--saveDir',
        type=str,
        default=None,
        help="Directory to save the generated heatmap. If None, the heatmap is not saved.",
    )

    args = parser.parse_args()

    # Execute groupCorrelation with the provided arguments
    groupCorrelation(
        adata=args.adata,
        groupBy=args.groupBy,
        condition=args.condition,
        normalize=args.normalize,
        subsetGroups=args.subsetGroups,
        orderRow=args.orderRow,
        orderColumn=args.orderColumn,
        clusterRows=args.clusterRows,
        clusterColumns=args.clusterColumns,
        cmap=args.cmap,
        figsize=args.figsize,
        overlayValues=args.overlayValues,
        fontSize=args.fontSize,
        fontColor=args.fontColor,
        fileName=args.fileName,
        saveDir=args.saveDir,
    )
