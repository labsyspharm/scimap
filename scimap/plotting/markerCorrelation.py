#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Mon Mar 18 08:48:29 2024
# @author: Ajit Johnson Nirmal


"""
!!! abstract "Short Description"
    The `sm.pl.markerCorrelation` function computes and visualizes the correlation among selected markers (genes, proteins, etc.) within an `AnnData` object.

## Function
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
import os
import anndata as ad
import warnings
import argparse


# plt.rcParams['figure.dpi'] = 100
# plt.rcParams['savefig.dpi']=300
# plt.rcParams['font.family']='sans serif'
# plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['pdf.fonttype'] = 42


def markerCorrelation(
    adata,
    layer='log',
    subsetMarkers=None,
    orderRow=None,
    orderColumn=None,
    clusterRows=True,
    clusterColumns=True,
    cmap='vlag',
    figsize=None,
    overlayValues=False,
    fontSize=10,
    fontColor='black',
    fileName='markerCorrelation.pdf',
    saveDir=None,
    **kwargs,
):
    """
    Parameters:
            adata (AnnData or str):
                An AnnData object containing the dataset, or a string path to an AnnData file to be loaded.

            layer (str, optional):
                Specifies the layer of `adata` to use for the heatmap. If None, the `.X` attribute is used. If you want to plot the raw data use `raw`

            subsetMarkers (list of str, optional):
                A list of marker names to include in the correlation analysis. If None, all markers are used.

            orderRow (list of str, optional):
                Specifies a custom order for the rows (markers) based on their names.

            orderColumn (list of str, optional):
                Specifies a custom order for the columns (markers) based on their names.

            clusterRows (bool, optional):
                Whether to apply hierarchical clustering to rows.

            clusterColumns (bool, optional):
                Whether to apply hierarchical clustering to columns.

            cmap (str, optional):
                The colormap for the heatmap.

            figsize (tuple of float, optional):
                The size of the figure to create. If None, the size is inferred based on the data.

            overlayValues (bool, optional):
                If True, overlays the actual correlation values on the heatmap.

            fontSize (int, optional):
                Font size for the overlay values.

            fontColor (str, optional):
                Color of the font used for overlay values.

            fileName (str, optional):
                Name of the file to save the heatmap. Relevant only if `saveDir` is not None.

            saveDir (str, optional):
                Directory to save the generated heatmap. If None, the heatmap is not saved.

    Returns:
            plot (matplotlib):
                Displays or saves a heatmap visualizing the correlation between specified markers.

    Example:
        ```python

            # Example 1: Basic usage with all markers and default parameters
            sm.pl.markerCorrelation(adata)

            # Example 2: With subset of markers, custom clustering, and overlaying correlation values
            sm.pl.markerCorrelation(adata, subsetMarkers=['Marker1', 'Marker2', 'Marker3'], clusterRows=False, overlayValues=True, fontSize=12)

            # Example 3: Saving the heatmap to a specific directory
            sm.pl.markerCorrelation(adata, fileName='myHeatmap.pdf', saveDir='/path/to/save')

        ```
    """

    # load adata
    if isinstance(adata, str):
        adata = ad.read_h5ad(adata)

    # subset the markers if user requests
    if subsetMarkers:
        subsetMarkers = (
            [subsetMarkers] if isinstance(subsetMarkers, str) else subsetMarkers
        )  # convert to list
        # isolate the data
        if layer == 'raw':
            matrix = adata[:, subsetMarkers].raw.X
        elif layer is None:
            matrix = adata[:, subsetMarkers].X
        else:
            matrix = adata[:, subsetMarkers].layers[layer]
    else:
        # take the whole data if the user does not subset anything
        if layer == 'raw':
            matrix = adata.raw.X
        elif layer is None:
            matrix = adata.X
        else:
            matrix = adata.layers[layer]

    # intialize the markers to be plotted
    if subsetMarkers is None:
        var_names = adata.var_names.tolist()
    else:
        var_names = subsetMarkers

    # run correlation
    corr_matrix = np.corrcoef(matrix.T)

    row_order = np.arange(corr_matrix.shape[0])
    col_order = np.arange(corr_matrix.shape[1])

    if orderRow:
        if clusterRows:
            warnings.warn(
                "Both orderRow and clusterRows were provided. Proceeding with orderRow and ignoring clusterRows."
            )
            clusterRows = False
        row_order = [var_names.index(name) for name in orderRow]

    if orderColumn:
        if clusterColumns:
            warnings.warn(
                "Both orderColumn and clusterColumns were provided. Proceeding with orderColumn and ignoring clusterColumns."
            )
            clusterColumns = False
        col_order = [var_names.index(name) for name in orderColumn]

    corr_matrix = corr_matrix[np.ix_(row_order, col_order)]

    if clusterRows:
        linkage_row = linkage(pdist(corr_matrix), method='average')
        row_order = dendrogram(linkage_row, no_plot=True)['leaves']
        corr_matrix = corr_matrix[row_order, :]

    if clusterColumns:
        linkage_col = linkage(pdist(corr_matrix.T), method='average')
        col_order = dendrogram(linkage_col, no_plot=True)['leaves']
        corr_matrix = corr_matrix[:, col_order]

    if figsize is None:
        base_size = 0.5  # Base size for each cell in inches
        figsize_width = max(10, len(corr_matrix) * base_size)
        figsize_height = max(8, len(corr_matrix) * base_size)
        figsize = (figsize_width, figsize_height)

    plt.figure(figsize=figsize)
    im = plt.imshow(corr_matrix, cmap=cmap, aspect='auto', **kwargs)
    plt.colorbar(im)

    if overlayValues:
        for i in range(corr_matrix.shape[0]):
            for j in range(corr_matrix.shape[1]):
                text = plt.text(
                    j,
                    i,
                    f"{corr_matrix[i, j]:.2f}",
                    ha="center",
                    va="center",
                    color=fontColor,
                    fontsize=fontSize,
                )

    plt.xticks(
        ticks=np.arange(len(col_order)),
        labels=np.array(var_names)[col_order],
        rotation=90,
    )
    plt.yticks(ticks=np.arange(len(row_order)), labels=np.array(var_names)[row_order])
    plt.tight_layout()

    # Saving the figure if saveDir and fileName are provided
    if saveDir:
        if not os.path.exists(saveDir):
            os.makedirs(saveDir)
        full_path = os.path.join(saveDir, fileName)
        plt.savefig(full_path, dpi=300)
        plt.close()
        print(f"Saved plot to {full_path}")
    else:
        plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compute and visualize the correlation among markers within an AnnData object.'
    )

    parser.add_argument(
        '--adata',
        type=str,
        required=True,
        help='Path to an AnnData object file containing the dataset to be visualized.',
    )
    parser.add_argument(
        '--layer',
        type=str,
        default='log',
        help="Specifies the layer of `adata` to use for the correlation analysis. Defaults to 'log'.",
    )
    parser.add_argument(
        '--subsetMarkers',
        type=str,
        nargs='+',
        default=None,
        help="A list of marker genes or features to include in the correlation analysis.",
    )
    parser.add_argument(
        '--orderRow',
        type=str,
        nargs='+',
        default=None,
        help="Custom order for the rows based on marker names.",
    )
    parser.add_argument(
        '--orderColumn',
        type=str,
        nargs='+',
        default=None,
        help="Custom order for the columns based on marker names.",
    )
    parser.add_argument(
        '--clusterRows',
        action='store_true',
        help="Whether to cluster rows. Defaults to True unless --no-clusterRows is specified.",
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
        help="Whether to cluster columns. Defaults to True unless --no-clusterColumns is specified.",
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
        help="The size of the figure to create. Specify width and height.",
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
        help="Font size for the overlay values. Defaults to 10.",
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
        default='markerCorrelation.pdf',
        help="Name of the file to save the heatmap. Defaults to 'markerCorrelation.pdf'.",
    )
    parser.add_argument(
        '--saveDir',
        type=str,
        default=None,
        help="Directory to save the generated heatmap. If None, the heatmap is not saved.",
    )

    args = parser.parse_args()

    # Execute markerCorrelation with the provided arguments
    markerCorrelation(
        adata=args.adata,
        layer=args.layer,
        subsetMarkers=args.subsetMarkers,
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
