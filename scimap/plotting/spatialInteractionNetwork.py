#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Tue Mar 19 09:11:11 2024
# @author: Ajit Johnson Nirmal

"""
!!! abstract "Short Description"
    The `sm.pl.spatialInteractionNetwork` function visualizes spatial interactions as a network graph, highlighting significant interactions between cell types.

## Function
"""

# Libs
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import anndata as ad
from scipy.stats import combine_pvalues
import matplotlib.pyplot as plt
import os
import argparse


# plt.rcParams['figure.dpi'] = 100
# plt.rcParams['savefig.dpi']=300
# plt.rcParams['font.family']='sans serif'
# plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['pdf.fonttype'] = 42


# Function


def spatialInteractionNetwork(
    adata,
    spatial_interaction='spatial_interaction',
    p_val=0.05,
    cmap='vlag',
    nodeColor='#22333b',
    nodeSize=None,
    alpha=0.9,
    figsize=None,
    fontSize=12,
    fontColor='white',
    subsetPhenotype=None,
    fileName='spatialInteractionNetwork.pdf',
    saveDir=None,
):
    """
    Parameters:
        adata (AnnData or str):
            An AnnData object or a path to an .h5ad file containing the dataset.

        spatial_interaction (str):
            Key in `adata.uns` for spatial interaction data.

        p_val (float):
            Threshold for significance of interactions to display.

        cmap (str):
            Colormap for the edges based on their z-scores.

        nodeColor (str):
            Color of the nodes.

        nodeSize (int or None):
            Size of the nodes. If None, size is dynamically adjusted.

        alpha (float):
            Opacity of the nodes.

        figsize (tuple or None):
            Figure size as (width, height). If None, size is dynamically calculated.

        fontSize (int):
            Font size for node labels.

        fontColor (str):
            Color of the node labels.

        subsetPhenotype (list of str or None):
            List of phenotypes to include. If None, all are included.

        fileName (str):
            Filename for saving the network plot.

        saveDir (str or None):
            Directory to save the plot file. If None, plot is shown and not saved.

    Returns:
            plot (matplotlib):
                Displays or saves a network plot visualizing the interactions between cell types.

    Example:
        ```python

        # Visualizes the network using the 'coolwarm' colormap.
        sm.pl.spatialInteractionNetwork(adata, cmap='coolwarm')

        # Filters for 'T cells' and 'B cells' interactions, saves the visualization as 'T_B_interaction.pdf' in './plots'.
        sm.pl.spatialInteractionNetwork(adata, subsetPhenotype=['T cells', 'B cells'], fileName='T_B_interaction.pdf', saveDir='./plots')
        ```
    """

    # Load adata if a path is provided
    if isinstance(adata, str):
        adata = ad.read_h5ad(adata)

    # create a copy of the distance measurement
    if spatial_interaction not in adata.uns:
        raise KeyError(
            f"{spatial_interaction} does not exist in adata.uns. Please check if the '{spatial_interaction}' column exists or run `sm.tl.spatial_interaction(adata)` to compute it."
        )

    # copy the data to plot
    df = adata.uns[spatial_interaction].copy()

    # subset
    if subsetPhenotype:
        if isinstance(subsetPhenotype, str):
            subsetPhenotype = [subsetPhenotype]
        df = df[
            df['phenotype'].isin(subsetPhenotype)
            & df['neighbour_phenotype'].isin(subsetPhenotype)
        ]
        # Convert to categorical if not already
        df['phenotype'] = df['phenotype'].astype('str').astype('category')
        df['neighbour_phenotype'] = (
            df['neighbour_phenotype'].astype('str').astype('category')
        )

    # now calculate a meta score across images
    # Automatically identify z-score and p-value columns
    z_score_columns = [
        col
        for col in df.columns
        if 'pvalue_' not in col and col not in ['phenotype', 'neighbour_phenotype']
    ]
    p_value_columns = [col for col in df.columns if 'pvalue_' in col]

    # Ensure there is a matching p-value column for each z-score column
    assert len(z_score_columns) == len(
        p_value_columns
    ), "The number of z-score columns does not match the number of p-value columns."

    def stouffers_method(z_scores):
        """Combines z-scores using Stouffer's method."""
        combined_z = np.sum(z_scores) / np.sqrt(len(z_scores))
        return combined_z

    def combine_z_scores(row):
        """Extracts and combines z-scores from the row."""
        z_scores = row[z_score_columns].values
        combined_z = stouffers_method(z_scores)
        return combined_z

    def combine_p_values(row, p_value_columns):
        """Extracts p-values from the row and combines them using Fisher's method."""
        p_values = [
            row[col] for col in p_value_columns
        ]  # Extract p-values for the specified columns

        # Convert to a NumPy array and ensure type float for NaN handling
        p_values_array = np.array(p_values, dtype=float)

        # Filter out NaN values before combining
        p_values_filtered = p_values_array[~np.isnan(p_values_array)]

        # Check if filtered array is empty
        if p_values_filtered.size == 0:
            return np.nan  # Return NaN if there are no valid p-values to combine

        # Combine p-values using Fisher's method
        _, combined_p = combine_pvalues(p_values_filtered)
        return combined_p

    # Combine z-scores and p-values for each interaction
    df['combined_z'] = df.apply(combine_z_scores, axis=1)
    df['combined_p'] = df.apply(
        lambda row: combine_p_values(row, p_value_columns), axis=1
    )

    # Create a consolidated DataFrame with the relevant columns
    df = df[['phenotype', 'neighbour_phenotype', 'combined_z', 'combined_p']]
    df.columns = ['phenotype', 'neighbour_phenotype', 'z score', 'p-value']

    # Filter out edges with p-value >= 0.05
    df_filtered = df[df['p-value'] < p_val]

    # Create a directed graph from the filtered DataFrame
    G = nx.from_pandas_edgelist(
        df_filtered,
        'phenotype',
        'neighbour_phenotype',
        ['z score', 'p-value'],
        create_using=nx.DiGraph(),
    )

    # Normalize z-scores for edge color
    z_scores = nx.get_edge_attributes(G, 'z score')
    min_z = min(z_scores.values())
    max_z = max(z_scores.values())
    # Apply normalization for coloring
    colors = [
        plt.cm.coolwarm((z_scores[edge] - min_z) / (max_z - min_z))
        for edge in G.edges()
    ]

    # Normalize p-values for edge thickness
    p_values = nx.get_edge_attributes(G, 'p-value')
    # Invert and normalize p-values to range for thickness: Higher values for lower p-values
    min_p = min(p_values.values())
    max_p = max(p_values.values())
    thicknesses = [
        10 * (1 - (p_values[edge] - min_p) / (max_p - min_p)) + 1 for edge in G.edges()
    ]

    # Use spring_layout considering the 'weight' for layout
    pos = nx.spring_layout(G, weight='weight')

    if figsize is None:
        # Adjust these scaling factors to suit your specific needs
        figsize_width_scale = 0.5
        figsize_height_scale = 0.5
        # Calculate width and height based on the DataFrame dimensions
        figsize_width = max(10, len(df.columns) * figsize_width_scale)
        figsize_height = max(8, len(df) * figsize_height_scale)
        figsize = (figsize_width, figsize_height)

    # Base node size that works well for a small number of nodes
    if nodeSize is None:
        node_count = G.number_of_nodes()
        nodeSize = 2000 / (
            node_count / 10
        )  # Example scaling, adjust the divisor as needed

    # Drawing the network graph
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    # Draw the network components
    nx.draw_networkx_nodes(
        G, pos, ax=ax, node_size=nodeSize, node_color=nodeColor, alpha=alpha
    )
    nx.draw_networkx_edges(
        G, pos, ax=ax, width=thicknesses, edge_color=colors, arrowstyle='->'
    )
    nx.draw_networkx_labels(
        G,
        pos,
        ax=ax,
        font_size=fontSize,
        font_family="sans-serif",
        font_color=fontColor,
    )

    # Setup the ScalarMappable for the colorbar reflecting z-scores
    sm = plt.cm.ScalarMappable(
        cmap=plt.get_cmap(cmap), norm=plt.Normalize(vmin=min_z, vmax=max_z)
    )
    # sm = plt.cm.ScalarMappable(cmap=plt.cm.coolwarm, norm=plt.Normalize(vmin=min_z, vmax=max_z))
    sm.set_array([])
    cbar = fig.colorbar(
        sm, ax=ax, orientation='vertical', fraction=0.046, pad=0.04, aspect=10
    )
    cbar.set_label('Z-score')

    # Since matplotlib's colorbar does not directly support displaying edge thickness, you might add a text or legend
    # describing the mapping of p-values to thickness if necessary.
    ax.axis('off')

    # Save or show the figure
    if saveDir and fileName:
        if not os.path.exists(saveDir):
            os.makedirs(saveDir)
        full_path = os.path.join(saveDir, fileName)
        plt.savefig(full_path, dpi=300)
        if not os.path.exists(saveDir):
            os.makedirs(saveDir)
        print(f"Saved network plot to {full_path}")
    else:
        plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Visualizes spatial interactions within single-cell data as a network graph.'
    )

    parser.add_argument(
        '--adata',
        type=str,
        required=True,
        help='Path to an AnnData object file containing the dataset.',
    )
    parser.add_argument(
        '--spatial_interaction',
        type=str,
        default='spatial_interaction',
        help="Key in `adata.uns` for spatial interaction data. Defaults to 'spatial_interaction'.",
    )
    parser.add_argument(
        '--p_val',
        type=float,
        default=0.05,
        help="P-value threshold for filtering interactions. Defaults to 0.05.",
    )
    parser.add_argument(
        '--cmap',
        type=str,
        default='vlag',
        help="Colormap for the edges based on their z-scores. Defaults to 'vlag'.",
    )
    parser.add_argument(
        '--nodeColor',
        type=str,
        default='#22333b',
        help="Color of the nodes. Defaults to '#22333b'.",
    )
    parser.add_argument(
        '--nodeSize',
        type=int,
        default=None,
        help="Size of the nodes. If None, size is dynamically adjusted. Defaults to None.",
    )
    parser.add_argument(
        '--alpha',
        type=float,
        default=0.9,
        help="Opacity of the nodes. Defaults to 0.9.",
    )
    parser.add_argument(
        '--figsize',
        type=float,
        nargs=2,
        default=None,
        help="Figure size as (width, height). If None, size is dynamically calculated. Defaults to None.",
    )
    parser.add_argument(
        '--fontSize',
        type=int,
        default=12,
        help="Font size for node labels. Defaults to 12.",
    )
    parser.add_argument(
        '--fontColor',
        type=str,
        default='white',
        help="Color of the node labels. Defaults to 'white'.",
    )
    parser.add_argument(
        '--subsetPhenotype',
        type=str,
        nargs='+',
        default=None,
        help="List of phenotypes to include. If None, all are included. Defaults to None.",
    )
    parser.add_argument(
        '--fileName',
        type=str,
        default='spatialInteractionNetwork.pdf',
        help="Filename for saving the network plot. Defaults to 'spatialInteractionNetwork.pdf'.",
    )
    parser.add_argument(
        '--saveDir',
        type=str,
        default=None,
        help="Directory to save the plot file. If None, plot is shown and not saved. Defaults to None.",
    )

    args = parser.parse_args()

    # Call spatialInteractionNetwork with the parsed arguments
    spatialInteractionNetwork(
        adata=args.adata,
        spatial_interaction=args.spatial_interaction,
        p_val=args.p_val,
        cmap=args.cmap,
        nodeColor=args.nodeColor,
        nodeSize=args.nodeSize,
        alpha=args.alpha,
        figsize=args.figsize,
        fontSize=args.fontSize,
        fontColor=args.fontColor,
        subsetPhenotype=args.subsetPhenotype,
        fileName=args.fileName,
        saveDir=args.saveDir,
    )
