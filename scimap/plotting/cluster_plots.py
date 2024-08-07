#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
!!! abstract "Short Description"
    `sm.pl.cluster_plots`: This versatile function streamlines the visualization 
    process by generating UMAP plots, heatmaps of the expression matrix, and 
    lists of ranked marker genes for each user-defined group, typically following 
    clustering analysis via `sm.tl.cluster`. It offers a comprehensive overview of 
    clustering results, facilitating the exploration of spatial patterns, 
    molecular profiles, and key markers distinguishing each cluster.

## Function
"""

# Import
import argparse
import sys
import pathlib
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns; sns.set(color_codes=True)
sns.set_style("white")

plt.rcParams['pdf.fonttype'] = 42

def main(argv=sys.argv):
    parser = argparse.ArgumentParser(
        description='Helper function that allows users to save the contents of the `scimap` object as a csv file.'
    )
    parser.add_argument(
        '--adata', required=True, 
        help='AnnData object loaded into memory or path to AnnData object.'
    )
    parser.add_argument(
        '--group_by', type=str, required=True,
        help='Name of the categorical column that contains the clustering results'
    )
    parser.add_argument(
        '--subsample', type=int, required=False, default='100000',
        help='Subsample number of observations.'
    )
    parser.add_argument(
        '--palette', type=str, required=False, default='tab10',
        help='Colors to use for plotting categorical annotation groups.'
    )
    parser.add_argument(
        '--size', type=int, required=False, default=None,
        help='Point size of Umap'
    )
    parser.add_argument(
        '--use_raw', type=int, required=False, default=False,
        help='Use .raw attribute of adata for coloring the matrixplot expression matrix'
    )
    parser.add_argument(
        '--output_dir', type=str, required=False, default=None,
        help='Path to output directory.'
    )
    args = parser.parse_args(argv[1:])
    print(vars(args))
    cluster_plots(**vars(args))

# Function
def cluster_plots (adata, 
                   group_by, 
                   subsample=100000, 
                   palette ='viridis', 
                   use_raw=False,
                   size=None, 
                   output_dir=None):
    """
Parameters:
        adata (anndata.AnnData):  
            The annotated data matrix.

        group_by (str):  
            The column name in `adata.obs` that contains the clustering labels to visualize.

        subsample (int, optional):  
            The number of cells to randomly subsample from the dataset for visualization to enhance performance. 
            Default is 100000. If set to None, no subsampling is performed.

        palette (str, optional):  
            The name of a matplotlib colormap to use for coloring clusters. Default is 'viridis'.

        use_raw (bool, optional):  
            If True, uses the `.raw` attribute of `adata` for extracting expression data for the matrix plot. 
            Default is False.

        size (int, optional):  
            The size of the points in the UMAP plot. Default is 40.

        output_dir (str, optional):  
            The directory where the plots should be saved. If not specified, plots are shown but not saved.

Returns:
        plots (matplotlib):  
            The function does not return a value but generates and optionally saves the specified plots.

Example:
    ```python
    
    # Generate cluster plots using default settings
    sm.pl.cluster_plots(adata, group_by='leiden')

    # Generate cluster plots with a custom palette and subsampling
    sm.pl.cluster_plots(adata, group_by='leiden', palette='plasma', subsample=50000)

    # Generate cluster plots without subsampling, using raw data, and save them to a directory
    sm.pl.cluster_plots(adata, group_by='leiden', subsample=None, use_raw=True, output_dir='./cluster_plots')
    
    ```
    """
    
    # Load the data 
    if isinstance(adata, str):
        imid = pathlib.Path(adata).stem
        adata = ad.read(adata)  
    else:
        adata = adata
        imid = ""
        
    # Subset data if needed
    if subsample is not None:
        if adata.shape[0] > subsample:
            sc.pp.subsample(adata, n_obs=subsample)
    
    
    # UMAP
    try:
        sc.pp.neighbors(adata) # Computing the neighborhood graph
        sc.tl.umap(adata)
        fig = sc.pl.umap(adata, color=group_by, palette = palette, size=size, return_fig=True, show=False) # View the clustering
        fig.tight_layout()
        # save figure
        if output_dir is not None:
            output_dir = pathlib.Path(output_dir)
            output_dir.mkdir(exist_ok=True, parents=True)
            #fig.savefig(output_dir / f"{imid}_umap.pdf")
            fig.savefig(pathlib.Path(output_dir) / f"{imid}_umap.pdf")
            
    except Exception as exc:
        print('UMAP could not be generated')
        print (exc)

    # Matrix plot
    try:
        mat_fig = sc.pl.matrixplot(adata, var_names=adata.var.index, groupby=group_by, use_raw=use_raw,
                         cmap='RdBu_r', dendrogram=True, title = group_by,
                         return_fig=True
                         )
        if output_dir is not None:
            #mat_fig.savefig(output_dir / 'matrixplot.pdf')
            mat_fig.savefig(pathlib.Path(output_dir) / f"{imid}_matrixplot.pdf")
            
    except Exception as exc:
        print('Heatmap could not be generated')
        print (exc)
    
    # Marker expression per group
    try:
        sc.tl.rank_genes_groups(adata, group_by, method='t-test')
     
        # find number of genes in dataset
        if len(adata.var.index) > 20:
            n_genes = 20
        else:
            n_genes = len(adata.var.index)
        
        if output_dir is not None:
            sc.pl.rank_genes_groups(adata, sharey=False, n_genes=n_genes, fontsize=12, show=False)
            plt.suptitle(group_by, fontsize=20)
            #plt.savefig(output_dir / 'ranked_markers_per_cluster.pdf')
            plt.savefig(pathlib.Path(output_dir) / f"{imid}_ranked_markers_per_cluster.pdf")
        else:
            sc.pl.rank_genes_groups(adata, sharey=False, n_genes=n_genes, fontsize=12)
            
    except Exception as exc:
        print('Finding differential markers per group cannot be completed')
        print (exc)


if __name__ == '__main__':
    main()  
