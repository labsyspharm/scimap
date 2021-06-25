#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
!!! abstract "Short Description"
    `sm.pl.cluster_plots`: A quick meta function that outputs umap plots and heatmap of the expression matrix
    corresponding to a grouping provided by the user (generally run after using the `sm.tl.cluster` function)

## Function
"""

# Import
import argparse
import sys
import anndata as ad
#import numba.core
import scanpy as sc
import seaborn as sns; sns.set(color_codes=True)
sns.set_style("white")

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
        help=''
    )
    parser.add_argument(
        '--subsample', type=int, required=False, default='100000',
        help=''
    )
    parser.add_argument(
        '--palette', type=str, required=False, default='tab10',
        help=''
    )
    parser.add_argument(
        '--size', type=int, required=False, default=None,
        help=''
    )
    parser.add_argument(
        '--output_dir', type=str, required=False, default=None,
        help='Path to output directory.'
    )
    args = parser.parse_args(argv[1:])
    print(vars(args))
    cluster_plots(**vars(args))

# Function
def cluster_plots (adata, group_by, subsample=100000, palette ='tab10', 
                   size=None, output_dir=None):
    
    # Load the data 
    if isinstance(adata, str):
        imid = str(adata.rsplit('/', 1)[-1])
        adata = ad.read(adata)
    else:
        adata = adata
        
    # Subset data if needed
    
    
    # UMAP
    sc.pp.neighbors(adata) # Computing the neighborhood graph
    sc.tl.umap(adata)
    fig = sc.pl.umap(adata, color=group_by, palette = palette, size=size, return_fig=True, show=False) # View the clustering
    # save figure
    if output_dir is not None:
        fig.savefig(output_dir + '/umap.pdf')
        
    # Matrix plot
    
    

if __name__ == '__main__':
    main()  

