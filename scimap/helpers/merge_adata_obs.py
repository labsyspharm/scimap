#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
!!! abstract "Short Description"
    `sm.pl.merge_adata_obs`: This function is designed to consolidate multiple AnnData 
    objects originating from the same image or dataset into a single cohesive unit. 
    It is particularly useful when each object contains distinct metadata in `.obs`, 
    such as results from various clustering algorithms executed in parallel. 
    By merging these objects, users can streamline analyses and comparisons within 
    a unified framework. It's important to note that this function is intended for 
    merging objects from the same source and may not be suitable for combining 
    data across different images or datasets.

## Function
"""

import argparse
import sys
import pathlib
import pandas as pd
import anndata as ad


def main(argv=sys.argv):
    parser = argparse.ArgumentParser(
        description='The function allows users to combine multiple anndata objects that are from the same image or dataset but differ in the stored metadata'
    )
    parser.add_argument(
        '--adata', required=True, nargs='*',
        help='AnnData object loaded into memory or path to AnnData object.'
    )
    parser.add_argument(
        '--output_dir', type=str, required=False, default=None,
        help='Path to output directory.'
    )
    args = parser.parse_args(argv[1:])
    print(vars(args))
    merge_adata_obs(**vars(args))
    
    

def merge_adata_obs (adata, 
                     output_dir=None,
                     verbose=True):
    
    
    """
Parameters:
        adata (list of anndata.AnnData or list of str):  
            A list containing AnnData objects to be merged or paths to AnnData files. 
            Each item in the list should either be an AnnData object already loaded into memory 
            or a string representing the path to an AnnData file.

        output_dir (str, optional):  
            The directory where the merged AnnData object should be saved. 
            If specified, the merged object is saved to this directory as 'merged_adata.h5ad'.
            If not specified, the merged AnnData object is not automatically saved to disk.
        
        verbose (bool, optional):  
            If True, prints messages about the renaming process. 

Returns:
        adata (anndata.AnnData):  
            A single AnnData object resulting from the merger of input AnnData objects or files.

Example:
        ```python
        
        # Example 1: Merge AnnData objects already loaded in memory
        combined_adata = sm.hl.merge_adata_obs(adata=[adata1, adata2])
    
        # Example 2: Merge AnnData objects from file paths
        combined_adata = sm.hl.merge_adata_obs(adata=['./data/adata1.h5ad', './data/adata2.h5ad'])
    
        # Example 3: Merge AnnData objects and save the combined object to a specified directory
        combined_adata = sm.hl.merge_adata_obs(adata=[adata1, adata2], output_dir='./merged_data')
        
        ```
    """
    
    #adata = ["/Users/aj/Downloads/mcmicro_output.h5ad", "/Users/aj/Downloads/mcmicro_output_1.h5ad"]
    
    
    # Convert to list of anndata objects
    if isinstance (adata, list):
        adata = adata
    else:
        adata = [adata]
    
    # Resoluve the OBS and UNS section
    if isinstance(adata[0], str):
        df = pd.DataFrame()
        uns_count = []
        for i in adata:
            tmp = ad.read_h5ad(i)
            # OBS
            tmp_df = tmp.obs
            df = pd.concat([df, tmp_df], axis=1)
            # UNS
            uns_count.append(len(tmp.uns))
        # Keep non-duplicate obs columns
        df = df.loc[:,~df.columns.duplicated()]
        # Highest UNS
        uns_index = [i for i, j in enumerate(uns_count) if j == max(uns_count)][0]
    else:
        df = pd.DataFrame()
        uns_count = []
        for i in adata:
            # OBS
            tmp_df = i.obs
            df = pd.concat([df, tmp_df], axis=1)
            # UNS
            uns_count.append(len(i.uns))
        # Keep non-duplicate obs columns
        df = df.loc[:,~df.columns.duplicated()]
        # Highest UNS
        uns_index = [i for i, j in enumerate(uns_count) if j == max(uns_count)][0]
    
    
    # create the final anndata object
    # Load the data 
    if isinstance(adata[0], str):
        final_adata = ad.read_h5ad(adata[uns_index])
    else:
        final_adata = adata[uns_index]
        
    # reindex
    df = df.reindex (final_adata.obs.index)
    
    # replace obs
    final_adata.obs = df
     
    
    # Save data if requested
    if output_dir is not None:
        output_dir = pathlib.Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        final_adata.write(output_dir / 'combined_adata.h5ad')
    else:    
        # Return data
        return final_adata
    
    

if __name__ == '__main__':
    main()  
    
