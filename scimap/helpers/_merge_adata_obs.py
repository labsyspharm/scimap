#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
!!! abstract "Short Description"
    `sm.pl.merge_adata_obs`: The function allows users to combine multiple anndata objects that are from the same image or dataset 
    but differ in the stored metadata `.obs`. For example, multiple clustering methods could have been run in parallel
    leading to generation of multiple anndata objects. This allows to merge them togeather into a single object. Please
    note this cannot be used to merge anndata objects from different images/ datasets. 

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
    
    

def merge_adata_obs (adata, output_dir=None):
    
    
    """
Parameters:
    adata : List of AnnData object loaded into memory or list of path to AnnData object.
        
    output_dir : string, optional  
        Path to output directory.

Returns:
    AnnData :   
        Combined anndata object.
        
Example:
```python
    # combining multiple anndata objects in memory
    adata = sm.hl.merge_adata_obs (adata = [bdata, cdata])
    
    # combining multiple anndata objects by providing their saved loaction
    adata = sm.hl.merge_adata_obs (adata = ['path to bdata.h5ad', 'path to cdata.h5ad'])
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
            tmp = ad.read(i)
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
        final_adata = ad.read(adata[uns_index])
    else:
        final_adata = adata[uns_index]
        
    # reindex
    df = df.reindex (final_adata.obs.index)
    
    # replace obs
    final_adata.obs = df
    
    # Find name of file
    image_path = pathlib.Path(adata[0])    
    
    # Save data if requested
    if output_dir is not None:
        output_dir = pathlib.Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        final_adata.write(output_dir / image_path.name)
    else:    
        # Return data
        return final_adata
    
    

if __name__ == '__main__':
    main()  
    
