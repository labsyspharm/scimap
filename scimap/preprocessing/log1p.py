# -*- coding: utf-8 -*-
#Created on Thu Nov 10 09:53:28 2022
#@author: Ajit Johnson Nirmal
#Function to incorporate log transformation

"""
!!! abstract "Short Description"
    `sm.pp.log1p` applies a log1p transformation to the raw data of an AnnData object and saves the transformed data in a specified layer. 
    If a string path to an AnnData file is provided, the file is loaded, transformed, and  overwritten with the 
    transformed data. The function ensures the raw data exists before applying the transformation and allows for verbose output.

## Function
"""

# import library
import anndata as ad
import numpy as np
from pathlib import Path
import argparse


def log1p (adata,
           layer='log',  
           verbose=True):
    
    """
Parameters:
    adata (str or anndata.AnnData): 
        The AnnData object or a string path to an AnnData file. If a path is provided, the function will load the AnnData object from the file. The transformation is applied to the raw data of this object.
        
    layer (str): 
        Name of the new layer where the log-transformed data will be stored. Default is 'log'. If the layer already exists, it will be overwritten with the new log-transformed data.
        
    verbose (bool): 
        If True, the function will print messages about its progress and actions, including warnings if the specified layer already exists and confirmation when the AnnData object is saved to a file. Default is True.

    
Returns:
    adata (anndata): 
        The modified AnnData object, only if the input is an AnnData object and not a file path. If a file path is provided and successfully processed, the function returns None.
    
Example:
    ```python
    
    # Using a path to an AnnData file:
    # will overwrite the original file with the transformed AnnData object.
    
    sm.pp.log1p("path/to/your/adata.h5ad", layer="log", verbose=True)
    
    # Using an AnnData object directly:
    # Assuming `adata` is your pre-loaded into memory
    # This will apply log1p transformation and the modified AnnData object
    
    adata = sm.pp.log1p(adata, layer="log", verbose=True)
        
    ```
"""
    
    # load adata
    if isinstance(adata, str):
        adata_path = Path(adata)
        if not adata_path.exists():
            raise FileNotFoundError(f"The file {adata} does not exist.")
        adata = ad.read_h5ad(adata_path)
    else:
        adata_path = None
    
    
    if layer in adata.layers:
        if verbose:
            print(f"Warning: Layer '{layer}' already exists. It will be overwritten with the new log-transformed data.")
    
    if adata.raw is None:
        raise AttributeError("adata.raw does not exist. Please assign RAW data to adata.raw before proceeding (e.g., adata.raw = adata).")
        
    # perform the operation
    adata.layers[layer] = np.log1p(adata.raw.X)
    
    # return 
    # Overwrite the original file with the modified AnnData object if requested
    if isinstance(adata_path, Path):
        adata.write_h5ad(adata_path)
        if verbose: 
            print(f"Modified AnnData object has been saved to {adata_path}")
    else:
        return adata
        
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Apply a log1p transformation to the raw data of an AnnData object and save the transformed data in a specified layer.')

    parser.add_argument('--adata', type=str, required=True, help='Path to an AnnData object file. The AnnData object will be loaded from this file.')
    parser.add_argument('--layer', type=str, default='log', help="Name of the new layer where the log-transformed data will be stored. Default is 'log'.")
    parser.add_argument('--verbose', action='store_true', help="If set, print messages about progress and actions. Default is False.")
    
    args = parser.parse_args()

    # Call log1p function with the provided arguments
    log1p(args.adata, layer=args.layer, verbose=args.verbose)