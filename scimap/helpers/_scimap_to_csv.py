# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 21:00:57 2020
@author: Ajit Johnson Nirmal
Helper function to save andata object as a csv file
"""

# Lib
import pandas as pd
import numpy as np


# Function

def scimap_to_csv (adata, data_type='raw'):
    """
    

    Parameters
    ----------
    adata : AnnData Object

    data_type : string, optional
        Three options are available:
        1) 'raw' - The raw data will be returned.
        2) 'log' - The raw data converted to log scale using `np.log1p` will be returned.
        3) 'scaled' - If you have scaled the data using the `sm.pp.rescale`, that will be
        returned. Please note, if you have not scaled the data, whatever is within
        `adata.X` will be returned.
        The default is 'raw'.

    Returns
    -------
    merged : DataFrame
        A single dataframe containing the expression and metadata will be returned.
        
    Example
    -------
    data = sm.hl.scimap_to_csv (adata, data_type='raw')

    """
    
    # Expression matrix
    if data_type is 'raw':
        data = pd.DataFrame(adata.raw.X, index=adata.obs.index, columns=adata.var.index)
    if data_type is 'log':
        data = pd.DataFrame(np.log1p(adata.raw.X), index=adata.obs.index, columns=adata.var.index)
    if data_type is 'scaled':
        data = pd.DataFrame(adata.X, index=adata.obs.index, columns=adata.var.index)
    
    # Metadata
    meta = pd.DataFrame(adata.obs)
    
    # Merge the two dataframes
    merged = pd.concat([data, meta], axis=1, sort=False)
    
    # Add a column to save cell-id
    merged['CellID'] = merged.index
    
    # reset index
    merged = merged.reset_index(drop=True)
    
    # return
    return merged

