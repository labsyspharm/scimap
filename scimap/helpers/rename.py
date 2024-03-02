#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Sun Mar 22 13:08:26 2020
# @author: Ajit Johnson Nirmal
"""
!!! abstract "Short Description"
    `sm.hl.rename`: This function offers a straightforward way to rename specific 
    categories within a chosen column of an AnnData object, with the new names 
    being stored in a separate column. It streamlines the process of updating or 
    consolidating category labels for enhanced data clarity and analysis.

## Function
"""
# Import
import functools 
import re

# Function
def rename (adata, 
            rename, 
            from_column='phenotype', 
            to_column='phenotype_renamed',
            verbose=True):
    """
Parameters:
        adata (anndata.AnnData):  
            Annotated data matrix or path to an AnnData object, containing spatial gene expression data.

        rename (dict):  
            A dictionary mapping existing category names (values) to new category names (keys). 
            Each key corresponds to the new name, and its value is a list of existing names to be consolidated under this new name.

        from_column (str, optional):  
            The name of the column in `adata.obs` where the categories to be renamed are located. Defaults to 'phenotype'.

        to_column (str, optional):  
            The name of the new column in `adata.obs` where the renamed categories will be stored. Defaults to 'phenotype_renamed'.

        verbose (bool, optional):  
            If True, prints messages about the renaming process. 

Returns:
        adata (anndata.AnnData):  
            The AnnData object after applying the renaming operation, with the newly named categories stored in the specified `adata.obs[to_column]`.

Example:
    ```python
    
    # Example 1: Simplify phenotype labels
    rename_dict = {'tumor': ['cd45 neg tumor', 'cd8 tumor', 'cd4 tumor'],
                   'macrophages': ['m1 macrophages', 'm2 macrophages']}
    adata = sm.hl.rename(adata, rename=rename_dict, from_column='phenotype', to_column='simplified_phenotype')

    # Example 2: Merge similar phenotypes under a common name
    merge_dict = {'immune cells': ['cd45+', 't-cells', 'b-cells']}
    adata = sm.hl.rename(adata, rename=merge_dict, from_column='cell_type', to_column='merged_cell_type')

    # Example 3: Rename and create a new column for easier identification
    new_names = {'activated': ['activated_tcells', 'activated_bcells'],
                 'resting': ['resting_tcells', 'resting_bcells']}
    adata = sm.hl.rename(adata, rename=new_names, from_column='status', to_column='status_simplified')
    
    ```
    """
    
    # Sanity check: if the values are not list convert them into list
    for i in rename:
        if isinstance(rename[i], str):
            rename[i] = [rename[i]]
    
    # Get the from_column
    rename_from = list(adata.obs[from_column].values)
    
    # Split multiple renaming events into independent events
    name = functools.reduce( lambda x,y: dict(x, **y), (dict(map(lambda x: (x,i), rename[i])) for i in rename))
     
    # Rename
    for i in name:
        if verbose:
            print ('Renaming ' + str(i) + ' to ' + str(name[i]))
        #rename_from = [x.replace(i, name[i]) for x in rename_from]
        s = str(i)
        s = s.replace('+', '\+')
        #rename_from = [re.sub(r'^\b%s\b$' % s,  name[i], j) for j in rename_from]
        rename_from = [re.sub(r'^\b%s$' % s,  name[i], j) for j in rename_from]
        
    
    # Append to adata as a new column
    adata.obs[to_column] = rename_from
    
    # Return
    return adata
