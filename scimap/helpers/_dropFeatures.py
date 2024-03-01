# -*- coding: utf-8 -*-
# Created on Mon Jul  4 09:27:32 2022
# @author: Ajit Johnson Nirmal
# Drop provided markers or cells from anndata object including the raw matrix


"""
!!! abstract "Short Description"
    `sm.hl.dropFeatures`: This versatile function streamlines the process of 
    refining an AnnData object by enabling users to selectively remove markers, 
    cells, metadata columns, and specific cell groups. It facilitates targeted 
    dataset curation, ensuring analyses are performed on relevant and clean data subsets.

## Function
"""

import numpy as np
import anndata as ad

def dropFeatures (adata, 
                  drop_markers=None, 
                  drop_cells=None, 
                  drop_meta_columns=None,
                  drop_groups=None, 
                  groups_column=None,
                  subset_raw=True,
                  verbose=True):
    """
Parameters:
        adata (anndata.AnnData):  
            Annotated data matrix or path to an AnnData object, containing spatial gene expression data.
        
        drop_markers (list, optional):  
            A list of gene or marker names to be removed from `adata.var`. 
        
        drop_cells (list, optional):  
            A list of cell identifiers (index names) to be removed from `adata.obs`. 
        
        drop_meta_columns (list, optional):  
            A list of metadata column names to be removed from `adata.obs`. 
        
        drop_groups (list, optional):  
            A list of category names to be removed based on the column specified by `groups_column`. 
        
        groups_column (str, optional):  
            The name of the column in `adata.obs` that contains the categorical data for `drop_groups`. 
        
        subset_raw (bool, optional):  
            If True, the same dropping operations are applied to `adata.raw`.
        
        verbose (bool, optional):  
            If True, print messages about the dropping process. 

Returns:
        adata (anndata.AnnData):  
            The AnnData object after the specified features have been removed.

Example:
        ```python
        # Example 1: Drop specific markers from the dataset
        adata = sm.hl.dropFeatures(adata, drop_markers=['CD3D', 'CD19'])
    
        # Example 2: Remove cells based on their identifiers
        adata = sm.hl.dropFeatures(adata, drop_cells=['cell_001', 'cell_002'])
    
        # Example 3: Remove metadata columns from adata.obs
        adata = sm.hl.dropFeatures(adata, drop_meta_columns=['Batch', 'Condition'])
    
        # Example 4: Exclude specific groups from a categorical column in adata.obs
        adata = sm.hl.dropFeatures(adata, drop_groups=['B cell', 'NK cell'], groups_column='Cell_Type')
        
        ```

    """
    
    # Drop Markers
    if drop_markers is not None:
        if isinstance(drop_markers, str):
            drop_markers = [drop_markers]
        # find the index of the given markers
        idx_markers = [adata.var.index.get_loc(x) for x in drop_markers]
        # remove from adata
        keep_markes = list(set(adata.var.index).difference(drop_markers))
        adata = adata[:, keep_markes]
        # remove from raw
        if subset_raw is True:
            raw = np.delete(adata.raw.X, idx_markers, axis=1)
            del adata.raw
            adata.raw = ad.AnnData (raw)
    
    # Drop cells
    if drop_cells is not None:
        if isinstance(drop_cells, str):
            drop_cells = [drop_cells]
        # find the index of the given markers
        idx_markers = [adata.obs.index.get_loc(x) for x in drop_cells]
        # remove from adata
        keep_markes = list(set(adata.obs.index).difference(drop_cells))
        adata = adata[keep_markes, :]
        # remove from raw
        if subset_raw is True:
            raw = np.delete(adata.raw.X, idx_markers, axis=1)
            del adata.raw
            adata.raw = ad.AnnData (raw)
    
    # Drop meta columns
    if drop_meta_columns is not None:
        if isinstance(drop_meta_columns, str):
            drop_meta_columns = [drop_meta_columns]
        # remove from adata
        adata.obs = adata.obs.drop(drop_meta_columns, axis=1)
    
    # Drop specific categories of cells
    if drop_groups is not None:
        if isinstance(drop_groups, str):
            drop_groups = [drop_groups]
        if isinstance(groups_column, list):
            groups_column = groups_column[0]
        # find the index of the given markers
        idx = adata[adata.obs[groups_column].isin(drop_groups)].obs.index
        idx_markers = [adata.obs.index.get_loc(x) for x in idx]
        # remove from raw
        if subset_raw is True:
            raw = np.delete(adata.raw.X, idx_markers, axis=0)
        # remove from adata
        adata = adata[~adata.obs[groups_column].isin(drop_groups)]
        # return adata raw
        if subset_raw is True:
            del adata.raw
            adata.raw = ad.AnnData (raw)

    
    # return
    return adata
    
    
        

        

       
        
        
        

