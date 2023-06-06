# -*- coding: utf-8 -*-
# Created on Mon Jul  4 09:27:32 2022
# @author: Ajit Johnson Nirmal
# Drop provided markers or cells from anndata object including the raw matrix


"""
!!! abstract "Short Description"
    `sm.hl.dropFeatures`:  A handy function that allows users to susbet the anndata object.  
    
    The function can be used to drop markers, cells, columns in metadata, groups of cells 
    belonging to a certain category. 

## Function
"""

import numpy as np
import anndata as ad

def dropFeatures (adata, drop_markers=None, drop_cells=None, 
                  drop_meta_columns=None,
                  drop_groups=None, groups_column=None,
                  subset_raw=True):
    """
Parameters  

    adata : AnnData Object  
    
    drop_markers (list):  
        Provide a list of markers to drop. The default is None.
        
    drop_cells (list):  
        Provide a list of cells (index name) to drop. The default is None.
        
    drop_meta_columns (list):  
        Provide a list of column names in `adata.obs` to drop. The default is None.
        
    drop_groups (list):  
        Provide a list of categorical names to drop. 
        Works in conjunction with `groups_column`. The default is None.
        
    groups_column (str):  
        Pass the column name of the column that contains the categories passed to 
        `drop_groups`. The default is None.
    
    subset_raw (bool):  
        Generally any subsetting of `AnnData` object does not affect the raw data 
        stored in `adata.raw`. Pass `True` to apply the same transformations to 
        `adata.raw` as well. The default is True.

Returns  

    adata : Modified `adata` object

Example

```python

# Drop certain genes 
drop_markers = ['ELANE', 'CD57']
adata = sm.hl.dropFeatures (adata, drop_markers=drop_markers)

# Drop a few cells
drop_cells = ['unmicst-exemplar-001_cell_1', 'unmicst-exemplar-001_cell_2','unmicst-exemplar-001_cell_3']
adata = sm.hl.dropFeatures (adata, drop_cells=drop_cells)

# Drop a few columns from `adata.obs`
drop_meta_columns = ['ROI', 'shapes']
adata = sm.hl.dropFeatures (adata, drop_meta_columns=drop_meta_columns)

# Drop two cell-types from the 'phenotype' column
drop_groups = ['Unknown', 'Treg']
adata = sm.hl.dropFeatures (adata, drop_groups=drop_groups, groups_column = 'phenotype')

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
    
    
        

        

       
        
        
        

