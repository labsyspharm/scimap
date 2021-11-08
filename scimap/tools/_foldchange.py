#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Wed Apr  7 13:22:29 2021
# @author: Ajit Johnson Nirmal
"""
!!! abstract "Short Description"
    `sm.tl.foldchange`: The function allows users to compute the foldchange (fc) in cell-type (phenotype) abundance 
    between samples or ROI's. 

    The reference sample or ROI needs to be passed via the `from_group` parameter. 
    The column name of `from_group` should be passed via `imageid`. The function computes the fc 
    to all other categories within the same `imageid` column. By default (can be turned off), the cell-abundance will
    be normalized for the total number of cells within the sample/ROI to account for difference in area.
    A `fisher-exact-test` is performed to compute the p-values.

    The results are stored in `.uns` section of the anndata object. 

## Function
"""

# library
import scipy.stats as stats
import pandas as pd
import numpy as np

# Function
def foldchange (adata, from_group, to_group=None, imageid='imageid', phenotype='phenotype',
                normalize=True, subset_phenotype=None, label='foldchange'):
    """
    

Parameters:
    adata : AnnData object

    from_group : list, required  
        Pass in the name of the sample or ROI that will serve as a reference for calculating fold change.
        If multiple sample names or ROI's are passed in as a list e.g. ['ROI1', 'ROI2''], please note that
        they will be combined for calculating the fold change. 

    to_group : list, optional  
        By default the reference sample/ROI passed via `from_group` will be compared to all other groups
        within the same column. However if users wish to restrict the comparision to a subset of
        samples/ROI's they can be passed though this paramenter as a list. e.g. ['ROI3', 'ROI4']. 

    imageid : string, optional  
        The column that contains the samples/ROI information.

    phenotype : string, optional  
        The column that contains the cell-type/phenotype information.

    normalize : bool, optional  
        Inorder to account for the sample/ROI area, the cellular abundance is first normalized
        to the total number of cells within the respective sample/ROI. Please note if you pass values in
        `subset_phenotype`, the abundance normalization is restricted to the total cells of the 
        cell types passed in via `subset_phenotype`.

    subset_phenotype : list, optional  
        If users are interested in only a subset of cell-types, the names of those can be passed in through
        this parameter. The data is subsetted to include only these cell types before computing foldchange.

    label : string, optional  
        Key for the returned data, stored in `adata.uns`. The foldchange and p-values 
        are returned seperately with the postfix `_fc` and `_pval`. 

Returns:
    adata : Updated anndata object
        Check `adata.uns['foldchange_fc']` and `adata.uns['foldchange_pval']` for results.
    
    
Example:
    ```python
    adata = sm.tl.foldchange (adata, from_group='image_1', to_group=None, 
                              imageid='imageid', phenotype='phenotype',
                              normalize=True, 
                              subset_phenotype=['Tcells','Bcells','Macs'], 
                              label='foldchange')
    
    ```


    """
    
    # prepare data
    data = adata.obs[[imageid,phenotype]]
    
    # convert from and to groups to list
    if isinstance(from_group, str):
        from_group = [from_group]
    if isinstance(to_group, str):
        to_group = [to_group]
        
    # subset phenotype of interest
    if subset_phenotype is not None:
        if isinstance (subset_phenotype, str):
            subset_phenotype = [subset_phenotype]
        data = data[data[phenotype].isin(subset_phenotype)]
    
    # subset data    
    from_data = data[data[imageid].isin(from_group)]
    from_data[imageid] = from_data[imageid].astype('str').astype('category')
    from_data[phenotype] = from_data[phenotype].astype('str').astype('category')
    if to_group is None:
        to_data = data[~data[imageid].isin(from_group)]
    else:
        to_data = data[data[imageid].isin(to_group)]
    to_data[imageid] = to_data[imageid].astype('str').astype('category')
    to_data[phenotype] = to_data[phenotype].astype('str').astype('category')
    
    # consolidated counts dataframe
    from_data_consolidated = pd.DataFrame(from_data.groupby([imageid,phenotype]).size()).unstack().fillna(0)
    from_data_consolidated.columns = np.unique(from_data_consolidated.columns.get_level_values(1))
    
    to_data_consolidated = pd.DataFrame(to_data.groupby([imageid,phenotype]).size()).unstack().fillna(0)
    to_data_consolidated.columns = np.unique(to_data_consolidated.columns.get_level_values(1))
    
    # make backup of the sample names
    from_b = list(from_data_consolidated.index)
    to_b = list(to_data_consolidated.index)
    
    # make sure from_data_consolidated and to_data_consolidated has the same columns
    x = from_data_consolidated.T
    x.columns = x.columns.astype(str)
    y = to_data_consolidated.T
    y.columns = y.columns.astype(str)
    consolidated = x.merge(y, how='outer', left_index=True, right_index=True).fillna(0)
    
    # split it back into from and to
    from_data_consolidated = consolidated[from_b].T
    to_data_consolidated = consolidated[to_b].T
    
    # create the total minus to and from tables
    from_data_total = abs(from_data_consolidated.sub( from_data_consolidated.sum(axis=1), axis=0))
    to_data_total = abs(to_data_consolidated.sub( to_data_consolidated.sum(axis=1), axis=0))
    
    # 
    print('calculating P values')
    p_vals = []
    for i in from_data_consolidated.columns:
        a = from_data_consolidated[i][0]
        c = from_data_total[i][0]
        for j in to_data_consolidated.index:
            b = to_data_consolidated[i][j]
            d = to_data_total[i][j]
            oddsratio, pvalue = stats.fisher_exact([[a, b], [c, d]])
            p_vals.append(pvalue)
    
    # replace 0 with a small number (1 cell) to avoind inf
    from_data_consolidated_zero = from_data_consolidated.replace(0, 1, inplace=False)
    to_data_consolidated_zero = to_data_consolidated.replace(0, 1, inplace=False)
    
    # normalize based on area i.e total cells if user requests
    if normalize is True:
        # Normalize for total cells
        from_data_ratio = from_data_consolidated_zero.div(from_data_consolidated_zero.sum(axis=1), axis=0)
        to_data_ratio = to_data_consolidated_zero.div(to_data_consolidated_zero.sum(axis=1), axis=0)   
    else:
        from_data_ratio = from_data_consolidated_zero
        to_data_ratio = to_data_consolidated_zero
        
    # foldchange
    fold_change = to_data_ratio.div(from_data_ratio.values,  axis=1)
    fold_change.index.name = '-'.join(from_group)
    
    # reshape the pvalues to the todata df
    p_values = np.reshape(p_vals, to_data_consolidated.shape)
    p_values = pd.DataFrame(p_values, columns = to_data_consolidated.columns, index= to_data_consolidated.index)
    
    # return data
    adata.uns[str(label)+'_pval'] = p_values
    adata.uns[str(label)+'_fc'] = fold_change
    
    return adata

    