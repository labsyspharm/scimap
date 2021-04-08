#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 13:22:29 2021
@author: Ajit Johnson Nirmal
Fisher-Exact test to test for significant changes in proportion of cell-types/ categories of interest
"""

# lib
import scipy.stats as stats
import pandas as pd
import numpy as np


from_group = 'Ton_192_n'

# Function
def prop_difference (adata, from_group, to_group=None, imageid='imageid', subset=None,
                     phenotype='phenotype',label='prop_difference'):
    
    # prepare data
    data = adata.obs[[imageid,phenotype]]
    
    # convert from and to groups to list
    if isinstance(from_group, str):
        from_group = [from_group]
    if isinstance(to_group, str):
        to_group = [to_group]
        
    # subset phenotype of interest
    if subset is not None:
        if isinstance (subset, str):
            subset = [subset]
        data = data[data[phenotype].isin(subset)]
    
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
    
    # replace 0 with a small number to avoind inf
    from_data_consolidated_zero = from_data_consolidated.replace(0, 1, inplace=False)
    to_data_consolidated_zero = to_data_consolidated.replace(0, 1, inplace=False)
    # Normalize for total cells
    from_data_ratio = from_data_consolidated_zero.div(from_data_consolidated_zero.sum(axis=1), axis=0)
    to_data_ratio = to_data_consolidated_zero.div(to_data_consolidated_zero.sum(axis=1), axis=0)    
    # foldchange
    fold_change = to_data_ratio.div(from_data_ratio.values,  axis=1)
    
    # reshape the pvalues to the todata df
    p_values = np.reshape(p_vals, to_data_consolidated.shape)
    p_values = pd.DataFrame(p_values, columns = to_data_consolidated.columns, index= to_data_consolidated.index)
    
    # return data
    adata.uns[str(label)+'_pval'] = p_values
    adata.uns[str(label)+'_foldchange'] = fold_change
    
    return adata

    

    

