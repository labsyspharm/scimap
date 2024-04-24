#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Wed Apr  7 13:22:29 2021
# @author: Ajit Johnson Nirmal
"""
!!! abstract "Short Description"
    The `sm.tl.foldchange` function computes the fold change in cell-type abundance between samples or ROIs, 
    using the `from_group` parameter to specify the reference group by its column name in `imageid`. 
    It normalizes cell abundance to the total cell count within each sample/ROI to adjust for size differences, 
    a feature that can be disabled. The function uses a Fisher exact test to calculate p-values, assessing the 
    statistical significance of the observed changes.  
    
    Results are stored in the `.uns` section of the Anndata object for easy access and further analysis.

## Function
"""

# library
import scipy.stats as stats
import pandas as pd
import numpy as np
import argparse

# Function
def foldchange (adata, 
                from_group, 
                to_group=None, 
                imageid='imageid', 
                phenotype='phenotype',
                normalize=True, 
                subset_phenotype=None, 
                label='foldchange',
                verbose=True):
    """
    

Parameters:
    adata (anndata.AnnData):
        The input AnnData object containing single-cell data for fold change analysis.

    from_group (list of str):  
        Specifies the reference sample(s) or Region of Interest (ROI) for calculating fold change. If multiple samples or ROIs are provided (e.g., ['ROI1', 'ROI2']), they will be aggregated to serve as a singular reference point for comparison.

    to_group (list of str, optional):  
        Defines a specific set of samples/ROIs to compare against the reference group specified in `from_group`. If not provided, the reference will be compared to all other groups within the `imageid` column. For example, ['ROI3', 'ROI4'].

    imageid (str):  
        The column in `adata.obs` that holds the sample/ROI identifiers.

    phenotype (str):  
        The column in `adata.obs` that contains cell type or phenotype information.

    normalize (bool):  
        If True, adjusts cell counts based on the total number of cells within each sample/ROI to account for differences in sample/ROI area. If `subset_phenotype` is provided, normalization considers only the total cells of the specified cell types.

    subset_phenotype (list of str, optional):  
        Limits the analysis to a particular subset of cell types. Only cell types listed here will be included in the fold change computation.

    label (str):   
        Designates the key under which the fold change results (both fold change values and p-values) are stored in `adata.uns`. The results will be accessible as `<label>_fc` for fold changes and `<label>_pval` for p-values.

    verbose (bool):  
        Enables the display of detailed progress updates and information about the execution steps when set to True.

Returns:
    adata (anndata.AnnData):
        The input `adata` object, updated with fold change analysis results. The fold change values and p-values can be found in `adata.uns['<label>_fc']` and `adata.uns['<label>_pval']`, respectively.
    

Example:
    ```python
    
    # Basic usage with automatic comparison to all other groups
    adata = sm.tl.foldchange(adata, from_group=['ROI1'], imageid='imageid', phenotype='phenotype', normalize=True, label='roi_comparison')
    
    # Specifying a subset of groups for comparison
    adata = sm.tl.foldchange(adata, from_group=['image_1'], to_group=['image_2', 'image_3'], imageid='imageid', phenotype='phenotype', normalize=True, label='specific_roi_comparison')
    
    # Focusing on specific cell types for fold change analysis
    adata = sm.tl.foldchange(adata, from_group=['ROI1'], to_group=['ROI3', 'ROI4'], subset_phenotype=['T cells', 'B cells', 'Macrophages'], label='subset_phenotype_comparison')
    
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
    if len(from_group) > 1:
        combined_name = '_'.join(from_group)
        #from_data[imageid] = combined_name
        from_data.loc[:, imageid] = combined_name


    from_data.loc[:, imageid] = from_data[imageid].astype('str').astype('category')
    from_data.loc[:, phenotype] = from_data[phenotype].astype('str').astype('category')
    if to_group is None:
        to_data = data[~data[imageid].isin(from_group)]
    else:
        to_data = data[data[imageid].isin(to_group)]
    to_data.loc[:, imageid] = to_data[imageid].astype('str').astype('category')
    to_data.loc[:, phenotype] = to_data[phenotype].astype('str').astype('category')
    

    if verbose:
        print('calculating foldchange')
    
    # consolidated counts dataframe
    from_data_consolidated = pd.DataFrame(from_data.groupby([imageid,phenotype],observed=False).size()).unstack().fillna(0)
    from_data_consolidated.columns = np.unique(from_data_consolidated.columns.get_level_values(1))
    
    to_data_consolidated = pd.DataFrame(to_data.groupby([imageid,phenotype],observed=False).size()).unstack().fillna(0)
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
    if verbose:
        print('calculating P values')
    p_vals = []
    for i in from_data_consolidated.columns:
        #a = from_data_consolidated[i][0]
        a = from_data_consolidated[i].iloc[0]
        #c = from_data_total[i][0]
        c = from_data_total[i].iloc[0]
        for j in to_data_consolidated.index:
            #b = to_data_consolidated[i][j]
            b = to_data_consolidated[i].loc[j]
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

    # Make the Function CLI compatable
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='The function allows users to compute the foldchange (fc) in cell-type (phenotype) abundance between samples or ROIs.')
    parser.add_argument('--adata', type=str, help='AnnData object')
    parser.add_argument('--fromgroup', type=str, help='Pass in the name of the sample or ROI that will serve as a reference for calculating fold change')
    parser.add_argument('--togroup', type=str, default=None, help='By default the reference sample/ROI passed via from_group will be compared to all other groupswithin the same column')
    parser.add_argument('--imageid', type=str, default='imageid', help='The column that contains the samples/ROI information.')
    parser.add_argument('--phenotype', type=str, default='phenotype', help='The column that contains the cell-type/phenotype information.')
    parser.add_argument('--normalize', type=bool, default=True, help='Inorder to account for the sample/ROI area, the cellular abundance is first normalized to the total number of cells within the respective sample/ROI')
    parser.add_argument('--subsetphenotype', type=str, default=None, help='If users are interested in only a subset of cell-types, the names of those can be passed in through this parameter')
    parser.add_argument('--label', type=str, default='foldchange', help='Key for the returned data, stored in `adata.uns`. The foldchange and p-values')
    parser.add_argument('--verbose', required=False, default=True, help='The function will print detailed messages about its progress.')
    args = parser.parse_args()
    
    foldchange(adata=args.adata,
                   from_group=args.fromgroup, 
                   to_group=args.togroup, 
                   imageid=args.imageid, 
                   phenotype=args.phenotype, 
                   normalize=args.normalize,
                   subset_phenotype=args.subsetphenotype,
                   label=args.label)
    
   