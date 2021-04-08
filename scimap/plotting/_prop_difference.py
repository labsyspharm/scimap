#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 17:46:29 2021
@author: Ajit Johnson Nirmal
Fold Change between cell types
"""


# lib
import seaborn as sns; sns.set(color_codes=True)
import matplotlib
import numpy as np
sns.set_style("white")


# Function

def prop_difference (adata, label='prop_difference', p_val=0.05, nonsig_color='grey',subset=None,
                     cmap = 'vlag', log=True,center=0, **kwargs):
    
    
    
    # set color for heatmap
    #cmap_updated = copy.copy(matplotlib.cm.get_cmap(cmap))
    cmap_updated = matplotlib.cm.get_cmap(cmap)
    cmap_updated.set_bad(color=nonsig_color)
    
    
    # get the data
    fc = adata.uns[str(label)+'_foldchange']
    p = adata.uns[str(label)+'_pval']
    
    #fold
    fold = fc.copy()
    
    # log
    if log is True:
        fold = np.log2(fold)
              
    
    # create a mask for non-sig values
    p_mask = p.copy()
    p_mask[p_mask > p_val] = np.nan
    
    # subset data
    if subset is not None:
        if isinstance (subset, str):
            subset = [subset]
        fold = fold [subset]
        p_mask = p_mask [subset]
    
    #mask
    mask = p_mask.isnull() # identify the NAN's for masking 
    
    # heatmap of the foldchange
    sns.clustermap(fold, cmap=cmap, mask=mask, center=center, **kwargs)
    

