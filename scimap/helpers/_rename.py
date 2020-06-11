#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 13:08:26 2020
@author: Ajit Johnson Nirmal
Rename any string within a column to anything
"""

import functools 
import re


def rename (adata, rename, from_column='phenotype', to_column='phenotype_renamed'):
    """
    
    Parameters
    ----------
    adata : AnnData object
    rename : dict
        Pass a dictionary with 'values' as elements that need to be altered and 'keys' as the elements that they need to be transformed into.
    from_column : string, required
        Column that need to be modified. The default is 'phenotype'.
    to_column : string, required
        Modified names will be stored in a new column with this name. The default is 'phenotype_renamed'.

    Returns
    -------
    adata : Modified AnnData Object
        DESCRIPTION.
    
    Example
    -------
    rename= {'tumor': ['cd45 neg tumor', 'cd8 tumor', 'cd4 tumor'],
             'macrophages': ['m1 macrophages', 'm2 macrophages']}
    Here we are renaming cd45 neg tumor, cd8 tumor and cd4 tumor into 'tumor' and 
    m1 macrophages and m2 macrophages into macrophages
    
    adata = sm.hl.rename_clusters (adata, rename, from_column='phenotype', to_column='phenotype_renamed')

    """
    
    # Get the from_column
    rename_from = list(adata.obs[from_column].values)
    
    # Split multiple renaming events into independent events
    name = functools.reduce( lambda x,y: dict(x, **y), (dict(map(lambda x: (x,i), rename[i])) for i in rename))
     
    # Rename
    for i in name:
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
