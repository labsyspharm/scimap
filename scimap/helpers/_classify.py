#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Mon Oct 26 12:04:17 2020
# @author: Ajit Johnson Nirmal
"""
!!! abstract "Short Description"
    `sm.hl.classify`: Helper function that allow users to annotate cells based on positivity/negativity 
    of defined markers. Users can classify the entire data or a subset of data that 
    has been previously phenotyped or clustered.

## Function
"""


#Library
import pandas as pd
import numpy as np


# Functions
def classify (adata, pos=None, neg=None, classify_label='passed_classify', failed_label='failed_classify',
              phenotype=None,subclassify_phenotype=None,threshold = 0.5,
              collapse_failed=True,label="classify",showPhenotypeLabel=True):
    
    """
Parameters:

    adata : AnnData object

    pos : list, optional  
        Pass a list of markers that should be expressed in the resultant cells.

    neg : list, optional  
        Pass a list of markers that should not be expressed in the resultant cells.

    classify_label : string, optional  
        Provide a name for the classified cells.

    failed_label : string, optional
        Provide a name for cells that did not pass classify.

    subclassify_phenotype : list, optional  
        If only a subset of phenotypes require to classified, pass the name of those phenotypes as a list
        through this argument.

    threshold: float, optional  
        Above or below the given value will be considered for positive and negative classification.
        If the data was scaled using the `sm.pp.rescale` function, 0.5 is the classification threshold.

    phenotype : string, required  
        Column name of the column containing the phenotype information. 
        This is important if `subclassify_phenotype` or `collapse_failed` arguments are used.

    collapse_failed : bool, optional  
        If set to true, the cells that were not classified based on the given criteria will be
        binned into a single category named 'failed_classify'. When False, the phenotype
        information for other cells will be borrowed from the `phenotype` argument.
        
    showPhenotypeLabel : bool, optional
        If set to True, returns the data under [phenotype]_[label] key,
        stored in `adata.obs`. Each cell's classification status will be appended to its phenotype,
        [phenotype]_[classified_label] or [phenotype]_[failed_label].

    label : string, optional  
        Key for the returned data, stored in `adata.obs`.

 Returns:

    adata : AnnData  
        Updated AnnData Object.
    
    
Example:
```python
    # Classify all cells with both pos and neg markers 
    # (Identify cytotoxic T-cells)
    adata = sm.hl.classify(adata, pos=['CD3D','CD8A'], neg=['ASMA'])
    
    # Classify specific sub-types of cells
    adata = sm.hl.classify(adata, pos=['CD3D','FOXP3'], 
    neg=['ASMA'], subclassify_phenotype=['T cells','Regulatory T cells'])
    
    # Classify specific sub-types of cells and borrow labels 
    # from another column
    adata = sm.hl.classify(adata, pos=['CD3D'], neg=['ASMA'], 
    subclassify_phenotype=['T cells'], collapse_failed=False, 
    phenotype='phenotype')
```
    """
    
    # clean the input
    if isinstance(pos, str):
        pos = [pos]
    if isinstance(neg, str):
        neg = [neg]
    if isinstance(subclassify_phenotype, str):
        subclassify_phenotype = [subclassify_phenotype]
    if (showPhenotypeLabel):
        phenotype_label=phenotype+"_"+label
    
    
    # Create a dataFrame with the necessary inforamtion
    data = pd.DataFrame(adata.X, index= adata.obs.index, columns = adata.var.index)
    
    # if user requests to subset a specific phenotype   
    if subclassify_phenotype is not None:
        meta = pd.DataFrame(adata.obs[phenotype])
        subset_index = meta[meta[phenotype].isin(subclassify_phenotype)].index
        data = data.loc[subset_index]
        
    # Subset cells that pass the pos criteria
    if pos is not None:
        for i in pos:
            data = data[data[i] >= threshold]
                
    # Subset cells that pass the neg criteria 
    if neg is not None and not data.empty:
        for j in neg:
            data = data[data[j] < threshold]
    
    # Cells that passed the classify criteria
    # Create classified, indexed by data.index, with label, phenotype_label columns
    if data.empty:
        raise TypeError("No cells were found to satisfy your `classify` criteria")
    else:
        # create new naming scheme for label and phenotype_label cols in classified
        non_summary = pd.DataFrame({phenotype: adata.obs[phenotype]}) # gets the index and phenotype
        non_summary[phenotype] = non_summary[phenotype].astype(str)

        classify_idx=data.index
        classified = pd.DataFrame(non_summary.loc[data.index]) #subsets phenotype rows to only classified cells
        if showPhenotypeLabel:
            classified[phenotype_label] = classified[phenotype]+"_"+classify_label # add phenotype_label col
        classified[label]=pd.DataFrame(np.repeat(classify_label, len(classify_idx)), index = classify_idx) # add label col
        classified.drop([phenotype], axis='columns', inplace=True) # drop phenotype col, for merge        



    if collapse_failed is True: 
        meta = non_summary # has index and phenotype col
        meta = meta.merge(classified, how='outer', left_index=True, right_index=True) # gain classified col(s) and NaNs for non-matches
        if showPhenotypeLabel is True:
            meta[phenotype_label]= meta[phenotype_label].fillna(meta[phenotype].astype(str)+"_"+failed_label)
            meta=meta[phenotype_label]
        else: 
            meta[label]=meta[label].fillna(failed_label)
            meta=meta[label]
    
        
    else:
        if phenotype is None:
            raise ValueError("Please pass a column name to the PHENOTYPE argument")
        
        if showPhenotypeLabel is True: 
            meta=non_summary # phenotype col
            classified=pd.DataFrame({phenotype: classified[phenotype_label]}) # takes phenotype_label col and renames to phenotype, ensures it's a df
            meta.update(classified) # updates with phenotype_label for only the classified cells
        else:
            meta= pd.DataFrame(adata.obs[phenotype])
            classified = pd.DataFrame(np.repeat(classify_label, len(classify_idx)), index = classify_idx, columns = [phenotype])
            meta.update(classified) # updates with label for only the classified cells
        
            
    # Add to Anndata 
    meta = meta.reindex(adata.obs.index)
    if showPhenotypeLabel is True:
        adata.obs[phenotype_label]=meta
    else:
        adata.obs[label]=meta 
            
    # return
    return adata

