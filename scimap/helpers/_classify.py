#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Mon Oct 26 12:04:17 2020
# @author: Ajit Johnson Nirmal
"""
!!! abstract "Short Description"
    `sm.hl.classify`: This utility function enables users to annotate cells by assessing 
    the presence or absence of specific markers. It offers flexibility to apply classifications 
    across the entire dataset or within previously defined subsets, such as phenotyped or 
    clustered cell groups, facilitating targeted analyses based on marker expression.

## Function
"""


#Library
import pandas as pd
import numpy as np


# Functions
def classify (adata, 
              pos=None, 
              neg=None, 
              classify_label='passed_classify', 
              failed_label='failed_classify',
              phenotype=None,
              subclassify_phenotype=None,
              threshold = 0.5,
              collapse_failed=True,
              label="classify",
              showPhenotypeLabel=False,
              verbose=True):
    
    """
Parameters:
        adata (anndata.AnnData):  
            The annotated data matrix for classification.
            
        pos (list, optional):  
            Markers that should be expressed in the cells of interest.
            
        neg (list, optional):  
            Markers that should not be expressed in the cells of interest.
            
        classify_label (str, optional):  
            Label for cells that meet the classification criteria.
            
        failed_label (str, optional):  
            Label for cells that do not meet the classification criteria.
            
        phenotype (str, required if subclassify_phenotype or collapse_failed is used):  
            Column in `adata.obs` containing the phenotype information.
            
        subclassify_phenotype (list, optional):  
            Phenotypes within which classification should be performed.
            
        threshold (float, optional):  
            Threshold for determining positive or negative expression.
            
        collapse_failed (bool, optional):  
            If True, unclassified cells are grouped under a single failed label.
            
        label (str, optional):  
            Key under which classification results are stored in `adata.obs`.
            
        showPhenotypeLabel (bool, optional):  
            If True, appends classification status to existing phenotype labels in the results. If True, classification
              results will instead be stored under "[phenotype]_[label]" key in  `adata.obs`
            
        verbose (bool, optional):  
            If True, prints progress and informational messages during the classification process.

Returns:
        adata (anndata.AnnData):  
            The input AnnData object, updated with classification results in `adata.obs[label]`.

Example:
    ```python
    
    # Basic classification with positive and negative markers
    adata = sm.hl.classify(adata, pos=['CD3D', 'CD8A'], neg=['PDGFRB'], label='T_cell_classification')
    
    # Classify specific phenotypes, preserving original phenotype labels for unclassified cells
    adata = sm.hl.classify(adata, pos=['CD19'], neg=['CD3D'], subclassify_phenotype=['B cells'],
                     phenotype='cell_type', collapse_failed=False, label='B_cell_subclassification')

    # Use showPhenotypeLabel to append classification status to existing phenotype labels
    adata = sm.hl.classify(adata, pos=['CD34'], neg=['CD45'], phenotype='cell_type',
                     showPhenotypeLabel=True, label='stem_cell_classification', verbose=True)
    
    ```
    """
    # clean the input
    if isinstance(pos, str):
        pos = [pos]
    if isinstance(neg, str):
        neg = [neg]
    if phenotype is not None:
        if isinstance(subclassify_phenotype, str):
            subclassify_phenotype = [subclassify_phenotype]
        if (showPhenotypeLabel):
            phenotype_label=phenotype+"_"+label
    elif phenotype is None:
         if isinstance(subclassify_phenotype, str) or (showPhenotypeLabel): 
            raise TypeError("You must pass a column name to the PHENOTYPE argument in order to use `subclassify_phenotype` or to set `showPhenotypeLabel = True`")
    
    
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
    if data.empty:
        raise TypeError("No cells were found to satisfy your `classify` criteria")
    else:
        # create new naming scheme for label and phenotype_label cols in classified
        classify_idx=data.index
        if showPhenotypeLabel is True:
            non_summary = pd.DataFrame({phenotype: adata.obs[phenotype]}) # gets the index and phenotype
            non_summary[phenotype] = non_summary[phenotype].astype(str)

            classified = pd.DataFrame(non_summary.loc[data.index]) #subsets phenotype rows to only classified cells
        
            classified[phenotype_label] = classified[phenotype]+"_"+classify_label # add phenotype_label col
            classified.drop([phenotype], axis='columns', inplace=True) # drop phenotype col, for merge        
        else:
            classified=pd.DataFrame(np.repeat(classify_label, len(classify_idx)),index= classify_idx, columns=[label]) # add label col



    if collapse_failed is True: 
        if showPhenotypeLabel is True:
            meta = non_summary # has index and phenotype col
        else:
            meta = pd.DataFrame(index= adata.obs.index)

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
            meta = meta[phenotype].astype("object")
            classified = pd.DataFrame(np.repeat(classify_label, len(classify_idx)), index = classify_idx, columns = [phenotype])
            classified = classified[phenotype].astype("object")
            meta.update(classified) # updates with label for only the classified cells
        
            
    # Add to Anndata 
    meta = meta.reindex(adata.obs.index)
    if showPhenotypeLabel is True:
        adata.obs[phenotype_label]=meta
    else:
        adata.obs[label]=meta 
    # return
    return adata
