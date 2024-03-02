#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Mon Mar  2 19:56:08 2020
# @author: Ajit Johnson Nirmal
"""
!!! abstract "Short Description"
    `sm.tl.phenotype_cells`: This function annotates each cell in the dataset with a phenotype based on `scaled data` and a predefined `phenotype workflow`. Before using this function, ensure the data is scaled with the `sm.tl.rescale` function.
    
    *Description of the Phenotype Workflow File:*  
    Find an example `phenotype_workflow.csv` [here](https://github.com/ajitjohnson/scimap/blob/master/scimap/tests/_data/phenotype_workflow.csv).
    
    The `phenotype_workflow` file outlines six gating strategies for cell phenotyping:
    
    - **allpos**: Requires all specified markers to be positive for a cell to be assigned the phenotype.
    - **allneg**: Requires all specified markers to be negative for a cell to be assigned the phenotype.
    - **anypos**: Requires at least one of the specified markers to be positive for a cell to be assigned the phenotype. For example, a macrophage could be identified if it is positive for any of the markers `CD68`, `CD163`, or `CD206`.
    - **anyneg**: Requires at least one of the specified markers to be negative for a cell to be assigned the phenotype.
    - **pos**: Specifies that a cell must be positive for the given marker(s) to be assigned the phenotype. If used for multiple markers, cells not meeting all criteria may still be classified as a potential phenotype, allowing for later refinement by the user. For instance, regulatory T cells could be defined as `CD4+` and `FOXP3+`; cells not fully meeting these criteria might be labeled as likely-regulatory T cells for further evaluation.
    - **neg**: Specifies that a cell must be negative for the given marker(s) to be assigned the phenotype.
    
    *Recommendation*: Prioritize using positive markers to define phenotypes whenever possible.

## Function
"""

# Library
import numpy as np
import pandas as pd
import argparse


def phenotype_cells (adata, 
                     phenotype, 
                     gate = 0.5, 
                     label="phenotype", 
                     imageid='imageid',
                     pheno_threshold_percent=None, 
                     pheno_threshold_abs=None,
                     verbose=True
                     ):
    """
    
Parameters:
    adata (anndata.AnnData):  
        The input AnnData object containing single-cell data for phenotyping.

    phenotype (pd.DataFrame):  
        A DataFrame specifying the gating strategy for cell phenotyping. It should outline the workflow for phenotype classification based on marker expression levels. An example workflow is available at [this GitHub link](https://github.com/ajitjohnson/scimap/blob/master/scimap/tests/_data/phenotype_workflow.csv).
        
    gate (float, optional):  
        The threshold value for determining positive cell classification based on scaled data. By convention, values above this threshold are considered to indicate positive cells. 
        
    label (str):  
        The name of the column in `adata.obs` where the final phenotype classifications will be stored. This label will be used to access the phenotyping results within the `AnnData` object.
        
    imageid (str, optional):  
        The name of the column in `adata.obs` that contains unique image identifiers. This is crucial for analyses that require differentiation of data based on the source image, especially when using phenotype threshold parameters (`pheno_threshold_percent` or `pheno_threshold_abs`).
        
    pheno_threshold_percent (float, optional):  
        A threshold value (between 0 and 100) specifying the minimum percentage of cells that must exhibit a particular phenotype for it to be considered valid. Phenotypes not meeting this threshold are reclassified as 'unknown'. This parameter is useful for minimizing the impact of low-frequency false positives. 
        
    pheno_threshold_abs (int, optional):  
        Similar to `pheno_threshold_percent`, but uses an absolute cell count instead of a percentage. Phenotypes with cell counts below this threshold are reclassified as 'unknown'. This can help in addressing rare phenotype classifications that may not be meaningful. 
    
    verbose (bool):  
        If set to `True`, the function will print detailed messages about its progress and the steps being executed.

Returns:
    adata (anndata.AnnData):  
        The input AnnData object, updated to include the phenotype classifications for each cell. The phenotyping results can be found in `adata.obs[label]`, where `label` is the name specified by the user for the phenotype column.

Example:    
    ```python
    
    # Load the phenotype workflow CSV file
    phenotype = pd.read_csv('path/to/csv/file/')  
    
    # Apply phenotyping to cells based on the specified workflow
    adata = sm.tl.phenotype_cells(adata, phenotype=phenotype, gate=0.5, label="phenotype")
    
    ```

    """

    # Create a dataframe from the adata object
    data = pd.DataFrame(adata.X, columns = adata.var.index, index= adata.obs.index)

    # Function to calculate the phenotype scores
    def phenotype_cells (data,phenotype,gate,group):

        # Subset the phenotype based on the group
        phenotype = phenotype[phenotype.iloc[:,0] == group]

        # Parser to parse the CSV file into four categories
        def phenotype_parser (p, cell):
            # Get the index and subset the phenotype row being passed in
            location = p.iloc[:,1] == cell
            idx = [i for i, x in enumerate(location) if x][0]
            phenotype = p.iloc[idx,:]
            # Calculate
            pos = phenotype[phenotype == 'pos'].index.tolist()
            neg = phenotype[phenotype == 'neg'].index.tolist()
            anypos = phenotype[phenotype == 'anypos'].index.tolist()
            anyneg = phenotype[phenotype == 'anyneg'].index.tolist()
            allpos = phenotype[phenotype == 'allpos'].index.tolist()
            allneg = phenotype[phenotype == 'allneg'].index.tolist()
            return {'pos': pos, 'neg': neg ,'anypos': anypos, 'anyneg': anyneg, 'allpos': allpos, 'allneg': allneg}
            #return pos, neg, anypos, anyneg

        # Run the phenotype_parser function on all rows
        p_list = phenotype.iloc[:,1].tolist()
        r_phenotype = lambda x: phenotype_parser(cell=x, p=phenotype) # Create lamda function
        all_phenotype = list(map(r_phenotype, p_list)) # Apply function
        all_phenotype = dict(zip(p_list, all_phenotype)) # Name the lists

        # Define function to check if there is any marker that does not satisfy the gate
        def gate_satisfation_lessthan (marker, data, gate):
            fail = np.where(data[marker] < gate, 1, 0) # 1 is fail
            return fail
        # Corresponding lamda function
        r_gate_satisfation_lessthan = lambda x: gate_satisfation_lessthan(marker=x, data=data, gate=gate)

        # Define function to check if there is any marker that does not satisfy the gate
        def gate_satisfation_morethan (marker, data, gate):
            fail = np.where(data[marker] > gate, 1, 0)
            return fail
        # Corresponding lamda function
        r_gate_satisfation_morethan = lambda x: gate_satisfation_morethan(marker=x, data=data, gate=gate)

        def prob_mapper (data, all_phenotype, cell, gate):
            
            if verbose:
                print("Phenotyping " + str(cell))

            # Get the appropriate dict from all_phenotype
            p = all_phenotype[cell]

            # Identiy the marker used in each category
            pos = p.get('pos')
            neg = p.get('neg')
            anypos = p.get('anypos')
            anyneg = p.get('anyneg')
            allpos = p.get('allpos')
            allneg = p.get('allneg')

            # Perform computation for each group independently
            # Positive marker score
            if len(pos) != 0:
                pos_score = data[pos].mean(axis=1).values
                pos_fail = list(map(r_gate_satisfation_lessthan, pos)) if len(pos) > 1 else []
                pos_fail = np.amax(pos_fail, axis=0) if len(pos) > 1 else []
            else:
                pos_score = np.repeat(0, len(data))
                pos_fail = []

            # Negative marker score
            if len(neg) != 0:
                neg_score = (1-data[neg]).mean(axis=1).values
                neg_fail = list(map(r_gate_satisfation_morethan, neg)) if len(neg) > 1 else []
                neg_fail = np.amax(neg_fail, axis=0) if len(neg) > 1 else []
            else:
                neg_score = np.repeat(0, len(data))
                neg_fail = []

            # Any positive score
            anypos_score = np.repeat(0, len(data)) if len(anypos) == 0 else data[anypos].max(axis=1).values

            # Any negative score
            anyneg_score = np.repeat(0, len(data)) if len(anyneg) == 0 else (1-data[anyneg]).max(axis=1).values

            # All positive score
            if len(allpos) != 0:
                allpos_score = data[allpos]
                allpos_score['score'] = allpos_score.max(axis=1)
                allpos_score.loc[(allpos_score < gate).any(axis = 1), 'score'] = 0
                allpos_score = allpos_score['score'].values + 0.01 # A small value is added to give an edge over the matching positive cell
            else:
                allpos_score = np.repeat(0, len(data))


            # All negative score
            if len(allneg) != 0:
                allneg_score = 1- data[allneg]
                allneg_score['score'] = allneg_score.max(axis=1)
                allneg_score.loc[(allneg_score < gate).any(axis = 1), 'score'] = 0
                allneg_score = allneg_score['score'].values + 0.01
            else:
                allneg_score = np.repeat(0, len(data))


            # Total score calculation
            # Account for differences in the number of categories used for calculation of the final score
            number_of_non_empty_features = np.sum([len(pos) != 0,
                                                len(neg) != 0,
                                                len(anypos) != 0,
                                                len(anyneg) != 0,
                                                len(allpos) != 0,
                                                len(allneg) != 0])

            total_score = (pos_score + neg_score + anypos_score + anyneg_score + allpos_score + allneg_score) / number_of_non_empty_features

            return {cell: total_score, 'pos_fail': pos_fail ,'neg_fail': neg_fail}
            #return total_score, pos_fail, neg_fail


        # Apply the fuction to get the total score for all cell types
        r_prob_mapper = lambda x: prob_mapper (data=data, all_phenotype=all_phenotype, cell=x, gate=gate) # Create lamda function
        final_scores = list(map(r_prob_mapper, [*all_phenotype])) # Apply function
        final_scores = dict(zip([*all_phenotype], final_scores)) # Name the lists

        # Combine the final score to annotate the cells with a label
        final_score_df = pd.DataFrame()
        for i in [*final_scores]:
            df = pd.DataFrame(final_scores[i][i])
            final_score_df= pd.concat([final_score_df, df], axis=1)
        # Name the columns
        final_score_df.columns = [*final_scores]
        final_score_df.index = data.index
        # Add a column called unknown if all markers have a value less than the gate (0.5)
        unknown = group + str('-rest')
        final_score_df[unknown] = (final_score_df < gate).all(axis=1).astype(int)

        # Name each cell
        labels = final_score_df.idxmax(axis=1)

        # Group all failed instances (i.e. when multiple markers were given
        # any one of the marker fell into neg or pos zones of the gate)
        pos_fail_all = pd.DataFrame()
        for i in [*final_scores]:
            df = pd.DataFrame(final_scores[i]['pos_fail'])
            df.columns = [i] if len(df) != 0 else []
            pos_fail_all= pd.concat([pos_fail_all, df], axis=1)
        pos_fail_all.index = data.index if len(pos_fail_all) != 0 else []
        # Same for Neg
        neg_fail_all = pd.DataFrame()
        for i in [*final_scores]:
            df = pd.DataFrame(final_scores[i]['neg_fail'])
            df.columns = [i] if len(df) != 0 else []
            neg_fail_all= pd.concat([neg_fail_all, df], axis=1)
        neg_fail_all.index = data.index if len(neg_fail_all) != 0 else []


        # Modify the labels with the failed annotations
        if len(pos_fail_all) != 0:
            for i in pos_fail_all.columns:
                labels[(labels == i) & (pos_fail_all[i] == 1)] = 'likely-' + i
        # Do the same for negative
        if len(neg_fail_all) != 0:
            for i in neg_fail_all.columns:
                labels[(labels == i) & (neg_fail_all[i] == 1)] = 'likely-' + i

        # Retun the labels
        return labels

    # Create an empty dataframe to hold the labeles from each group
    phenotype_labels = pd.DataFrame()

    # Loop through the groups to apply the phenotype_cells function
    for i in phenotype.iloc[:,0].unique():

        if phenotype_labels.empty:
            phenotype_labels = pd.DataFrame(phenotype_cells(data = data, group = i, phenotype=phenotype, gate=gate))
            phenotype_labels.columns = [i]

        else:
            # Find the column with the cell-type of interest
            column_of_interest = [] # Empty list to hold the column name
            try:
                column_of_interest = phenotype_labels.columns[phenotype_labels.eq(i).any()]
            except:
                pass
            # If the cell-type of interest was not found just add NA
            if len(column_of_interest) == 0:
                phenotype_labels[i] = np.nan
            else:
                #cells_of_interest = phenotype_labels[phenotype_labels[column_of_interest] == i].index
                cells_of_interest = phenotype_labels[phenotype_labels[column_of_interest].eq(i).any(axis=1)].index
                d = data.loc[cells_of_interest]
                if verbose:
                    print("-- Subsetting " + str(i))
                phenotype_l = pd.DataFrame(phenotype_cells(data = d, group = i, phenotype=phenotype, gate=gate), columns = [i])
                phenotype_labels = phenotype_labels.merge(phenotype_l, how='outer', left_index=True, right_index=True)

    # Rearrange the rows back to original
    phenotype_labels = phenotype_labels.reindex(data.index)
    phenotype_labels = phenotype_labels.replace('-rest', np.nan, regex=True)

    if verbose:
        print("Consolidating the phenotypes across all groups")
    phenotype_labels_Consolidated = phenotype_labels.fillna(method='ffill', axis = 1)
    phenotype_labels[label] = phenotype_labels_Consolidated.iloc[:,-1].values

    # replace nan to 'other cells'
    phenotype_labels[label] = phenotype_labels[label].fillna('Unknown')

    # Apply the phenotype threshold if given
    if pheno_threshold_percent or pheno_threshold_abs is not None:
        p = pd.DataFrame(phenotype_labels[label])
        q = pd.DataFrame(adata.obs[imageid])
        p = q.merge(p, how='outer', left_index=True, right_index=True)

        # Function to remove phenotypes that are less than the given threshold
        def remove_phenotype(p, ID, pheno_threshold_percent, pheno_threshold_abs):
            d = p[p[imageid] == ID]
            x = pd.DataFrame(d.groupby([label]).size())
            x.columns = ['val']
            # FInd the phenotypes that are less than the given threshold
            if pheno_threshold_percent is not None:
                fail = list(x.loc[x['val'] < x['val'].sum() * pheno_threshold_percent/100].index)
            if pheno_threshold_abs is not None:
                fail = list(x.loc[x['val'] < pheno_threshold_abs].index)
            d[label] = d[label].replace(dict(zip(fail, np.repeat('Unknown',len(fail)))))
            # Return
            return d

        # Apply function to all images
        r_remove_phenotype = lambda x: remove_phenotype (p=p, ID=x,
                                                         pheno_threshold_percent=pheno_threshold_percent,
                                                         pheno_threshold_abs=pheno_threshold_abs) # Create lamda function
        final_phrnotypes= list(map(r_remove_phenotype, list(p[imageid].unique()))) # Apply function

        final_phrnotypes = pd.concat(final_phrnotypes, join='outer')
        phenotype_labels = final_phrnotypes.reindex(adata.obs.index)


    # Return to adata
    adata.obs[label] = phenotype_labels[label]

    #for i in phenotype_labels.columns:
    #    adata.obs[i] = phenotype_labels[i]

    return adata

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='The phenotyping function takes in the `scaled data` and a prior knowledge based `phenotype workflow` file to assign phenotype annotation to each cell in the dataset')
    parser.add_argument('--adata', type=str, help='anndata object')
    parser.add_argument('--phenotype', type=str, default='phenotype', help='A gating strategy for phenotyping the cells. An example `workflow` provided [here](https://github.com/ajitjohnson/scimap/blob/master/scimap/tests/_data/phenotype_workflow.csv).')
    parser.add_argument('--gate', type=int, default=0.5, help='By default rescale function, scales the data such that values above 0.5 are considered positive cells.')
    parser.add_argument('--label', type=str, default='phenotype', help='Name the column underwhich the final phenotype calling will be saved.')
    parser.add_argument('--imageid', type=str, default='imageid', help='Name of the column that contains the unique imageid')
    parser.add_argument('--pheno_threshold_percent', type=float, default=True, help='Accepts values between (0-100). If any particular phenotype is below the user defined threshold, it is recategorised as unknown')
    parser.add_argument('--pheno_threshold_abs', type=int, default=None, help='Serves the same purpose as that of pheno_threshold_percent. However, an absolute number can be passed')
    parser.add_argument('--verbose', required=False, default=True, help='The function will print detailed messages about its progress.')
    args = parser.parse_args()
    
    phenotype_cells(adata=args.adata,
                   phenotype=args.phenotype, 
                   gate=args.gate,
                   label=args.label, 
                   imageid=args.imageid, 
                   phenotype_threshhold_percent=args.pheno_threshold_percent, 
                   pheno_threshold_abs=args.pheno_threshold_abs)
    