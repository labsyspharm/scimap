#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 19:56:08 2020
@author: Ajit Johnson Nirmal
Cell Phenotyping
"""

# Library
import numpy as np
import pandas as pd


def phenotype_cells (adata, 
                     phenotype, 
                     gate = 0.5, 
                     label="phenotype", 
                     imageid='imageid',
                     pheno_threshold_percent=None, 
                     pheno_threshold_abs=None):
    """
    
    !!! example "Function Call"
        `sm.tl.phenotype_cells` (
          adata, phenotype,  
          gate = 0.5,  
          label="phenotype",  
          imageid='imageid',  
          pheno_threshold_percent=None,  
          pheno_threshold_abs=None)
    
    
    Short description
    ----------
    The phenotyping function takes in the `scaled data` and a prior knowledge based `phenotype workflow` 
    file to assign phenotype annotation to each cell in the dataset. Use the `sm.tl.rescale` function to
    rescale the data first. 
    
    *Phenotype workflow file description:*  
    An example of the `phenotype_workflow.csv` can be found [here](https://github.com/ajitjohnson/scimap/blob/master/scimap/tests/_data/phenotype_workflow.csv).  
    
    The `phenotype_workflow` accepts six categories of gating strategy for performing phenotyping.
    
    - allpos
    - allneg
    - anypos
    - anyneg
    - pos
    - neg
    
    `allpos`- All of the defined markers should be positive.  
    `allneg`- All of the defined markers should be negative.  
    `anypos`- Any one of the defined marker is sufficient to be positive. (e.g) For defining macrophages, one could use a strategy in which a cell is defined as a macrophage if any of `CD68, CD163 or CD206` is positive.  
    `anyneg`- Any of the defined marker is negative.  
    `pos`- A given marker is positive. If this argument is passed to multiple markers. (e.g) If regulatory T cell is defined as `CD4+`, `FOXP3+` by passing `pos` to each the markers and the algorithm finds that for a few cells one of the two is not, the algorithm will assign the cell as likely-regulatory T cell and will allow the user to make the decision later.  
    `neg`- A given marker is negative.  
    *It is always advised to use positive markers over negative markers*  

    
    Parameters
    ----------
    `adata` : anndata object
    
    `phenotype` : dataframe, required
        A gating strategy for phenotyping the cells. An example `workflow` provided [here](https://github.com/ajitjohnson/scimap/blob/master/scimap/tests/_data/phenotype_workflow.csv).
        
    `gate` : int, optional
        By default rescale function, scales the data such that values above 0.5 are considered positive cells.
        The default is `0.5`.
        
    `label` : string, optional
        Name the column underwhich the final phenotype calling will be saved. The default is `phenotype`.
        
    `imageid` : string, optional
        Name of the column that contains the unique imageid. This is only utilized
        when the user uses `pheno_threshold_percent` or `pheno_threshold_abs` parameters.
        The default is `imageid`.
        
    `pheno_threshold_percent` : float, optional
        Accepts values between (0-100). If any particular phenotype is below the user defined threshold,
        it is recategorised as 'unknown. Generally used to deal with low background false positives.
        The default is `None`.
        
    `pheno_threshold_abs` : int, optional
        Serves the same purpose as that of pheno_threshold_percent. However, an absolute
        number can be passed. For example, if user passes in 10- any phenotype that contains
        less than 10 cells will be recategorized as unknown. The default is `None`.

    Returns
    -------
    `adata`
        Updated AnnData object with the phenotype calls for each cell. Check `adata.obs['phenotype']` for results.

    Example
    -------
    
    ```
    # load the workflow csv file
    phenotype = pd.read_csv('path/to/csv/file/')  
    # phenotype the cells based on the workflow provided
    adata = sm.tl.phenotype_cells (adata, phenotype=phenotype, gate = 0.5, label="phenotype")
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
                print("-- Subsetting " + str(i))
                phenotype_l = pd.DataFrame(phenotype_cells(data = d, group = i, phenotype=phenotype, gate=gate), columns = [i])
                phenotype_labels = phenotype_labels.merge(phenotype_l, how='outer', left_index=True, right_index=True)

    # Rearrange the rows back to original
    phenotype_labels = phenotype_labels.reindex(data.index)
    phenotype_labels = phenotype_labels.replace('-rest', np.nan, regex=True)

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
