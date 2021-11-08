#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Mon Oct 19 20:03:01 2020
# @author: Ajit Johnson Nirmal

"""
!!! abstract "Short Description"
    `sm.tl.spatial_interaction`: The function allows users to computes how likely celltypes are found next to each another
    compared to random background. 

## Function
"""

# Import library
import pandas as pd
from sklearn.neighbors import BallTree
import numpy as np
from joblib import Parallel, delayed
import scipy
from functools import reduce


# Function
def spatial_interaction (adata,x_coordinate='X_centroid',y_coordinate='Y_centroid',
                         phenotype='phenotype',
                         method='radius', radius=30, knn=10,
                         permutation=1000,
                         imageid='imageid',subset=None,
                         pval_method='zscore',
                         label='spatial_interaction'):
    """
Parameters:
    adata : AnnData object

    x_coordinate : float, required  
        Column name containing the x-coordinates values.

    y_coordinate : float, required  
        Column name containing the y-coordinates values.

    phenotype : string, required  
        Column name of the column containing the phenotype information. 
        It could also be any categorical assignment given to single cells.

    method : string, optional  
        Two options are available: a) 'radius', b) 'knn'.
        a) radius - Identifies the neighbours within a given radius for every cell.
        b) knn - Identifies the K nearest neigbours for every cell.

    radius : int, optional  
        The radius used to define a local neighbhourhood.

    knn : int, optional  
        Number of cells considered for defining the local neighbhourhood.

    permutation : int, optional  
        The number of permutations to be performed for calculating the P-Value.

    imageid : string, optional  
        Column name of the column containing the image id.

    subset : string, optional  
        imageid of a single image to be subsetted for analyis.

    pval_method : string, optional  
        Two options are available: a) 'histocat', b) 'zscore'.  
        a) P-values are calculated by subtracting the permuted mean from the observed mean
        divided by the number of permutations as described in the histoCAT manuscript (Denis et.al, Nature Methods 2017)  
        b) zscores are calculated from the mean and standard deviation and further p-values are
        derived by fitting the observed values to a normal distribution. The default is 'histocat'.

    label : string, optional  
        Key for the returned data, stored in `adata.obs`. The default is 'spatial_interaction'.

Returns:
    adata : AnnData object  
        Updated AnnData object with the results stored in `adata.obs['spatial_aggregate']`.
    
Example:
```python
    # Using the radius method to identify local neighbours and histocat to compute P-values
    adata = sm.tl.spatial_interaction(adata, method='radius', radius=30,pval_method='histocat',
                                      imageid='ImageId',x_coordinate='X_position',y_coordinate='Y_position')
    
    
    # Using the KNN method to identify local neighbours and zscore to compute P-values
    adata = sm.tl.spatial_interaction(adata, method='knn', radius=30,pval_method='zscore',
                                      imageid='ImageId',x_coordinate='X_position',y_coordinate='Y_position')
```
    """
    
    
    def spatial_interaction_internal (adata_subset,x_coordinate,y_coordinate,
                                      phenotype,
                                      method, radius, knn,
                                      permutation, 
                                      imageid,subset,
                                      pval_method):
        
        print("Processing Image: " + str(adata_subset.obs[imageid].unique()))
        
        # Create a dataFrame with the necessary inforamtion
        data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate], 'phenotype': adata_subset.obs[phenotype]})
        
        
        # Identify neighbourhoods based on the method used
        # a) KNN method
        if method == 'knn':
            print("Identifying the " + str(knn) + " nearest neighbours for every cell")
            tree = BallTree(data[['x','y']], leaf_size= 2)
            ind = tree.query(data[['x','y']], k=knn, return_distance= False)
            neighbours = pd.DataFrame(ind.tolist(), index = data.index) # neighbour DF
            neighbours.drop(0, axis=1, inplace=True) # Remove self neighbour
            
        # b) Local radius method
        if method == 'radius':
            print("Identifying neighbours within " + str(radius) + " pixels of every cell")
            kdt = BallTree(data[['x','y']], metric='euclidean') 
            ind = kdt.query_radius(data[['x','y']], r=radius, return_distance=False)
            for i in range(0, len(ind)): ind[i] = np.delete(ind[i], np.argwhere(ind[i] == i))#remove self
            neighbours = pd.DataFrame(ind.tolist(), index = data.index) # neighbour DF
            
        # Map Phenotypes to Neighbours
        # Loop through (all functionized methods were very slow)
        phenomap = dict(zip(list(range(len(ind))), data['phenotype'])) # Used for mapping
        print("Mapping phenotype to neighbors")
        for i in neighbours.columns:
            neighbours[i] = neighbours[i].dropna().map(phenomap, na_action='ignore')
            
        # Drop NA
        neighbours = neighbours.dropna(how='all')
        
        # Collapse all the neighbours into a single column
        n = pd.DataFrame(neighbours.stack(), columns = ["neighbour_phenotype"])
        n.index = n.index.get_level_values(0) # Drop the multi index
        
        # Merge with real phenotype
        n = n.merge(data['phenotype'], how='inner', left_index=True, right_index=True)
        
        # Permutation
        print('Performing '+ str(permutation) + ' permutations')
    
        def permutation_pval (data):
            data = data.assign(neighbour_phenotype=np.random.permutation(data['neighbour_phenotype']))
            #data['neighbour_phenotype'] = np.random.permutation(data['neighbour_phenotype'])
            data_freq = data.groupby(['phenotype','neighbour_phenotype']).size().unstack()
            data_freq = data_freq.fillna(0).stack().values 
            return data_freq
        
        # Apply function
        final_scores = Parallel(n_jobs=-1)(delayed(permutation_pval)(data=n) for i in range(permutation)) 
        perm = pd.DataFrame(final_scores).T
        
        # Consolidate the permutation results
        print('Consolidating the permutation results')
        # Calculate P value
        # real
        n_freq = n.groupby(['phenotype','neighbour_phenotype']).size().unstack().fillna(0).stack() 
        # permutation
        mean = perm.mean(axis=1)
        std = perm.std(axis=1)
        # P-value calculation
        if pval_method == 'histocat':
            # real value - prem value / no of perm 
            p_values = abs(n_freq.values - mean) / (permutation+1)
            p_values = p_values[~np.isnan(p_values)].values
        if pval_method == 'zscore':
            z_scores = (n_freq.values - mean) / std        
            z_scores[np.isnan(z_scores)] = 0
            p_values = scipy.stats.norm.sf(abs(z_scores))*2
            p_values = p_values[~np.isnan(p_values)]
            
        # Compute Direction of interaction (interaction or avoidance)
        direction = ((n_freq.values - mean) / abs(n_freq.values - mean)).fillna(1)

        # Normalize based on total cell count
        k = n.groupby(['phenotype','neighbour_phenotype']).size().unstack().fillna(0)
        # add neighbour phenotype that are not present to make k a square matrix
        columns_to_add = dict.fromkeys(np.setdiff1d(k.index,k.columns), 0)
        k = k.assign(**columns_to_add)

        total_cell_count = data['phenotype'].value_counts()
        total_cell_count = total_cell_count[k.columns].values # keep only cell types that are present in the column of k
        # total_cell_count = total_cell_count.reindex(k.columns).values # replaced by above
        k_max = k.div(total_cell_count, axis = 0)
        k_max = k_max.div(k_max.max(axis=1), axis=0).stack()
        
        # DataFrame with the neighbour frequency and P values
        count = (k_max.values * direction).values # adding directionallity to interaction
        neighbours = pd.DataFrame({'count': count,'p_val': p_values}, index = k_max.index)
        #neighbours.loc[neighbours[neighbours['p_val'] > p_val].index,'count'] = np.NaN
        #del neighbours['p_val']
        neighbours.columns = [adata_subset.obs[imageid].unique()[0], 'pvalue_' + str(adata_subset.obs[imageid].unique()[0])]
        neighbours = neighbours.reset_index()
        #neighbours = neighbours['count'].unstack()
        
        # return
        return neighbours
          
      
    # subset a particular subset of cells if the user wants else break the adata into list of anndata objects
    if subset is not None:
        adata_list = [adata[adata.obs[imageid] == subset]]
    else:
        adata_list = [adata[adata.obs[imageid] == i] for i in adata.obs[imageid].unique()]
    
    
    # Apply function to all images and create a master dataframe
    # Create lamda function 
    r_spatial_interaction_internal = lambda x: spatial_interaction_internal (adata_subset=x, x_coordinate=x_coordinate, y_coordinate=y_coordinate, phenotype=phenotype, method=method,  radius=radius, knn=knn, permutation=permutation, imageid=imageid,subset=subset,pval_method=pval_method) 
    all_data = list(map(r_spatial_interaction_internal, adata_list)) # Apply function 
    

    # Merge all the results into a single dataframe    
    df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['phenotype', 'neighbour_phenotype'], how='outer'), all_data)
    

    # Add to anndata
    adata.uns[label] = df_merged
    
    # return
    return adata


