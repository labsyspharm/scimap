#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 09:23:07 2020
@author: Ajit Johnson Nirmal
This function will allow users to calculate the average shortest between every cell to phenotypes/clusters
of interest.  
"""

# Import library
import pandas as pd
from sklearn.neighbors import BallTree
from joblib import Parallel, delayed
import itertools

# Function

def spatial_distance (adata,x_coordinate='X_centroid',y_coordinate='Y_centroid',
                      phenotype='phenotype',subset=None,imageid='imageid',
                      label='spatial_distance'):
    """
    

    Parameters
    ----------
    adata : AnnData object

    x_coordinate : float, required
        Column name containing the x-coordinates values. The default is 'X_centroid'.
    y_coordinate : float, required
        Column name containing the y-coordinates values. The default is 'Y_centroid'.
    phenotype : string, required
        Column name of the column containing the phenotype information. 
        It could also be any categorical assignment given to single cells. The default is 'phenotype'.
    subset : string, optional
        imageid of a single image to be subsetted for analyis. The default is None.
    imageid : string, optional
        Column name of the column containing the image id. The default is 'imageid'.
    label : string, optional
        Key for the returned data, stored in `adata.obs`. The default is 'spatial_distance'.

    Returns
    -------
    adata : AnnData object
        Updated AnnData object with the results stored in `adata.uns['spatial_distance']`.
        
    Example
    -------
    adata = sm.tl.spatial_distance (adata,x_coordinate='X_position',y_coordinate='Y_position',imageid='ImageId')
        

    """
    
    
    def spatial_distance_internal (adata_subset,x_coordinate,y_coordinate,
                                   phenotype,subset,imageid,label):
        
        print("Processing Image: " + str(adata_subset.obs[imageid].unique()[0]))
        # Create a dataFrame with the necessary inforamtion
        data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate], 'phenotype': adata_subset.obs[phenotype]})

        # Function to identify shortest distance for each phenotype of interest
        def distance (pheno):
            pheno_interest = data[data['phenotype'] == pheno]
            # Build the ball-tree for search space
            tree = BallTree(pheno_interest[['x','y']], metric='euclidean') 
            # Calculate shortest distance
            dist, ind = tree.query(data[['x','y']], k=1, return_distance= True)
            dist = list(itertools.chain.from_iterable(dist))
            return dist
        
        # Run in parallel for all phenotypes
        phenotype_list = list(data['phenotype'].unique())
        # Apply function
        final_dist = Parallel(n_jobs=-1)(delayed(distance)(pheno=i) for i in phenotype_list)     
        final_dist = pd.DataFrame(final_dist, index = phenotype_list, columns = data.index).T
        
        return final_dist
    
    # subset a particular subset of cells if the user wants else break the adata into list of anndata objects
    if subset is not None:
        adata_list = [adata[adata.obs[imageid] == subset]]
    else:
        adata_list = [adata[adata.obs[imageid] == i] for i in adata.obs[imageid].unique()]
        
    # Apply function to all images and create a master dataframe
    # Create lamda function 
    r_spatial_distance_internal = lambda x: spatial_distance_internal (adata_subset=x,x_coordinate=x_coordinate,y_coordinate=y_coordinate, phenotype=phenotype,subset=subset,imageid=imageid,label=label) 
    all_data = list(map(r_spatial_distance_internal, adata_list)) # Apply function 
    
    # Merge all the results into a single dataframe    
    result = []
    for i in range(len(all_data)):
        result.append(all_data[i])
    result = pd.concat(result, join='outer')  
    
    
    # Add to anndata
    adata.uns[label] = result
    
    # return
    return adata
    
        

