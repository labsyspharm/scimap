#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Wed Oct 14 09:23:07 2020
# @author: Ajit Johnson Nirmal
"""
!!! abstract "Short Description"
    `sm.tl.spatial_distance`: This function computes the average shortest distance 
    between specified phenotypes or clusters, supporting analysis of both 2D and 3D spatial data. 
    It facilitates the quantitative assessment of spatial relationships among cellular phenotypes 
    or clusters within tissue sections or 3D cultures.

## Function
"""

# Import library
import pandas as pd
from sklearn.neighbors import BallTree
from joblib import Parallel, delayed
import itertools

# Function
def spatial_distance (adata,
                      x_coordinate='X_centroid',
                      y_coordinate='Y_centroid',
                      z_coordinate=None,
                      phenotype='phenotype',
                      subset=None,
                      imageid='imageid',
                      verbose=True,
                      label='spatial_distance'):
    """
    
Parameters:
        adata (anndata.AnnData):  
            Annotated data matrix with spatial information.
        
        x_coordinate (str, required):  
            Column name containing x-coordinates.
        
        y_coordinate (str, required):  
            Column name containing y-coordinates.
        
        z_coordinate (str, optional):  
            Column name containing z-coordinates, for 3D spatial data analysis.
        
        phenotype (str, required):  
            Column name containing the phenotype information or any categorical cell classification.
        
        subset (str, optional):  
            Identifier for a subset of data to analyze, typically an image ID.
        
        imageid (str, optional):  
            Column name containing image identifiers, useful for analyzing distances within specific images.
        
        verbose (bool, optional):  
            If True, prints progress and informational messages during the calculation.
        
        label (str, optional):  
            Custom label for storing results in `adata.uns`.

Returns:
        adata (anndata.AnnData):  
            The input `adata` object, updated with the spatial distance results stored in `adata.uns[label]`.

Example:
    ```python
    
    # Calculate spatial distance in 2D
    adata = sm.tl.spatial_distance(adata, x_coordinate='X_centroid', y_coordinate='Y_centroid',
                             phenotype='cell_type', label='2D_distance')

    # Calculate spatial distance in 3D
    adata = sm.tl.spatial_distance(adata, x_coordinate='X_centroid', y_coordinate='Y_centroid',
                             z_coordinate='Z_centroid', phenotype='cell_type', label='3D_distance')

    # Calculate spatial distance for a specific image subset
    adata = sm.tl.spatial_distance(adata, x_coordinate='X_centroid', y_coordinate='Y_centroid',
                             phenotype='cell_type', imageid='image_id', subset='image01',
                             label='distance_image01')
    
    ```
    """
    
    
    def spatial_distance_internal (adata_subset,x_coordinate,y_coordinate,z_coordinate,
                                   phenotype,subset,imageid,label):
        
        if verbose:
            print("Processing Image: " + str(adata_subset.obs[imageid].unique()[0]))
        # Create a dataFrame with the necessary inforamtion
        if z_coordinate is not None:
            if verbose:
                print("Including Z -axis")
            data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate], 'z': adata_subset.obs[z_coordinate], 'phenotype': adata_subset.obs[phenotype]})
        else:
            data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate], 'phenotype': adata_subset.obs[phenotype]})

        # Function to identify shortest distance for each phenotype of interest
        def distance (pheno):
            pheno_interest = data[data['phenotype'] == pheno]
            # Build the ball-tree for search space
            tree = BallTree(pheno_interest[['x','y']], metric='euclidean') 
            # Calculate shortest distance (if statement to account for K)
            if len(pheno_interest) > 1:
                dist, ind = tree.query(data[['x','y']], k=2, return_distance= True)
                dist = pd.DataFrame(dist)
                dist.loc[dist[0] == 0, 0]  = dist[1]
                dist = dist[0].values
            else:
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
        adata_list = [adata[adata.obs[imageid].isin(subset)]]
    else:
        adata_list = [adata[adata.obs[imageid] == i] for i in adata.obs[imageid].unique()]
        
    # Apply function to all images and create a master dataframe
    # Create lamda function 
    r_spatial_distance_internal = lambda x: spatial_distance_internal (adata_subset=x,
                                                                       x_coordinate=x_coordinate,y_coordinate=y_coordinate, z_coordinate=z_coordinate,
                                                                       phenotype=phenotype,subset=subset,imageid=imageid,label=label) 
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
    
        

