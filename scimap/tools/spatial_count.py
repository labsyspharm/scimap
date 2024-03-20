#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Wed Aug 19 15:00:39 2020
# @author: Ajit Johnson Nirmal
"""
!!! abstract "Short Description"
    `sm.tl.spatial_count` computes a neighborhood matrix from spatial data using categorical variables, 
    such as cell types, to identify local cell clusters. It offers two neighborhood definition methods:

    - **Radius Method**: Identifies neighbors within a specified radius for each cell, allowing for 
    the exploration of spatial relationships based on physical proximity.
    - **KNN Method**: Determines neighbors based on the K nearest neighbors, focusing on the closest 
    spatial associations irrespective of physical distance.
    
    The generated neighborhood matrix is stored in `adata.uns`, providing a basis for further analysis. 
    To uncover Recurrent Cellular Neighborhoods (RCNs) that share similar spatial patterns, users can 
    cluster the neighborhood matrix using the `spatial_cluster` function. This approach enables the 
    identification of spatially coherent cell groups, facilitating insights into the cellular 
    architecture of tissues.

## Function
"""

# Import library
import pandas as pd
import numpy as np
from sklearn.neighbors import BallTree

# Function
def spatial_count (adata,
                   x_coordinate='X_centroid',
                   y_coordinate='Y_centroid',
                   z_coordinate=None,
                   phenotype='phenotype',
                   method='radius',
                   radius=30,knn=10,
                   imageid='imageid',
                   subset=None,
                   verbose=True,
                   label='spatial_count'):
    """
Parameters:
        adata (anndata.AnnData):  
            Annotated data matrix with spatial information.
        
        x_coordinate (str, required):  
            Column name containing x-coordinates.
        
        y_coordinate (str, required):  
            Column name containing y-coordinates.
        
        z_coordinate (str, optional):  
            Column name containing z-coordinates, for 3D spatial data.
        
        phenotype (str, required):  
            Column name containing phenotype or any categorical cell classification.
        
        method (str, optional):  
            Neighborhood definition method: 'radius' for fixed distance, 'knn' for K nearest neighbors.
        
        radius (int, optional):  
            Radius used to define neighborhoods (applicable when method='radius').
        
        knn (int, optional):  
            Number of nearest neighbors to consider (applicable when method='knn').
        
        imageid (str, optional):  
            Column name containing image identifiers, for analyses limited to specific images.
        
        subset (str, optional):  
            Specific image identifier for subsetting data before analysis.
        
        verbose (bool, optional):  
            If True, prints progress and informational messages.
        
        label (str, optional):  
            Key for storing results in `adata.uns`.

Returns:
        adata (anndata.AnnData):  
            Updated AnnData object with the neighborhood matrix stored in `adata.uns[label]`.

Example:
    ```python
    
    # Analyze spatial relationships using the radius method
    adata = sm.tl.spatial_count(adata, x_coordinate='X_centroid', y_coordinate='Y_centroid',
                          phenotype='phenotype', method='radius', radius=50,
                          label='neighborhood_radius50')

    # Explore spatial neighborhoods with KNN
    adata = sm.tl.spatial_count(adata, x_coordinate='X_centroid', y_coordinate='Y_centroid',
                          phenotype='phenotype', method='knn', knn=15,
                          label='neighborhood_knn15')

    # 3D spatial analysis using a radius method
    adata = sm.tl.spatial_count(adata, x_coordinate='X_centroid', y_coordinate='Y_centroid',
                          z_coordinate='Z_centroid', phenotype='phenotype', method='radius', radius=30,
                          label='neighborhood_3D_radius30')
    
    ```
    """

    def spatial_count_internal (adata_subset,x_coordinate,y_coordinate,z_coordinate,phenotype,method,radius,knn,
                                imageid,subset,label):
        
        # Create a dataFrame with the necessary inforamtion
        if z_coordinate is not None:
            if verbose:
                print("Including Z -axis")
            data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate], 'z': adata_subset.obs[z_coordinate], 'phenotype': adata_subset.obs[phenotype]})
        else:
            data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate], 'phenotype': adata_subset.obs[phenotype]})


        # Create a DataFrame with the necessary inforamtion
        #data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate], 'phenotype': adata_subset.obs[phenotype]})
        
        # Identify neighbourhoods based on the method used
        # a) KNN method
        if method == 'knn':
            if verbose:
                print("Identifying the " + str(knn) + " nearest neighbours for every cell")
            if z_coordinate is not None:
                tree = BallTree(data[['x','y','z']], leaf_size= 2)
                ind = tree.query(data[['x','y','z']], k=knn, return_distance= False)
            else:
                tree = BallTree(data[['x','y']], leaf_size= 2)
                ind = tree.query(data[['x','y']], k=knn, return_distance= False)
            neighbours = pd.DataFrame(ind.tolist(), index = data.index) # neighbour DF
            neighbours.drop(0, axis=1, inplace=True) # Remove self neighbour
        
        # b) Local radius method
        if method == 'radius':
            if verbose:
                print("Identifying neighbours within " + str(radius) + " pixels of every cell")
            if z_coordinate is not None:
                kdt = BallTree(data[['x','y','z']], metric='euclidean') 
                ind = kdt.query_radius(data[['x','y','z']], r=radius, return_distance=False)
            else:
                kdt = BallTree(data[['x','y']], metric='euclidean') 
                ind = kdt.query_radius(data[['x','y']], r=radius, return_distance=False)
                
            for i in range(0, len(ind)): ind[i] = np.delete(ind[i], np.argwhere(ind[i] == i))#remove self
            neighbours = pd.DataFrame(ind.tolist(), index = data.index) # neighbour DF
            
        # Map phenotype
        phenomap = dict(zip(list(range(len(ind))), data['phenotype'])) # Used for mapping
        
        # Loop through (all functionized methods were very slow)
        for i in neighbours.columns:
            neighbours[i] = neighbours[i].dropna().map(phenomap, na_action='ignore')
        
        # Drop NA
        #n_dropped = neighbours.dropna(how='all')
           
        # Collapse all the neighbours into a single column
        n = pd.DataFrame(neighbours.stack(), columns = ["neighbour_phenotype"])
        n.index = n.index.get_level_values(0) # Drop the multi index
        n = pd.DataFrame(n)
        n['order'] = list(range(len(n)))
        
        # Merge with real phenotype
        n_m = n.merge(data['phenotype'], how='inner', left_index=True, right_index=True)
        n_m['neighbourhood'] = n_m.index
        n = n_m.sort_values(by=['order'])
        
        # Normalize based on total cell count
        k = n.groupby(['neighbourhood','neighbour_phenotype']).size().unstack().fillna(0)
        k = k.div(k.sum(axis=1), axis=0)
        
        # return the normalized neighbour occurance count
        return k
    
    # Subset a particular image if needed
    if subset is not None:
        adata_list = [adata[adata.obs[imageid] == subset]]
    else:
        adata_list = [adata[adata.obs[imageid] == i] for i in adata.obs[imageid].unique()]
    
    # Apply function to all images and create a master dataframe
    # Create lamda function 
    r_spatial_count_internal = lambda x: spatial_count_internal(adata_subset=x,x_coordinate=x_coordinate,
                                                   y_coordinate=y_coordinate,
                                                   z_coordinate=z_coordinate,
                                                   phenotype=phenotype,
                                                   method=method,radius=radius,knn=knn,
                                                   imageid=imageid,subset=subset,label=label) 
    all_data = list(map(r_spatial_count_internal, adata_list)) # Apply function 
    
    
    # Merge all the results into a single dataframe    
    result = []
    for i in range(len(all_data)):
        result.append(all_data[i])
    result = pd.concat(result, join='outer')  
    
    # Reindex the cells
    result = result.reindex(adata.obs.index)
    result = result.fillna(0)
    
    # Add to adata
    adata.uns[label] = result
    
    # Return        
    return adata