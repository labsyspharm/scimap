#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 15:00:39 2020
@author: Ajit Johnson Nirmal
Function to compute the proportion of phenotypes or any single cell annotation 
within the local neighbourhood of each cell. 
"""

# Import library
import pandas as pd
import numpy as np
from sklearn.neighbors import BallTree


# Function
def spatial_count (adata,x_coordinate='X_centroid',y_coordinate='Y_centroid',
                   phenotype='phenotype',method='radius',radius=30,knn=10,
                   imageid='imageid',subset=None,label='spatial_count'):
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
    method : string, optional
        Two options are available: a) 'radius', b) 'knn'.
        a) radius - Identifies the neighbours within a given radius for every cell.
        b) knn - Identifies the K nearest neigbours for every cell.
        The default is 'radius'.
    radius : int, optional
        The radius used to define a local neighbhourhood. The default is 30.
    knn : int, optional
        Number of cells considered for defining the local neighbhourhood. The default is 10.
    imageid : string, optional
        Column name of the column containing the image id. The default is 'imageid'.
    subset : string, optional
        imageid of a single image to be subsetted for analyis. The default is None.
    label : string, optional
        Key for the returned data, stored in `adata.uns`. The default is 'spatial_count'.

    Returns
    -------
    adata : AnnData object
        Updated AnnData object with the results stored in `adata.uns['spatial_count']`.
    
    Example
    -------
    # Running the radius method
    adata = spatial_count (adata,x_coordinate='X_centroid',y_coordinate='Y_centroid',
                           phenotype='phenotype',method='radius',radius=30,
                           imageid='imageid',subset=None,label='spatial_count_radius')
    # Running the knn method
    adata = spatial_count (adata,x_coordinate='X_centroid',y_coordinate='Y_centroid',
                           phenotype='phenotype',method='knn',knn=10,
                           imageid='imageid',subset=None,label='spatial_count_knn')

    """

    def spatial_count_internal (adata_subset,x_coordinate,y_coordinate,phenotype,method,radius,knn,
                                imageid,subset,label):

        # Create a DataFrame with the necessary inforamtion
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
                                                   y_coordinate=y_coordinate,phenotype=phenotype,
                                                   method=method,radius=radius,knn=knn,
                                                   imageid=imageid,subset=subset,label=label) 
    all_data = list(map(r_spatial_count_internal, adata_list)) # Apply function 
    
    
    # Merge all the results into a single dataframe    
    result = []
    for i in range(len(all_data)):
        result.append(all_data[i])
    result = pd.concat(result, join='outer')  
    
    # Reindex the cells
    result = result.fillna(0)
    result = result.reindex(adata.obs.index)
    
    # Add to adata
    adata.uns[label] = result
    
    # Return        
    return adata