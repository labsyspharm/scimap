#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 15:00:39 2020
@author: Ajit Johnson Nirmal
Function to compute a proximity based weighted expression scoring. The function generates
a neighbourhood for each cell and computes a expression score for all markers
for each cell based on its proximity to other cells within the neighbourhood.
"""

# Import library
import pandas as pd
import numpy as np
from sklearn.neighbors import BallTree
from scipy.sparse import csr_matrix

# Function
def spatial_expression (adata, x_coordinate='X_centroid',y_coordinate='Y_centroid',
                        method='radius', radius=30, knn=10, imageid='imageid', 
                        use_raw=True,subset=None,label='spatial_expression'):
    """
    

    Parameters
    ----------
    adata : AnnData object

    x_coordinate : float, required
        Column name containing the x-coordinates values. The default is 'X_centroid'.
    y_coordinate : float, required
        Column name containing the y-coordinates values. The default is 'Y_centroid'.
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
    use_raw : boolian, optional
        Argument to denote whether to use the raw data or scaled data after applying `sm.pp.rescale`.
        If `True`, the log of raw data is used. The default is True.
    label : string, optional
        Key for the returned data, stored in `adata.uns`. The default is 'spatial_count'.

    Returns
    -------
    adata : AnnData object
        Updated AnnData object with the results stored in `adata.uns['spatial_expression']`.
        
        
    Example
    -------
    # Running the radius method
    adata = sm.tl.spatial_expression (adata, x_coordinate='X_centroid',y_coordinate='Y_centroid',
                                method='radius', radius=30, imageid='imageid', 
                                use_raw=True,subset=None,label='spatial_expression_radius')
    # Running the knn method
    adata = sm.tl.spatial_expression (adata, x_coordinate='X_centroid',y_coordinate='Y_centroid',
                        method='knn', knn=10, imageid='imageid', 
                        use_raw=True,subset=None,label='spatial_expression_knn')

    """
    
    # Error statements
    if use_raw is False:
        if all(adata.X[0] < 1) is False:
            raise ValueError('Please run `sm.pp.rescale` first if you wish to use `use_raw = False`')
        
     
    def spatial_expression_internal (adata_subset, x_coordinate, y_coordinate,
                                     method, radius, knn, imageid, use_raw, subset,label):
         
        # Create a DataFrame with the necessary inforamtion
        data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate]})
        
        # Identify neighbourhoods based on the method used
        # a) KNN method
        if method == 'knn':
            print("Identifying the " + str(knn) + " nearest neighbours for every cell")
            tree = BallTree(data, leaf_size= 2)
            dist, ind = tree.query(data, k=knn, return_distance= True)
            # Remove self from ind
            neighbours = pd.DataFrame(ind.tolist(), index = data.index) # neighbour DF
            neighbours.drop(0, axis=1, inplace=True) # Remove self neighbour
            # Remove self from dist
            neighbours_dist = pd.DataFrame(dist.tolist(), index = data.index) # neighbour DF
            neighbours_dist.drop(0, axis=1, inplace=True) # Remove self neighbour
            # Merge the indeces and distance into a single dictionary
            ind_dict = lambda x,y: dict(zip(x, y)) # Function to create a dict of indeces and distances
            ind_dist = list(map(ind_dict, neighbours.values, neighbours_dist.values))
            
            
        # b) Local radius method
        if method == 'radius':
            print("Identifying neighbours within " + str(radius) + " pixels of every cell")
            kdt = BallTree(data, metric='euclidean')
            ind, dist = kdt.query_radius(data, r=radius, return_distance= True)
            # Remove self from ind
            for i in range(0, len(ind)): ind[i] = np.delete(ind[i], np.argwhere(ind[i] == i))#remove self
            for i in range(0, len(ind)): dist[i] =  dist[i][dist[i] != 0]#remove self distance
            # Merge the indeces and distance into a single dictionary
            ind_dict = lambda x,y: dict(zip(x, y)) # Function to create a dict of indeces and distances
            ind_dist = list(map(ind_dict, ind.tolist(),dist.tolist()))
        
        # Modify the dictionary to custom range for normalization 
        def range_normalize (data, upper_limit, lower_limit): # Function to normalize based on custom range
            
            if len(data) >= 2:
                # calculate data range
                x = (max (data.values()) - min (data.values())) / len (data)
                # calculate normalization factor range
                y = (upper_limit - lower_limit) / len (data)
                # Step-2: For each data point calculate the normalization factor
                xij = []
                for i in data:
                    xij = np.append(xij, ((max (data.values()) - data[i]) * y) / x)
                # Step-3: Compute the normalized values from the factors determined
                yij = []
                for j in xij:
                    yij = xij + lower_limit
                # modify the data object to reflect the new normlized values
                modified_data = dict(zip(data.keys(), yij))
            elif len(data) == 1:
                data[list(data.keys())[0]] = upper_limit
                modified_data = data
            elif data == {}:
                modified_data = {}
            
            return modified_data
        
        # Define range
        n = lambda x: range_normalize(x, 1, 0)
        # Apply function
        normalized_dist = list(map(n, ind_dist))
        
        # Normalize the spatial weights based on total number of neighbours per cell
        def norm (data):
            norm_data = {k: v/sum(data.values()) for k, v in data.items()}
            return norm_data
        
        # Run the function on all cell neighbourhoods  
        cell_number_normalized_dist = list(map(norm, normalized_dist))
        
        ## Spatial Lag calculation
        # Convert the dictionary to a martix
        idx = range(len(cell_number_normalized_dist))
        wn_matrix = pd.DataFrame(cell_number_normalized_dist, columns=idx, index=idx).fillna(0)
        wn_matrix.columns = data.index # add the cell name
        wn_matrix.index = data.index  # add the cell name
        
        # Convert the matrix to a sparce matrix
        wn_matrix_sparse = csr_matrix(wn_matrix)
        
        # Calculation of spatial lag
        if use_raw==True:
            spatial_lag = pd.DataFrame(wn_matrix_sparse * np.log1p(adata_subset.raw.X), columns = adata_subset.var.index, index=adata_subset.obs.index)
        else:
            spatial_lag = pd.DataFrame(wn_matrix_sparse * adata_subset.X, columns = adata_subset.var.index, index=adata_subset.obs.index)
        
        # return value
        return spatial_lag
    
    # Subset a particular image if needed
    if subset is not None:
        adata_list = [adata[adata.obs[imageid] == subset]]
    else:
        adata_list = [adata[adata.obs[imageid] == i] for i in adata.obs[imageid].unique()]
    
    # Apply function to all images and create a master dataframe
    # Create lamda function 
    r_spatial_expression_internal = lambda x: spatial_expression_internal(adata_subset=x, 
                                                                x_coordinate=x_coordinate, 
                                                                y_coordinate=y_coordinate, 
                                                                method=method, radius=radius, 
                                                                knn=knn, imageid=imageid, 
                                                                use_raw=use_raw, subset=subset,
                                                                label=label) 
    all_data = list(map(r_spatial_expression_internal, adata_list)) # Apply function 
    
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

