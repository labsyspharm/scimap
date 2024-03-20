#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Wed Aug 19 15:00:39 2020
# @author: Ajit Johnson Nirmal
"""
!!! abstract "Short Description"
    `sm.tl.spatial_aggregate`: This function identifies spatial clusters of phenotypically similar cells 
    within specified regions. By adjusting the `purity` parameter, users can specify the minimum 
    percentage of similarity required among cells within a defined `radius` or nearest neighbors, 
    enabling precise delineation of cellular aggregates.

## Function
"""

# Import library
import pandas as pd
import numpy as np
from sklearn.neighbors import BallTree

# Function
def spatial_aggregate (adata, 
                       x_coordinate='X_centroid',
                       y_coordinate='Y_centroid',
                       z_coordinate= None,
                       purity = 60, 
                       phenotype='phenotype', 
                       method='radius', 
                       radius=30, 
                       knn=10, 
                       imageid='imageid',
                       subset=None,
                       verbose=True,
                       label='spatial_aggregate'):
    """
    
Parameters:
        adata (anndata.AnnData):  
            The annotated data matrix of shape (n_obs, n_vars), where rows correspond to cells and columns to genes, used for spatial analysis.

        x_coordinate (str, required):  
            The column name in `adata` containing the x-coordinates of cells.

        y_coordinate (str, required):  
            The column name in `adata` containing the y-coordinates of cells.
        
        z_coordinate (str, required):  
            The column name in `adata` containing the z-coordinates of cells.

        purity (int, optional):  
            The minimum percentage (between 1 and 100) of cells with a similar phenotype required in a neighborhood for it to be considered a cluster. 

        phenotype (str, required):  
            The column name in `adata` representing cell phenotype information or any other categorical classification of cells.

        method (str, optional):  
            The neighborhood definition method: 'radius' for radial distance-based neighborhoods or 'knn' for k-nearest neighbors-based neighborhoods. 

        radius (int, optional):  
            The radius used to define a neighborhood around each cell, applicable when `method='radius'`. Measured in the same units as x and y coordinates.

        knn (int, optional):  
            The number of nearest neighbors to consider for defining a neighborhood around each cell, applicable when `method='knn'`.

        imageid (str, optional):  
            The column name in `adata` containing identifiers for different images, allowing for analysis within specific images.

        subset (str, optional):  
            The identifier of a specific image to restrict the analysis to. If provided, analysis will only be performed on this subset.

        label (str, optional):  
            The key under which to store the results in `adata.obs`, allowing for customized labeling of the output.
        
        verbose (bool):  
        If set to `True`, the function will print detailed messages about its progress and the steps being executed.

Returns:
        adata (anndata.AnnData):    
            The input AnnData object updated with the results stored under `adata.obs[label]`, where `label` is the specified output label.
            
Example:
    ```python
    
    # Analyze spatial aggregation using the radius method
    adata = sm.tl.spatial_aggregate(adata, x_coordinate='X_centroid', y_coordinate='Y_centroid',
                                    phenotype='phenotype', method='radius', radius=30, purity=60,
                                    imageid='imageid', subset=None, label='spatial_aggregate_radius')
    
    # Analyze spatial aggregation using the knn method
    adata = sm.tl.spatial_aggregate(adata, x_coordinate='X_centroid', y_coordinate='Y_centroid',
                                    phenotype='phenotype', method='knn', knn=10, purity=60,
                                    imageid='imageid', subset=None, label='spatial_aggregate_knn')
    
    # Subset analysis to a specific image using the radius method
    adata = sm.tl.spatial_aggregate(adata, x_coordinate='X_centroid', y_coordinate='Y_centroid',
                                    phenotype='phenotype', method='radius', radius=30, purity=60,
                                    imageid='imageid', subset='image_01', label='spatial_aggregate_image_01')
    
    ```
    """
    
    # Error statements
    #if purity < 51:
    #    raise ValueError('purity should be set to a value greater than 50')
        
    def spatial_aggregate_internal (adata_subset, x_coordinate,y_coordinate,z_coordinate,phenotype,purity,
                                    method,radius,knn,imageid,subset,label):
    

        # Create a DataFrame with the necessary inforamtion
        if z_coordinate is not None:
            if verbose:
                print("Including Z -axis")
            data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate], 'z': adata_subset.obs[z_coordinate], 'phenotype': adata_subset.obs[phenotype]})
        else:
            data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate], 'phenotype': adata_subset.obs[phenotype]})

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
            
            
# =============================================================================
#         if method == 'knn':
#             if verbose:
#                 print("Identifying the " + str(knn) + " nearest neighbours for every cell")
#             tree = BallTree(data[['x','y']], leaf_size= 2)
#             ind = tree.query(data[['x','y']], k=knn, return_distance= False)
#             neighbours = pd.DataFrame(ind.tolist(), index = data.index) # neighbour DF
#             neighbours.drop(0, axis=1, inplace=True) # Remove self neighbour
#         
#         # b) Local radius method
#         if method == 'radius':
#             if verbose:
#                 print("Identifying neighbours within " + str(radius) + " pixels of every cell")
#             kdt = BallTree(data[['x','y']], leaf_size= 2) 
#             ind = kdt.query_radius(data[['x','y']], r=radius, return_distance=False)
#             for i in range(0, len(ind)): ind[i] = np.delete(ind[i], np.argwhere(ind[i] == i))#remove self
#             neighbours = pd.DataFrame(ind.tolist(), index = data.index) # neighbour DF
# =============================================================================
            
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
        
        # Count the neighbours
        k = n.groupby(['neighbourhood','neighbour_phenotype']).size().unstack().fillna(0)
        k = k.div(k.sum(axis=1), axis=0)
        
        # Iteratte over all rows and find the column which passes the purity test
        #def col_name_mapper (row_data, purity):
        #    p = row_data[row_data >= purity/100]
        #    #phenotype_name = 'non-significant' if len(p.index) == 0 else p.index[0]
        #    phenotype_name = 'non-significant' if len(p.index) == 0 else p.idxmax()
        #    return phenotype_name
        # Apply the iteration function
        #aggregate_pheno = pd.DataFrame(k.apply(lambda x: col_name_mapper(row_data=x,purity=purity), axis=1))


        # Within the spatial_aggregate_internal function
        # Create an empty DataFrame to hold the results
        aggregate_pheno = pd.DataFrame(index=k.index, columns=[0])
        
        # Iterate over rows in DataFrame k
        for idx, row in k.iterrows():
            filtered_row = row[row >= purity / 100]  # Apply purity threshold
            if not filtered_row.empty:  # Check if the filtered row is not empty
                # If not empty, find the index of the maximum value
                max_idx = filtered_row.idxmax()
            else:
                # If empty, set to 'non-significant'
                max_idx = 'non-significant'
            # Store the result
            aggregate_pheno.at[idx, 0] = max_idx
        
        
        #aggregate_pheno = pd.DataFrame(k[k>=purity/100].idxmax(axis=1).fillna('non-significant'))
        aggregate_pheno.columns = ['spatial_aggregate']
        
        # Return 
        return aggregate_pheno
    

    # Subset a particular image if needed
    if subset is not None:
        adata_list = [adata[adata.obs[imageid] == subset]]
    else:
        adata_list = [adata[adata.obs[imageid] == i] for i in adata.obs[imageid].unique()]
    
    # Apply function to all images and create a master dataframe
    # Create lamda function 
    r_spatial_aggregate_internal = lambda x: spatial_aggregate_internal(adata_subset=x,
                                                                          x_coordinate=x_coordinate,
                                                                          y_coordinate=y_coordinate,
                                                                          z_coordinate=z_coordinate,
                                                                          phenotype=phenotype,
                                                                          method=method,
                                                                          radius=radius,knn=knn,
                                                                          imageid=imageid,subset=subset,
                                                                          purity=purity,
                                                                          label=label) 
    all_data = list(map(r_spatial_aggregate_internal, adata_list)) # Apply function 
    
    # Merge all the results into a single dataframe    
    result = []
    for i in range(len(all_data)):
        result.append(all_data[i])
    result = pd.concat(result, join='outer')  
    
    # Reindex the cells
    result = result.reindex(adata.obs.index)
    result = result.fillna('non-significant')
    
    # Add to adata
    adata.obs[label] = result
    
    # Return        
    return adata

