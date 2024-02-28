#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Wed Aug 19 15:00:39 2020
# @author: Ajit Johnson Nirmal
"""
!!! abstract "Short Description"
    The `sm.tl.spatial_count` function allows users to compute a neighbourhood matrix 
    using any categorical variable (e.g. cell-types) as input.

    The function supports two methods to define a local neighbourhood <br>
    **Radius method**: Can be used to identifies the neighbours within a user defined radius for every cell.
    **KNN method**: Can be used to identifies the neighbours based on K nearest neigbours for every cell
        
    The resultant neighbourhood matrix is saved with `adata.uns`. 

    This can be further clustered to identify similar neighbourhoods. 
    Use the [spatial_cluster] function to further group the neighbourhoods into 
    Reccurent Cellular Neighbourhoods (RCNs)

## Function
"""

# Import library
import pandas as pd
import numpy as np
import argparse
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
                   label='spatial_count'):
    """
Parameters:
    adata : anndata object

    x_coordinate (float):  
        Column name containing the x-coordinates values.
        
    y_coordinate (float):   
        Column name containing the y-coordinates values.
        
    phenotype (string):   
        Column name of the column containing the phenotype information. 
        It could also be any categorical assignment given to single cells.
        
    method (string):   
        Two options are available: a) `radius`, b) `knn`.  
        a) radius - Identifies the neighbours within a given radius for every cell.  
        b) knn - Identifies the K nearest neigbours for every cell.  
        
    radius (int):   
        The radius used to define a local neighbhourhood.
        
    knn (int):   
        Number of cells considered for defining the local neighbhourhood.
        
    imageid (string):   
        Column name of the column containing the image id.
        
    subset (string):   
        imageid of a single image to be subsetted for analyis.
        
    label (string):   
        Key for the returned data, stored in `adata.uns`.

Returns:
    adata : anndata object  
        Updated AnnData object with the results stored in `adata.uns ['spatial_count']`.
    
Example:
    ```python
    # Running the radius method
    adata = sm.tl.spatial_count (adata,x_coordinate='X_centroid',
                                 y_coordinate='Y_centroid',
                                 phenotype='phenotype',
                                 method='radius',radius=30,
                                 imageid='imageid',subset=None,
                                 label='spatial_count_radius')
    
    # Running the knn method
    adata = sm.tl.spatial_count (adata,x_coordinate='X_centroid',
                                 y_coordinate='Y_centroid',
                                 phenotype='phenotype',method='knn',
                                 knn=10, imageid='imageid',
                                 subset=None,label='spatial_count_knn')
    ```
    """

    def spatial_count_internal (adata_subset,x_coordinate,y_coordinate,z_coordinate,phenotype,method,radius,knn,
                                imageid,subset,label):
        
        # Create a dataFrame with the necessary inforamtion
        if z_coordinate is not None:
            print("Including Z -axis")
            data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate], 'z': adata_subset.obs[z_coordinate], 'phenotype': adata_subset.obs[phenotype]})
        else:
            data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate], 'phenotype': adata_subset.obs[phenotype]})


        # Create a DataFrame with the necessary inforamtion
        #data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate], 'phenotype': adata_subset.obs[phenotype]})
        
        # Identify neighbourhoods based on the method used
        # a) KNN method
        if method == 'knn':
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
    result = result.fillna(0)
    result = result.reindex(adata.obs.index)
    
    # Add to adata
    adata.uns[label] = result
    
    # Return        
    return adata

if __name__ == '__main__':
    # Create argparse parser
    parser = argparse.ArgumentParser(description='Perform spatial counting analysis.')

    # Add arguments
    parser.add_argument('--adata',type=str ,help='Path to the AnnData object file.')
    parser.add_argument('--x_coordinate',type=str, default='X_centroid', help='Column name for x-coordinates.')
    parser.add_argument('--y_coordinate',type=str, default='Y_centroid', help='Column name for y-coordinates.')
    parser.add_argument('--phenotype',type=str, default='phenotype', help='Column name for phenotype information.')
    parser.add_argument('--method', type=str, default='radius', help='Method for identifying neighbors.')
    parser.add_argument('--radius', type=int, default=30, help='Radius used to define a local neighborhood.')
    parser.add_argument('--knn', type=int, default=10, help='Number of cells considered for defining the local neighborhood.')
    parser.add_argument('--imageid', type=str, default='imageid', help='Column name for image ID.')
    parser.add_argument('--subset', type=str, help='Image ID of a single image to be subsetted for analysis.')
    parser.add_argument('--label', type=str, default='spatial_count', help='Key for the returned data.')

    # Parse the command line arguments
    args = parser.parse_args()

    # Call the spatial_count function with the parsed arguments
    spatial_count(args.adata, args.x_coordinate, args.y_coordinate, args.phenotype,
                  args.method, args.radius, args.knn, args.imageid,
                  args.subset, args.label)