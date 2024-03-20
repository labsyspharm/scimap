#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Thu Jan 28 22:55:06 2021
# @author: Ajit Johnson Nirmal

"""
!!! abstract "Short Description"
    sm.tl.spatial_pscore: This function introduces a refined scoring system to quantify the proximity between 
    specified cell types within spatial data, including support for 3D datasets. It calculates two distinct scores:
        
    - **Proximity Density**: Total number of interactions identified divided by the total number of 
    cells of the cell-types that were used for interaction analysis.  
    - **Proximity Volume**: Total number of interactions identified divided by the total number of all cells in the data.  
      
    Interaction sites are cataloged and accessible in `adata.obs`. Both scores are stored in `adata.uns`.

## Functions
"""


# Import library
import pandas as pd
from sklearn.neighbors import BallTree
import numpy as np

# Function
def spatial_pscore (adata,proximity, 
                    score_by='imageid', 
                    x_coordinate='X_centroid',
                    y_coordinate='Y_centroid',
                    z_coordinate= None,
                    phenotype='phenotype',
                    method='radius',
                    radius=20,
                    knn=3,
                    imageid='imageid',
                    subset=None, 
                    verbose= True,
                    label='spatial_pscore'):
    """
Parameters:
        adata (anndata.AnnData):  
            Annotated data matrix with spatial information.

        proximity (list):  
            List of cell types to calculate proximity scores for, e.g., ['CellType-A', 'CellType-B'].

        score_by (str, optional):  
            Column name for ROI comparison. Scores are computed within these regions if specified.

        x_coordinate (str, required):  
            Column name for x-coordinates.

        y_coordinate (str, required):  
            Column name for y-coordinates.

        z_coordinate (str, optional):  
            Column name for z-coordinates, for 3D spatial data.

        phenotype (str, required):  
            Column name indicating cell phenotype or classification.

        method (str, optional):  
            Neighborhood definition method: 'radius' for fixed distance, 'knn' for K nearest neighbors.

        radius (int, optional):  
            Radius defining local neighborhoods (applicable for 'radius' method).

        knn (int, optional):  
            Number of nearest neighbors to consider (applicable for 'knn' method).

        imageid (str, optional):  
            Column name specifying image identifiers, for analyses within specific images.

        subset (str, optional):  
            Identifier for subset analysis, typically an image ID.
        
        verbose (bool, optional):  
            If True, enables progress and informational messages.

        label (str, optional):  
            Custom label for storing results in `adata.obs` and `adata.uns`.

Returns:
        adata (anndata.AnnData):  
            Updated `adata` object with proximity scores stored in both `adata.obs[label]` and `adata.uns[label]`.

Example:
        ```python
        
        # Compute proximity scores between two cell types across all images
        adata = sm.tl.spatial_pscore(adata, proximity=['CellType-A', 'CellType-B'],
                               method='radius', radius=20, label='proximity_score_all')
    
        # Compute proximity scores within a specific image subset
        adata = sm.tl.spatial_pscore(adata, proximity=['CellType-C', 'CellType-D'],
                               method='knn', knn=3, imageid='imageid', subset='image_02',
                               label='proximity_score_image_02')
    
        # 3D data proximity score calculation
        adata = sm.tl.spatial_pscore(adata, proximity=['CellType-E', 'CellType-F'],
                               x_coordinate='X_centroid', y_coordinate='Y_centroid', z_coordinate='Z_centroid',
                               method='radius', radius=30, label='proximity_score_3D')
        
        ```
    """
    
    
    # Start
    def spatial_pscore_internal (adata_subset,proximity,x_coordinate,y_coordinate,z_coordinate,phenotype,method,radius,knn,
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
            neighbours_ind = neighbours.copy()

        if method == 'radius':
            if verbose:
                print("Identifying neighbours within " + str(radius) + " pixels of every cell")
            if z_coordinate is not None:
                kdt = BallTree(data[['x','y','z']], metric='euclidean') 
                ind = kdt.query_radius(data[['x','y','z']], r=radius, return_distance=False)
            else:
                kdt = BallTree(data[['x','y']], metric='euclidean') 
                ind = kdt.query_radius(data[['x','y']], r=radius, return_distance=False)
            neighbours = pd.DataFrame(ind.tolist(), index = data.index) # neighbour DF
            neighbours_ind = neighbours.copy() # neighbour DF
                            
# =============================================================================
#             
#         if method == 'knn':
#             print("Identifying the " + str(knn) + " nearest neighbours for every cell")
#             tree = BallTree(data[['x','y']], leaf_size= 2)
#             ind = tree.query(data[['x','y']], k=knn, return_distance= False)
#             neighbours = pd.DataFrame(ind.tolist(), index = data.index) # neighbour DF
#             neighbours_ind = neighbours.copy() # neighbour DF
#             #neighbours.drop(0, axis=1, inplace=True) # Remove self neighbour
#         
#         # b) Local radius method
#         if method == 'radius':
#             print("Identifying neighbours within " + str(radius) + " pixels of every cell")
#             kdt = BallTree(data[['x','y']], metric='euclidean') 
#             ind = kdt.query_radius(data[['x','y']], r=radius, return_distance=False)
#             #for i in range(0, len(ind)): ind[i] = np.delete(ind[i], np.argwhere(ind[i] == i))#remove self
#             neighbours = pd.DataFrame(ind.tolist(), index = data.index) # neighbour DF
#             neighbours_ind = neighbours.copy() # neighbour DF
#             
#             
# =============================================================================
        # Map phenotype
        phenomap = dict(zip(list(range(len(ind))), data['phenotype'])) # Used for mapping
        phenomap_ind = dict(zip(list(range(len(ind))), data.index)) # Used for mapping cell_nme
        
        # Loop through (all functionized methods were very slow)
        for i in neighbours.columns:
            neighbours[i] = neighbours[i].dropna().map(phenomap, na_action='ignore')
        # do the same index and cell name
        for i in neighbours_ind.columns:
            neighbours_ind[i] = neighbours_ind[i].dropna().map(phenomap_ind, na_action='ignore')
            
            
        # Idetify all the neighbourhoods that contains the user defined proximity phenotypes
        #for i in proximity:
        #    print (str('Finding neighbourhoods with ') + str(i))
        #    nn = neighbours[neighbours.isin([i])].dropna(how='all').index
        #    neighbours = neighbours.loc[nn]
        matches = np.ones(len(neighbours), bool)
        for v in proximity:
            matches &= (neighbours == v).any(axis=1)
        neighbours = neighbours[matches]
        

        # Identify all the cells that was part of the neighbourhood in this analysis
        neighbours_ind = neighbours_ind.loc[neighbours.index]
        neighbours_ind_unique = pd.unique(neighbours_ind.values.ravel())
        
        # subset the neighbourhood cells to include only the cells in the user defined list
        cleaned_neighbours_ind_unique = [x for x in neighbours_ind_unique if str(x) != 'nan']
        d = data.loc[cleaned_neighbours_ind_unique]
        d = d[d['phenotype'].isin(proximity)].index
        
        # return neighbours for score and image_neighbours for plotting on image
        return {'neighbours': neighbours.index, 'image_neighbours': d }
        
        
    # Subset a particular image if needed
    if subset is not None:
        adata_list = [adata[adata.obs[imageid] == subset]]
    else:
        adata_list = [adata[adata.obs[imageid] == i] for i in adata.obs[imageid].unique()]
    
    # Apply function to all images and create a master dataframe
    # Create lamda function 
    r_spatial_pscore_internal = lambda x: spatial_pscore_internal(adata_subset=x,proximity=proximity,
                                                   x_coordinate=x_coordinate,
                                                   y_coordinate=y_coordinate,
                                                   z_coordinate=z_coordinate,
                                                   phenotype=phenotype,
                                                   method=method,radius=radius,knn=knn,
                                                   imageid=imageid,subset=subset,label=label) 
    all_data = list(map(r_spatial_pscore_internal, adata_list)) # Apply function
    
    
    # Merge all the results into a single dataframe    
    proximity_site_cells =  np.concatenate([d['image_neighbours'] for d in all_data], axis=0)
    
    # Add it to the AnnData Object
    adata.obs[label] = np.where(adata.obs.index.isin(proximity_site_cells), '_'.join(proximity), "other")
        
    ##### SCORING #####
    proximity_neigh = np.concatenate([d['neighbours'] for d in all_data], axis=0)
    wh_d = adata.obs.copy()
    wh_d[label] = np.where(wh_d.index.isin(proximity_neigh), '_'.join(proximity), "other")
    
    # Define a scoring system
    name = '_'.join(proximity)
    whole_data = wh_d[[score_by, label, phenotype]]
    
    # proximity volume
    p_v = whole_data.groupby([score_by, label], observed=False).size().unstack().fillna(0)
    p_v ['All Cells'] = p_v[name] + p_v["other"]
    p_v['Proximity Volume'] = p_v[name] / p_v['All Cells']
    p_v = p_v.fillna(0) # replace NA
    p_v = p_v.replace([np.inf, -np.inf], 0)
    p_v = p_v.drop(columns = 'other')
    
    # subset the phenotypes of interest
    w_d = whole_data[whole_data[phenotype].isin(proximity)]
    # proximity density
    p_d = w_d.groupby([score_by, label], observed=False).size().unstack().fillna(0)
    p_d ['Celltype of interest'] = p_d[name] + p_d["other"]
    p_d['Proximity Density'] = p_d[name] / p_d['Celltype of interest']
    p_d = p_d.fillna(0) # replace NA
    p_d = p_d.replace([np.inf, -np.inf], 0)
    p_d = p_d.drop(columns = ['other', name])
    
    # Merge Promimity volumne and density
    proximity_score = pd.merge(p_v, p_d, left_index=True, right_index=True)
    
    # Add it to the anndata object
    adata.uns[label] = proximity_score
    
    # Print
    if verbose:
        print("Please check:\nadata.obs['" + str(label) + "'] &\nadata.uns['"+ str(label) + "'] for results")
      
    # Return 
    return adata

         