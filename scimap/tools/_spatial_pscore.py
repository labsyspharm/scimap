#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 22:55:06 2021
@author: Ajit Johnson Nirmal
A system to score user defined interation between cell types.
The function generates two scores and saved at `adata.uns`: 
A) Proximity Density: Total number of interactions identified divided by the total number of 
cells of the cell-types that were used for interaction analysis.
B) Proximity Volume: Total number of interactions identified divided by the total number of all cells in the data.
The interaction sites are also recorded and saved in `adata.obs`
"""

# Import library
import pandas as pd
from sklearn.neighbors import BallTree
import numpy as np

# Function
def spatial_pscore (adata,proximity, score_by='imageid', x_coordinate='X_centroid',y_coordinate='Y_centroid',
                    phenotype='phenotype',method='radius',radius=20,knn=3,
                    imageid='imageid',subset=None, label='spatial_pscore'):
    """
    

    Parameters
    ----------
    adata : AnnData object

    proximity : list
        Pass a list of cell-types for which the proximity score needs to calculated. e.g. ['CellType-A', 'CellType-B']
    score_by : string, optional
        If the scores need to compared across region's of interest, the column name containing the ROI's
        should be passed. By default the score is calculated across the entire image. The default is 'imageid'.
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
        The radius used to define a local neighbhourhood. The default is 20.
    knn : int, optional
        Number of cells considered for defining the local neighbhourhood. The default is 3.
    imageid : string, optional
        Column name of the column containing the image id. The default is 'imageid'.
    subset : string, optional
        imageid of a single image to be subsetted for analyis. The default is None.
    label : string, optional
        Key for the returned data, stored in `adata.obs` and `adata.uns`. The default is 'spatial_pscore'.

    Returns
    -------
    adata : AnnData object
        Updated AnnData object with the results stored in `adata.obs['spatial_pscore']` and `adata.uns['spatial_pscore']`.

    Example
    -------
    # Calculate the score for proximity between `Tumor CD30+` cells and `M2 Macrophages`
    adata =  sm.tl.spatial_pscore (adata,proximity= ['Tumor CD30+', 'M2 Macrophages'], score_by = 'ImageId',
                             x_coordinate='X_position',y_coordinate='Y_position',
                             phenotype='phenotype',method='radius',radius=20,knn=3,
                             imageid='ImageId',subset=None, label='spatial_pscore')
    

    """
    
    
    # Start
    def spatial_pscore_internal (adata_subset,proximity,x_coordinate,y_coordinate,phenotype,method,radius,knn,
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
            neighbours_ind = neighbours.copy() # neighbour DF
            #neighbours.drop(0, axis=1, inplace=True) # Remove self neighbour
        
        # b) Local radius method
        if method == 'radius':
            print("Identifying neighbours within " + str(radius) + " pixels of every cell")
            kdt = BallTree(data[['x','y']], metric='euclidean') 
            ind = kdt.query_radius(data[['x','y']], r=radius, return_distance=False)
            #for i in range(0, len(ind)): ind[i] = np.delete(ind[i], np.argwhere(ind[i] == i))#remove self
            neighbours = pd.DataFrame(ind.tolist(), index = data.index) # neighbour DF
            neighbours_ind = neighbours.copy() # neighbour DF
            
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
        for i in proximity:
            print (str('Finding neighbourhoods with ') + str(i))
            nn = neighbours[neighbours.isin([i])].dropna(how='all').index
            neighbours = neighbours.loc[nn]
        
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
                                                   y_coordinate=y_coordinate,phenotype=phenotype,
                                                   method=method,radius=radius,knn=knn,
                                                   imageid=imageid,subset=subset,label=label) 
    all_data = list(map(r_spatial_pscore_internal, adata_list)) # Apply function 
    
    
    # Merge all the results into a single dataframe    
    proximity_site_cells =  np.concatenate([d['image_neighbours'] for d in all_data], axis=0)
    
    # Add it to the AnnData Object
    adata.obs[label] = np.where(adata.obs.index.isin(proximity_site_cells), '_'.join(proximity), "other")
    
    
    ##### SCORING #####
    proximity_neigh = np.concatenate([d['neighbours'] for d in all_data], axis=0)
    wh_d = adata.obs
    wh_d[label] = np.where(wh_d.index.isin(proximity_neigh), '_'.join(proximity), "other")
    
    # Define a scoring system
    name = '_'.join(proximity)
    whole_data = wh_d[[score_by, label, phenotype]]
    
    # proximity volume
    p_v = whole_data.groupby([score_by, label]).size().unstack().fillna(0)
    p_v ['All Cells'] = p_v[name] + p_v["other"]
    p_v['Proximity Volume'] = p_v[name] / p_v['All Cells']
    p_v = p_v.fillna(0) # replace NA
    p_v = p_v.replace([np.inf, -np.inf], 0)
    p_v = p_v.drop(columns = 'other')
    
    # subset the phenotypes of interest
    w_d = whole_data[whole_data[phenotype].isin(proximity)]
    # proximity density
    p_d = w_d.groupby([score_by, label]).size().unstack().fillna(0)
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
    print("Please check:\nadata.obs['" + str(label) + "'] &\nadata.uns['"+ str(label) + "'] for results")
      
    # Return 
    return adata

         