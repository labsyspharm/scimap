#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 15:00:39 2020
@author: Ajit Johnson Nirmal
Function to find regions of aggregration of similar cells
"""

# Import library
import pandas as pd
import numpy as np
from sklearn.neighbors import BallTree

# Function
def spatial_aggregate (adata, x_coordinate='X_centroid',y_coordinate='Y_centroid',
                       purity = 60, phenotype='phenotype', method='radius', radius=30, knn=10, 
                       imageid='imageid',subset=None,label='spatial_aggregate'):
    """
    
    Parameters
    ----------
    adata : AnnData object

    x_coordinate : float, required
        Column name containing the x-coordinates values. The default is 'X_centroid'.
    y_coordinate : float, required
        Column name containing the y-coordinates values. The default is 'Y_centroid'.
    purity : int, optional
        Supply a value between 1 to 100. It is the percent purity of neighbouring cells.
        For e.g. if 60 is chosen, every neighbourhood is tested such that if a 
        particular phenotype makes up greater than 60% of the total 
        population it is annotated to be an aggregate of that particular phenotype. The default is 60.
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
        Key for the returned data, stored in `adata.obs`. The default is 'spatial_aggregate'.

    Returns
    -------
    adata : AnnData object
        Updated AnnData object with the results stored in `adata.obs['spatial_aggregate']`.
        
        
    Example
    -------
    # Running the radius method
    adata = sm.tl.spatial_aggregate (adata, x_coordinate='X_centroid',y_coordinate='Y_centroid',
                        phenotype='phenotype', method='radius', radius=30, purity = 60,
                        imageid='imageid',subset=None,label='spatial_aggregate_radius')
    # Running the knn method
    adata =  sm.tl.spatial_aggregate (adata, x_coordinate='X_centroid',y_coordinate='Y_centroid',
                        phenotype='phenotype', method='knn', knn=10, purity = 60,
                        imageid='imageid',subset=None,label='spatial_aggregate_knn')

    """
    
    # Error statements
    #if purity < 51:
    #    raise ValueError('purity should be set to a value greater than 50')
        
    def spatial_aggregate_internal (adata_subset, x_coordinate,y_coordinate,phenotype,purity,
                                    method,radius,knn,imageid,subset,label):
    

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
            kdt = BallTree(data[['x','y']], leaf_size= 2) 
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
        aggregate_pheno = pd.DataFrame(k[k>=purity/100].idxmax(axis=1).fillna('non-significant'))
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
    result = result.fillna(0)
    result = result.reindex(adata.obs.index)
    
    # Add to adata
    adata.obs[label] = result
    
    # Return        
    return adata

