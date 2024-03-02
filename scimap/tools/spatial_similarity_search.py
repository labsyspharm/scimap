# -*- coding: utf-8 -*-
# Created on Thu May 12 10:01:20 2022
# @author: Ajit Johnson Nirmal

"""
!!! abstract "Short Description"
    `sm.tl.spatial_similarity_search`: This function enables the discovery of regions 
    within spatial datasets that resemble a specified region of interest (ROI) 
    in terms of molecular characteristics. Outcomes are stored in `adata.obs`, 
    facilitating subsequent visualization with `sm.pl.image_viewer`. 
    Users can iteratively adjust the `similarity_threshold` to refine results, 
    optimizing the identification of molecularly similar areas.

## Function
"""

# Import Libs
import pandas as pd
import numpy as np
from sklearn.neighbors import BallTree
import scipy
import anndata
import pathlib
import numba
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import RobustScaler


# Function

def spatial_similarity_search (adata,
                               ROI_column,
                               x_coordinate='X_centroid',
                               y_coordinate='Y_centroid',
                               z_coordinate=None,
                               similarity_threshold=0.5,
                               ROI_subset = None,
                               method='radius', 
                               radius=30, 
                               knn=10, 
                               imageid='imageid', 
                               use_raw=True, 
                               subset=None,
                               label='spatial_similarity_search',
                               reuse_similarity_matrix=None,
                               morphological_features=None,
                               use_only_morphological_features=False,
                               verbose=True,
                               output_dir=None):
    """
Parameters:
        adata (anndata.AnnData):  
            Annotated data matrix or path to an AnnData object, containing spatial gene expression data.
            
        x_coordinate (str, required):  
            Column name in `adata` for the x-coordinates.

        y_coordinate (str, required):  
            Column name in `adata` for the y-coordinates.

        z_coordinate (str, optional):  
            Column name in `adata` for the z-coordinates, for 3D spatial data analysis.
            
        ROI_column (str, required):  
            Column name in `adata.obs` indicating the ROIs or regions for similarity analysis.
            
        similarity_threshold (float, optional):  
            Threshold to adjust the strictness of similarity, with a default of 0.5.
            
        ROI_subset (list, optional):  
            List of specific ROIs within `ROI_column` for targeted similarity search.
                    
        method (str, optional):  
            'radius' or 'knn' for neighborhood definition, with a focus on proximity or nearest neighbors.
            
        radius, knn (int, optional):  
            Parameters defining the spatial neighborhood for the chosen method.
            
        imageid (str, optional):  
            Column name in `adata.obs` for image identifiers, facilitating analysis within specific images.
            
        use_raw (bool, optional):  
            Whether to utilize raw data for analysis.
            
        subset (str, optional):  
            Specific image identifier for targeted analysis, typically an image ID.
            
        label (str, optional):  
            Custom label for storing results in `adata.obs` and `adata.uns`.
            
        reuse_similarity_matrix (str, optional):  
            Reuse a previously computed similarity matrix to adjust thresholds without recalculating.
            
        morphological_features (list, optional):  
            List of morphological features stored in `adata.obs` to include in similarity calculations.
            
        use_only_morphological_features (bool, optional):  
            If true, calculates similarity based solely on morphological features.
        
        verbose (bool):  
            If set to `True`, the function will print detailed messages about its progress and the steps being executed.

Returns:
        adata (anndata.AnnData):  
            Updated `adata` object with spatial similarity search results stored under `adata.obs[{label}_{ROIname}]`.

Example:
        ```python
        
        # Basic spatial similarity search within a specific ROI using radius method
        adata = sm.tl.spatial_similarity_search(adata, ROI_column='Tumor_ROI', similarity_threshold=0.6,
                                          ROI_subset=['Tumor_ROI_1'], method='radius', radius=40,
                                          label='tumor_similarity')
    
        # Adjusting similarity threshold and reusing similarity matrix
        adata = sm.tl.spatial_similarity_search(adata, ROI_column='Tumor_ROI', similarity_threshold=0.8,
                                          reuse_similarity_matrix='tumor_similarity', label='tumor_similarity_adjusted')
    
        # Incorporating morphological features in similarity search
        adata = sm.tl.spatial_similarity_search(adata, ROI_column='Immune_ROI', similarity_threshold=0.5,
                                          morphological_features=['Area', 'MajorAxisLength'],
                                          use_only_morphological_features=True, label='immune_similarity_morph')
        
        # visulaize the results in napari
        image_path = "/Users/aj/Documents/exemplar-001/registration/exemplar-001.ome.tif"
        sm.pl.image_viewer (image_path, adata, subset = 'unmicst-exemplar-001_cell', 
                            overlay='spatial_similarity_search_blood_vessel', point_color='White')
        
        ```
    """
    
    
    #x_coordinate='X_centroid'; y_coordinate='Y_centroid'; method='radius'; radius=30; knn=10; imageid='imageid'; 
    #use_raw=True ; log=True; subset=None; label='spatial_similarity_search'; output_dir=None; ROI_column='ASMA'; ROI_subset = None; similarity_threshold=0.5
    # morphological_features = ['Area', 'MajorAxisLength','MinorAxisLength', 'Eccentricity', 'Solidity', 'Extent', 'Orientation']
    #adata_subset = adata.copy()
    
    # Load the andata object    
    if isinstance(adata, str):
        imid = str(adata.rsplit('/', 1)[-1])
        adata = anndata.read(adata)
    else:
        imid = "adata_spatial_similarity_search"
        adata = adata
 
    # Error statements
    if use_raw is False:
        if all(adata.X[0] < 1) is False:
            raise ValueError('Please run `sm.pp.rescale` first if you wish to use `use_raw = False`')
            
    # Function to calculate the distance between two vectors
    @numba.jit(nopython=True, parallel=True, cache=True)
    def euclidian_score(query_neighbourhood):
        return 1.0 / ((np.sqrt(np.sum((spatial_lag_array - query_neighbourhood) ** 2, axis=1))) + 1.0)
    
    
    def spatial_expression_internal (adata_subset, x_coordinate, y_coordinate,z_coordinate,
                                     method, radius, knn, imageid, use_raw,
                                     morphological_features, use_only_morphological_features):
        
        # Create a DataFrame with the necessary inforamtion
        #data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate]})
        
        # Create a dataFrame with the necessary inforamtion
        if z_coordinate is not None:
            if verbose:
                print("Including Z -axis")
            data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate], 'z': adata_subset.obs[z_coordinate] })
        else:
            data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate] })

        
        # Identify neighbourhoods based on the method used
        # a) KNN method
        if method == 'knn':
            if verbose:
                print("Identifying the " + str(knn) + " nearest neighbours for every cell")
            if z_coordinate is not None:
                tree = BallTree(data, leaf_size= 2)
                dist, ind = tree.query(data, k=knn, return_distance= True)
            else:
                tree = BallTree(data, leaf_size= 2)
                dist, ind = tree.query(data, k=knn, return_distance= True)

        # b) Local radius method
        if method == 'radius':
            if verbose:
                print("Identifying neighbours within " + str(radius) + " pixels of every cell")
            if z_coordinate is not None:
                kdt = BallTree(data, metric='euclidean') 
                ind = kdt.query_radius(data, r=radius, return_distance=True)
            else:
                kdt = BallTree(data, metric='euclidean') 
                ind = kdt.query_radius(data, r=radius, return_distance=True)

# =============================================================================
#         
#         if method == 'knn':
#             print("Identifying the " + str(knn) + " nearest neighbours for every cell")
#             tree = BallTree(data, leaf_size= 2)
#             dist, ind = tree.query(data, k=knn, return_distance= True)
# 
#             
#         # b) Local radius method
#         if method == 'radius':
#             print("Identifying neighbours within " + str(radius) + " pixels of every cell")
#             kdt = BallTree(data, metric='euclidean')
#             ind, dist = kdt.query_radius(data, r=radius, return_distance= True)
#             
# =============================================================================
            
            
        # Normalize range (0-1) and account for total number of cells 
        d = scipy.sparse.lil_matrix((len(data), len(data)))
        for row, (columns, values) in enumerate(zip(ind, dist)):
            # Drop self-distance element.
            idx = columns != row
            columns = columns[idx]
            values = values[idx]
            if len(values) == 1:
                values = [1.0]
            elif len(values) > 1:
                # Normalize distances.
                values = (values.max() - values) / (values.max() - values.min())
                values /= values.sum()
            # Assign row to matrix.
            d[row, columns] = values
        
        # convert to csr sparse matrix
        wn_matrix_sparse = d.tocsr()
        
        # to dense matrix
        #dense = pd.DataFrame(wn_matrix_sparse.todense())
        
        
        # normalize data
        molecular_matrix = pd.DataFrame(adata_subset.raw.X, columns = adata_subset.var.index, index=adata_subset.obs.index)
        # clip outliers
        def clipping (x):
                clip = x.clip(lower =np.percentile(x,0.01), upper=np.percentile(x,99.99)).tolist()
                return clip
        gmm_data = molecular_matrix.apply(clipping)
        # log trasform
        normalised_data = np.log1p(gmm_data)
        # scale data
        #transformer = RobustScaler().fit(n_log)
        #normalised_data = pd.DataFrame(transformer.transform(n_log), columns = adata_subset.var.index, index=adata_subset.obs.index)
        #normalised_data = n_log

        #### Calculation of spatial lag
        
        # a) use only morphological features?
        if morphological_features is not None:
            if isinstance(morphological_features, str):
                morphological_features = [morphological_features]
            morph_f = adata_subset.obs[morphological_features]
            
            transformer = RobustScaler().fit(morph_f)
            morph_f = pd.DataFrame(transformer.transform(morph_f), columns = morph_f.columns, index=morph_f.index)
            
            if use_only_morphological_features is True:
                spatial_lag = pd.DataFrame(wn_matrix_sparse * morph_f, columns = morph_f.columns, index=morph_f.index)
        
        
        # b) use morphological features and molecular features?
        if morphological_features is not None and use_only_morphological_features is False:
            if use_raw==True:
                combined_matrix = pd.concat([normalised_data, morph_f], axis=1)
                spatial_lag = pd.DataFrame(wn_matrix_sparse * combined_matrix, columns = combined_matrix.columns, index=combined_matrix.index)     
            else:
                molecular_matrix = pd.DataFrame(adata_subset.X, columns = adata_subset.var.index, index=adata_subset.obs.index)
                combined_matrix = pd.concat([molecular_matrix, morph_f], axis=1)
                spatial_lag = pd.DataFrame(wn_matrix_sparse * combined_matrix, columns = combined_matrix.columns, index=combined_matrix.index) 
                    
        # c) use only molecular features
        if morphological_features is None:
            if use_raw==True:
                spatial_lag = pd.DataFrame(wn_matrix_sparse * normalised_data, columns = adata_subset.var.index, index=adata_subset.obs.index)
            else:
                spatial_lag = pd.DataFrame(wn_matrix_sparse * adata_subset.X, columns = adata_subset.var.index, index=adata_subset.obs.index)
                      
        # return value
        return spatial_lag


    # check if the user wants to reuse the spatial lag vector that was previously calculated
    if reuse_similarity_matrix is None:
        # Subset a particular image if needed
        if subset is not None:
            if isinstance(subset, str):
                subset = [subset]
            adata_list = [adata[adata.obs[imageid].isin(subset)]]
        else:
            adata_list = [adata[adata.obs[imageid] == i] for i in adata.obs[imageid].unique()]
    
    
        # Apply function to all images and create a master dataframe
        # Create lamda function 
        r_spatial_expression_internal = lambda x: spatial_expression_internal(adata_subset=x, 
                                                                    x_coordinate=x_coordinate, 
                                                                    y_coordinate=y_coordinate, 
                                                                    z_coordinate=z_coordinate,
                                                                    method=method, radius=radius, 
                                                                    knn=knn, imageid=imageid, 
                                                                    use_raw=use_raw,
                                                                    morphological_features=morphological_features, 
                                                                    use_only_morphological_features=use_only_morphological_features) 
        
        all_data = list(map(r_spatial_expression_internal, adata_list)) # Apply function 
    
        # Merge all the results into a single dataframe    
        result = []
        for i in range(len(all_data)):
            result.append(all_data[i])
        result = pd.concat(result, join='outer')  
        
        # save the results in adata object
        adata.uns[label] = result
    else:
        result = adata.uns[reuse_similarity_matrix]


    ### Identify the ROI's of interest and then correlate it with all spatial lag
    
    # calculate the distance between queri ROI and all neighbourhoods
    # figure out the roi's that need to be processed
    if ROI_subset is None:
        ROI_subset = list(adata.obs[ROI_column].unique())
        ROI_subset = [ x for x in ROI_subset if x != 'Other' ]
    else:
        if isinstance(ROI_subset, str):
            ROI_subset = [ROI_subset]
            

    # Figure out all the cells that fall within the user defined ROI's
    query_neigh = adata[adata.obs[ROI_column].isin(ROI_subset)].obs[[ROI_column]]
        
    #result = spatial_lag
    #np.max(all_roi_scores)
    #np.min(all_roi_scores)
    
    # for each ROI calculate the median spatial lag
    median_spatial_lag = pd.merge(result.loc[query_neigh.index], query_neigh, left_index=True, right_index=True, how='outer')
    median_spatial_lag = median_spatial_lag.groupby(ROI_column).median()
        
    # apply the distance function to each defined ROI's
    spatial_lag_array = np.array(result)
    median_spatial_lag_array = np.array(median_spatial_lag)
    
    # func
    #all_roi_scores = np.array([distance.euclidean(median_spatial_lag_array[0],x) for x in spatial_lag_array])
    #all_roi_scores = pd.DataFrame(all_roi_scores, columns=median_spatial_lag.index, index = result.index)
    all_roi_scores = np.array([euclidian_score(x) for x in median_spatial_lag_array])
    #all_roi_scores = median_spatial_lag.apply(euclidian_score, axis = 1) when spatial lag is a df
    # convert that to a df for returning
    all_roi_scores = pd.DataFrame(all_roi_scores, index=median_spatial_lag.index, columns = result.index).T
    
    # rescale the scores
    scaler = MinMaxScaler(feature_range=(0, 1))
    s = scaler.fit_transform(all_roi_scores)
    all_roi_scores = pd.DataFrame(s, columns = all_roi_scores.columns, index= all_roi_scores.index)
        
    ### Threshold the results to identify neighbourhoods that are similar to the 
    all_roi_scores_threshold = pd.DataFrame(np.where(all_roi_scores >= similarity_threshold, 'similar_to_ROI', 'other'), index = all_roi_scores.index, columns = all_roi_scores.columns)
    
    # rename columns of the 
    A = list(all_roi_scores.columns)
    column_names = [label + "_" + str(s) for s in A]
    all_roi_scores_threshold.columns = column_names
    
    # delete the result columns from adata if they already exist
    adata_obs = adata.obs
    if any(x in adata_obs.columns for x in column_names):
        adata_obs = adata_obs.drop(column_names, axis = 1)
    
    # Merge the results with adata.obs
    final_results = pd.merge(adata_obs, all_roi_scores_threshold, left_index=True, right_index=True, how='outer')
    
    # Reindex the cells
    final_results = final_results.replace(np.nan, 'not_computed')
    final_results = final_results.reindex(adata.obs.index)
    
    # return the data
    adata.obs = final_results
    
    # Save data if requested
    if output_dir is not None:
        output_dir = pathlib.Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        adata.write(output_dir / imid)
    else:    
        # Return data
        return adata


        
