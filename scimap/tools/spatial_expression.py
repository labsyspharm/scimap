#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Created on Wed Aug 19 15:00:39 2020
# @author: Ajit Johnson Nirmal

"""
!!! abstract "Short Description"
    `sm.tl.spatial_expression`: This function generates a neighborhood weighted matrix from spatial data, 
    integrating expression values to assess local cellular environments. 
      
    It employs two approaches for neighborhood definition:

    - **Radius Method**: Identifies neighbors within a specified radius, enabling analyses based on physical proximity.
    - **KNN Method**: Determines neighbors based on the K nearest neighbors, focusing on immediate spatial relationships.
    
    The output, a proportion matrix reflecting local expression patterns, is stored in `adata.uns`. 
    This matrix can be further analyzed using the `spatial_cluster` function to identify Recurrent 
    Cellular Neighborhoods (RCNs), facilitating the exploration of spatial expression dynamics and 
    neighborhood-specific gene expression patterns.

## Function
"""

# Import library
import pandas as pd
import numpy as np
from sklearn.neighbors import BallTree
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
import scipy
import argparse
import sys
import anndata
import pathlib


def main(argv=sys.argv):
    parser = argparse.ArgumentParser(
        description='The function allows users to compute a neighbourhood weighted matrix based on the expression values.'
    )
    parser.add_argument(
        '--adata', required=True, 
        help='AnnData object loaded into memory or path to AnnData object.'
    )
    parser.add_argument(
        '--x_coordinate', type=str, required=False, default='X_centroid',
        help='Column name containing the x-coordinates values.'
    )
    parser.add_argument(
        '--y_coordinate', type=str, required=False, default='Y_centroid',
        help='Column name containing the y-coordinates values.'
    )
    parser.add_argument(
        '--z_coordinate', type=str, required=False, default=None,
        help='Column name containing the z-coordinates values.'
    )
    parser.add_argument(
        '--method', type=str, required=False, default='radius',
        help='Two options are available: a) `radius`, b) `knn`.'
    )
    parser.add_argument(
        '--radius', type=int, required=False, default=30,
        help='The radius used to define a local neighbhourhood.'
    )
    parser.add_argument(
        '--knn', type=int, required=False, default=10,
        help='Number of cells considered for defining the local neighbhourhood.'
    )
    parser.add_argument(
        '--imageid', type=str, required=False, default='imageid',
        help='Column name of the column containing the image id.'
    )
    parser.add_argument(
        '--use_raw', type=bool, required=False, default=True,
        help='Argument to denote whether to use the raw data or scaled data after applying `sm.pp.rescale`.'
    )
    parser.add_argument(
        '--log', type=bool, required=False, default=True,
        help='If `True`, the log of raw data is used. Set use_raw = `True` for this to take effect.'
    )
    parser.add_argument(
        '--subset', type=str, required=False, default=None,
        help='imageid of a single image to be subsetted for analyis.'
    )
    parser.add_argument(
        '--label', type=str, required=False, default='spatial_expression',
        help='Key for the returned data, stored in `adata.uns`.'
    )
    parser.add_argument(
        '--verbose', required=False, default=True,
        help='The function will print detailed messages about its progress.'
    )
    parser.add_argument(
        '--output_dir', type=str, required=False, default=None,
        help='Path to output directory.'
    )
    args = parser.parse_args(argv[1:])
    print(vars(args))
    spatial_expression(**vars(args))

# Function
def spatial_expression (adata, 
                        x_coordinate='X_centroid',
                        y_coordinate='Y_centroid',
                        z_coordinate=None,
                        method='radius', 
                        radius=30, 
                        knn=10, 
                        imageid='imageid', 
                        use_raw=True, 
                        log=True, 
                        subset=None,
                        label='spatial_expression',
                        verbose=True,
                        output_dir=None):
    """
Parameters:
        adata (anndata.AnnData):  
            The AnnData object, containing spatial gene expression data.

        x_coordinate (str, required):  
            Column name in `adata` for the x-coordinates.

        y_coordinate (str, required):  
            Column name in `adata` for the y-coordinates.
        
        z_coordinate (str, optional):  
            Column name in `adata` for the z-coordinates, for 3D spatial data analysis.
        
        method (str, optional):  
            Method for defining neighborhoods: 'radius' for fixed distance, 'knn' for K nearest neighbors.
        
        radius (int, optional):  
            Radius for defining local neighborhoods when using the 'radius' method.
        
        knn (int, optional):  
            Number of nearest neighbors to consider when using the 'knn' method.
        
        imageid (str, optional):  
            Column name in `adata` for image identifiers, useful for analyses within specific images.
        
        use_raw (bool, optional):  
            Whether to use raw or processed data for calculation. Log transformation is applied if `log=True`.
        
        log (bool, optional):  
            Apply log transformation to the data (requires `use_raw=True`).
        
        subset (str, optional):  
            Identifier for a subset of data, typically an image ID, for targeted analysis.
        
        label (str, optional):  
            Custom label for storing the weighted matrix in `adata.uns`.
        
        verbose (bool, optional):  
            If True, enables progress and informational messages.
        
        output_dir (str, optional):  
            Directory path for saving output files.

Returns:
        adata (anndata.AnnData):  
            The input `adata` object, updated with the spatial expression results in `adata.uns[label]`.

Example:
        ```python
        
        # Calculate spatial expression using a 30-pixel radius
        adata = sm.tl.spatial_expression(adata, x_coordinate='X_centroid', y_coordinate='Y_centroid',
                                   method='radius', radius=30, 
                                   label='expression_radius_30')
    
        # Calculate spatial expression using 10 nearest neighbors
        adata = sm.tl.spatial_expression(adata, x_coordinate='X_centroid', y_coordinate='Y_centroid',
                                   method='knn', knn=10, use_raw=True,
                                   label='expression_knn_10')
    
        # Analyze spatial expression within a specific image using radius method
        adata = sm.tl.spatial_expression(adata, x_coordinate='X_centroid', y_coordinate='Y_centroid',
                                   method='radius', radius=50, imageid='imageid', subset='specific_image',
                                   label='expression_specific_image')
        
        ```
    """
    
    # Load the andata object    
    if isinstance(adata, str):
        imid = str(adata.rsplit('/', 1)[-1])
        adata = anndata.read(adata)
    else:
        adata = adata

    
    # Error statements
    if use_raw is False:
        if all(adata.X[0] < 1) is False:
            raise ValueError('Please run `sm.pp.rescale` first if you wish to use `use_raw = False`')
        
     
    def spatial_expression_internal (adata_subset, x_coordinate, y_coordinate,z_coordinate,log,
                                     method, radius, knn, imageid, use_raw, subset,label):
         
        
        # Create a dataFrame with the necessary inforamtion
        if z_coordinate is not None:
            if verbose:
                print("Including Z -axis")
            data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate], 'z': adata_subset.obs[z_coordinate] })
        else:
            data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate] })


        # Create a DataFrame with the necessary inforamtion
        #data = pd.DataFrame({'x': adata_subset.obs[x_coordinate], 'y': adata_subset.obs[y_coordinate]})
        
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
                ind, dist = kdt.query_radius(data, r=radius, return_distance=True)
            else:
                kdt = BallTree(data, metric='euclidean') 
                ind, dist = kdt.query_radius(data, r=radius, return_distance=True)
        
# =============================================================================
#         # a) KNN method
#         if method == 'knn':
#             if verbose:
#                 print("Identifying the " + str(knn) + " nearest neighbours for every cell")
#             tree = BallTree(data, leaf_size= 2)
#             dist, ind = tree.query(data, k=knn, return_distance= True)
# 
#         # b) Local radius method
#         if method == 'radius':
#             if verbose:
#                 print("Identifying neighbours within " + str(radius) + " pixels of every cell")
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
        
        
        # Calculation of spatial lag
        if use_raw==True:
            if log is True:
                spatial_lag = pd.DataFrame(wn_matrix_sparse * np.log1p(adata_subset.raw.X), columns = adata_subset.var.index, index=adata_subset.obs.index)
            else:
                spatial_lag = pd.DataFrame(wn_matrix_sparse * adata_subset.raw.X, columns = adata_subset.var.index, index=adata_subset.obs.index)
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
                                                                z_coordinate=z_coordinate,
                                                                method=method, radius=radius, 
                                                                knn=knn, imageid=imageid, 
                                                                use_raw=use_raw, subset=subset,
                                                                log=log,
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
    
    # Save data if requested
    if output_dir is not None:
        output_dir = pathlib.Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        adata.write(output_dir / imid)
    else:    
        # Return data
        return adata


if __name__ == '__main__':
    main()
