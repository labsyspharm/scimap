#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Mon Nov 16 08:34:04 2020
# @author: Ajit Johnson Nirmal and Yuan Chen
"""
!!! abstract "Short Description"
    `sm.hl.add_roi_omero`: The function allows users to add annotations that have been 
    extracted from Omero using the following 
    script: https://gist.github.com/Yu-AnChen/58754f960ccd540e307ed991bc6901b0.

## Function
"""

# Library
import pandas as pd
import numpy as np
import re
import matplotlib.patches as mpatches
import scipy.spatial.distance as sdistance
from joblib import Parallel, delayed
import argparse


def add_roi_omero (adata, roi, x_coordinate='X_centroid',y_coordinate='Y_centroid',
                   imageid='imageid', subset=None, overwrite=True,
                   label='ROI', n_jobs=-1, verbose=False):
    
    """
Parameters:

    adata : AnnData object

    roi : DataFrame  
        Pandas dataframe of ROI's that have been extracted from Omero using the following script: https://gist.github.com/Yu-AnChen/58754f960ccd540e307ed991bc6901b0.
        Please note that the function currently does not handle overlapping ROI's and so make sure the ROI's are mutually exclusive.

    x_coordinate (float):  
        Column name containing the x-coordinates values.

    y_coordinate (float):  
        Column name containing the y-coordinates values.
    
    imageid (string):  
        In the event that the adata object contains multiple images, it is
        important that ROIs are added to each image seperately. Pass the name 
        of the column that contains the `imageid` and use it in conjunction with
        the `subset` parameter to add ROI's to a specific image.

    subset (list):  
        Name of the image to which the ROI is to be added. Note if you have multiple images in the 
        adata object, you will need to add ROI's to each image one after the other independently. 
    
    overwrite (bool):  
        In the event you have multiple images in the adata object, ROI can be added to each image
        independently using the `imageid` and `subset` parameter. If you wish the results to be
        all saved with in the same column set this parameter to `False`. By default, the 
        function will overwrite the `label` column. 

    label (string):  
        Key for the returned data, stored in `adata.obs`.
        
    n_jobs (int):  
        Number of cores to use. Default is to use all available cores.

Returns:
    adata
        Modified AnnData object. Check `adata.obs` for an additional column.
    
Example:
```python
    roi = pd.read_csv('ROI/ROI_Z147_1_750.csv')
    adata = sm.hl.add_roi_omero (adata, roi, label='aj_ROI')
```
    """
    
    # create data matrix that has the co-ordinates
    data = pd.DataFrame(adata.obs)[[x_coordinate, y_coordinate, imageid]]
    
    # subset the data if needed
    if subset is not None:
        # convert string to list
        if isinstance(subset, str): 
            subset = [subset]
        # subset data
        sub_data = data[data['imageid'].isin(subset)]
    else:
        sub_data = data
    
    
    def parse_roi_points(all_points):
        return np.array(
            re.findall(r'-?\d+\.?\d+', all_points), dtype=float
        ).reshape(-1, 2)

    def ellipse_points_to_patch(
        vertex_1, vertex_2,
        co_vertex_1, co_vertex_2
    ):
        """
        Parameters
        ----------
        vertex_1, vertex_2, co_vertex_1, co_vertex_2: array like, in the form of (x-coordinate, y-coordinate)

        """
        v_and_co_v = np.array([
            vertex_1, vertex_2,
            co_vertex_1, co_vertex_2
        ])
        centers = v_and_co_v.mean(axis=0)

        d = sdistance.cdist(v_and_co_v, v_and_co_v, metric='euclidean')
        width = d[0, 1]
        height = d[2, 3]

        vector_2 = v_and_co_v[1] - v_and_co_v[0]
        vector_2 /= np.linalg.norm(vector_2)

        angle = np.degrees(np.arccos([1, 0] @ vector_2))

        ellipse_patch = mpatches.Ellipse(
            centers, width=width, height=height, angle=angle        
        )
        return ellipse_patch

    def get_mpatch(roi):
        points = parse_roi_points(roi['all_points'])

        roi_type = roi['type']
        if roi_type in ['Point', 'Line']:
            roi_mpatch = mpatches.Polygon(points, closed=False)
        elif roi_type in ['Rectangle', 'Polygon', 'Polyline']:
            roi_mpatch = mpatches.Polygon(points, closed=True)
        elif roi_type == 'Ellipse':
            roi_mpatch = ellipse_points_to_patch(*points)
        else:
            raise ValueError
        return roi_mpatch

    def add_roi_internal (roi_id):
        roi_subset = roi[roi['Id'] == roi_id].iloc[0]
        
        roi_mpatch = get_mpatch(roi_subset)
        inside = sub_data[roi_mpatch.contains_points(sub_data[[x_coordinate, y_coordinate]])]
        inside['ROI_internal'] = roi_subset['Name']

        # return
        return inside
    
    # Apply function to all rows of the ROI dataframe
    roi_list = roi['Id'].unique()
    final_roi = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(add_roi_internal)(roi_id=i) for i in roi_list)   
    
    # Merge all into a single DF
    final_roi = pd.concat(final_roi)[['ROI_internal']]
    
    # Add the list to obs
    result = pd.merge(data, final_roi, left_index=True, right_index=True, how='outer')
    
    # Reindex
    result = result.reindex(adata.obs.index)
    
    # check if adata already has a column with the supplied label
    # if avaialble overwrite or append depending on users choice
    if label in adata.obs.columns:
        if overwrite is False:
            # Append
            # retreive the ROI information
            old_roi = adata.obs[label]
            combined_roi = pd.merge(result, old_roi, left_index=True, right_index=True, how='outer')
            combined_roi['ROI_internal'] = combined_roi['ROI_internal'].fillna(combined_roi[label])
        else:
            # Over write
            combined_roi = result.copy()
            combined_roi['ROI_internal'] = combined_roi['ROI_internal'].fillna('Other')     
    else:
        # if label is not present just create a new one
        combined_roi = result.copy()
        combined_roi['ROI_internal'] = combined_roi['ROI_internal'].fillna('Other') 
    
    
    # Add to adata
    adata.obs[label] = combined_roi['ROI_internal']
    
    # return
    return adata

parser = argparse.ArgumentParser(description='Add ROI to AnnData object from Omero')

parser.add_argument('adata', help='Path to the AnnData object')
parser.add_argument('roi', help='Path to the ROI DataFrame')
parser.add_argument('--x_coordinate', default='X_centroid', help='Column name for x-coordinates')
parser.add_argument('--y_coordinate', default='Y_centroid', help='Column name for y-coordinates')
parser.add_argument('--imageid', default='imageid', help='Column name for image IDs')
parser.add_argument('--subset', help='List of image names to add ROIs')
parser.add_argument('--overwrite', action='store_true', help='Overwrite the label column')
parser.add_argument('--label', default='ROI', help='Key for the returned data')
parser.add_argument('--n_jobs', type=int, default=-1, help='Number of cores to use')
parser.add_argument('--verbose', action='store_true', help='Enable verbose mode')

args = parser.parse_args()

    # Call the add_roi_omero function
result = add_roi_omero(args.adata,args.roi, x_coordinate=args.x_coordinate, y_coordinate=args.y_coordinate,
                           imageid=args.imageid, subset=args.subset, overwrite=args.overwrite,
                           label=args.label, n_jobs=args.n_jobs, verbose=args.verbose)

    # Save the modified AnnData object
result.write('modified_anndata.h5ad')

if __name__ == '__main__':
    add_roi_omero_command_line()