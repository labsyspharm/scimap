#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Mon Nov 16 08:34:04 2020
# @author: Ajit Johnson Nirmal
"""
!!! abstract "Short Description"
    `sm.hl.voronoi`: The function allows users to add annotations that have been 
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


def add_roi_omero (adata, roi, x_coordinate='X_centroid',y_coordinate='Y_centroid',label='ROI'):
    
    """
Parameters:

    adata : AnnData object

    roi : DataFrame  
        Pandas dataframe of ROI's that have been extracted from Omero using the following script: https://gist.github.com/Yu-AnChen/58754f960ccd540e307ed991bc6901b0.
        Please note that the function currently does not handle overlapping ROI's and so make sure the ROI's are mutually exclusive.

    x_coordinate : float, required  
        Column name containing the x-coordinates values.

    y_coordinate : float, required  
        Column name containing the y-coordinates values.

    label : string, optional  
        Key for the returned data, stored in `adata.obs`.

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
    data = pd.DataFrame(adata.obs)[[x_coordinate, y_coordinate]]
    
    def parse_roi_points(all_points):
        return np.array(
            re.findall(r'\d+\.?\d+', all_points), dtype=float
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
        inside = data[roi_mpatch.contains_points(data)]
        inside['ROI'] = roi_subset['Name']

        # return
        return inside
    
    # Apply function to all rows of the ROI dataframe
    roi_list = roi['Id'].unique()
    final_roi = Parallel(n_jobs=-1)(delayed(add_roi_internal)(roi_id=i) for i in roi_list)   
    
    # Merge all into a single DF
    final_roi = pd.concat(final_roi)[['ROI']]
    
    # Add the list to obs
    result = pd.merge(data, final_roi, left_index=True, right_index=True, how='outer')
    
    # Reindex
    result = result.reindex(adata.obs.index)
    result['ROI'] = result['ROI'].fillna('Other')
    
    # Add to adata
    adata.obs[label] = result['ROI']
    
    # return
    return adata
