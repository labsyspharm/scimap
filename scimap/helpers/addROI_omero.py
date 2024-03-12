#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Mon Nov 16 08:34:04 2020
# @author: Ajit Johnson Nirmal and Yuan Chen
"""
!!! abstract "Short Description"
    `sm.hl.add_roi_omero`: This function seamlessly integrates Regions of Interest (
        ROIs) extracted from Omero into AnnData objects, enriching spatial datasets 
    with precise, user-defined annotations for focused analysis.  
      
    The function allows users to add annotations that have been 
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


def addROI_omero (adata, 
                  roi, 
                  x_coordinate='X_centroid',
                  y_coordinate='Y_centroid',
                  imageid='imageid', 
                  subset=None, 
                  overwrite=True,
                  label='ROI', 
                  n_jobs=-1, 
                  verbose=False):
    
    """
Parameters:
        adata (anndata.AnnData):  
            Annotated data matrix or path to an AnnData object, containing spatial gene expression data.

        roi (DataFrame):  
            DataFrame containing ROI coordinates and identifiers, obtained from Omero.

        x_coordinate (str, required):  
            Column name in `adata` for the x-coordinates.

        y_coordinate (str, required):  
            Column name in `adata` for the y-coordinates.

        imageid (str):  
            Column name in `adata.obs` identifying distinct images in the dataset.

        subset (list, optional):  
            Specific image identifier(s) for targeted ROI addition.

        overwrite (bool):  
            If True, replaces existing ROI data; if False, appends new ROI data without overwriting.

        label (str):  
            Label under which ROI data will be stored in `adata.obs`.

        n_jobs (int):  
            Number of jobs for parallel processing. Defaults to -1, using all cores.

        verbose (bool):  
            If set to `True`, the function will print detailed messages about its progress and the steps being executed.

Returns:
        adata (anndata.AnnData):  
            The updated `adata` object with ROI annotations added to `adata.obs[label]`.

Example:
        ```python
        
        # load the ROI's
        roi_df = pd.read_csv('path/to/roi.csv')
        
        # Add ROIs to a single image dataset
        adata = sm.hl.addROI_omero(adata, roi=roi_df, label='Sample_ROI')
    
        # Add ROIs to a specific image in a dataset containing multiple images
        adata = sm.hl.addROI_omero(adata, roi=roi_df, imageid='image_column', subset=['image_01'], label='Image1_ROI')
    
        # Append multiple ROI datasets to the same column in AnnData without overwriting
        adata = sm.hl.addROI_omero(adata, roi=roi_df_first_set, label='Combined_ROI', overwrite=False)
        adata = sm.hl.addROI_omero(adata, roi=roi_df_second_set, label='Combined_ROI', overwrite=False)
        
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
