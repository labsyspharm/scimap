#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 08:34:04 2020
@author: Ajit Johnson Nirmal
Omero ROI to Scimap
"""

# Library
import pandas as pd
import numpy as np
import re
import matplotlib.path as mpath
from joblib import Parallel, delayed


def add_roi (adata, roi, x_coordinate='X_centroid',y_coordinate='Y_centroid',label='ROI'):
    
    """
    
    Parameters
    ----------
    adata : AnnData object

    roi : DataFrame
        Pandas dataframe of ROI's that have been extracted from Omero using the following script: https://gist.github.com/Yu-AnChen/58754f960ccd540e307ed991bc6901b0.
        Please note that the function currently does not handle overlapping ROI's and so make sure the ROI's are mutually exclusive.
    x_coordinate : float, required
        Column name containing the x-coordinates values. The default is 'X_centroid'.
    y_coordinate : float, required
        Column name containing the y-coordinates values. The default is 'Y_centroid'.
    label : string, optional
        Key for the returned data, stored in `adata.obs`. The default is 'ROI'.

    Returns
    -------
    adata
        Modified AnnData object. Check `adata.obs` for an additional column.
    
    Example
    -------
    roi = pd.read_csv('ROI/ROI_Z147_1_750.csv')
    adata = sm.hl.add_roi (adata, roi, label='aj_ROI')

    """
    
    # create data matrix that has the co-ordinates
    data = pd.DataFrame(adata.obs)[[x_coordinate, y_coordinate]]
    
    def add_roi_internal (roi_id):
        roi_subset = roi[roi['Id'] == roi_id]
        
        # Make a list of all ROI's from the ROI table
        all_points = [
        np.array(re.findall(r'\d+\.\d+', s))
            .astype(np.float64)
            .reshape(-1, 2)
        for s in roi_subset['all_points']
        ]
        
        # Identify the cells within the ROI
        path = mpath.Path(all_points[0], closed=True)
        
        # the mask to query cells within one ROI
        inside = data[path.contains_points(data)]
        inside['ROI'] = roi_subset['Name'].iloc[0]

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
    
