#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Mon Nov 16 08:34:04 2020
# @author: Ajit Johnson Nirmal and Yuan Chen
"""
!!! abstract "Short Description"
    `sm.hl.add_roi_omero`: This function integrates Regions of Interest (ROIs)
    extracted from Omero into AnnData objects, enriching spatial datasets with
    precise, user-defined annotations. An optional buffering parameter allows
    the creation of a ring-like border region for specified ROIs to aid in spatial analyses.
      
    The function allows users to add annotations extracted from Omero via:
    https://gist.github.com/Yu-AnChen/58754f960ccd540e307ed991bc6901b0.
    
## Function
"""

# Library imports
import pandas as pd
import numpy as np
import re
import matplotlib.patches as mpatches
import scipy.spatial.distance as sdistance
from joblib import Parallel, delayed

# Importing shapely for geometric operations (buffering, union, difference)
from shapely.geometry import Polygon, Point
from shapely.ops import unary_union

def addROI_omero(adata, 
                 roi, 
                 x_coordinate='X_centroid',
                 y_coordinate='Y_centroid',
                 imageid='imageid', 
                 naming_column='Name',
                 subset=None, 
                 overwrite=True,
                 label='ROI',
                 buffer_roi=0,  # Buffer distance in pixels for border creation (e.g., 40)
                 buffer_regions=None,  # If None, apply to all; else a str or list specifying ROI names for buffering
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
        
        naming_column (str):  
            The column name in the `roi` file that contains the ROI names—typically, this is `Name` or `Text`.

        subset (list, optional):  
            Specific image identifier(s) for targeted ROI addition.

        overwrite (bool):  
            If True, replaces existing ROI data; if False, appends new ROI data without overwriting.

        label (str):  
            Label under which ROI data will be stored in `adata.obs`.

        buffer_roi (int or float, optional):  
            A numeric value (in pixels) specifying the total buffering distance. When positive,
            a border ROI is computed by contracting the ROI inward by half the buffer and expanding it
            outward by the remaining half. The new annotation will be named "[ROI_name]_border". The function
            ensures that the buffered borders do not overlap other ROIs by enforcing a minimum 1-pixel gap.

        buffer_regions (None, str, or list, optional):  
            If set to None (default), the buffering is applied to all ROIs. Otherwise, it may be provided as a string or list of strings
            specifying the ROI name(s) (as found in the `naming_column`) to which the buffering should be applied.

        n_jobs (int):  
            Number of jobs for parallel processing. Defaults to -1, using all cores.

        verbose (bool):  
            If set to True, the function will print detailed progress messages and statistics.

Returns:
        adata (anndata.AnnData):  
            The updated AnnData object with ROI annotations added to adata.obs[label].

Example:
        ```python
        # Load the ROI CSV file
        roi_df = pd.read_csv('path/to/roi.csv')
        
        # Apply composite labeling for overlapping ROIs and buffered borders to all ROIs with a total buffer of 40 pixels
        adata = sm.hl.addROI_omero(adata, roi=roi_df, label='Sample_ROI', buffer_roi=40, verbose=True)
    
        # Apply buffering only to ROI named 'Tumor 1' with a 40-pixel buffer
        adata = sm.hl.addROI_omero(adata, roi=roi_df, label='Sample_ROI', buffer_roi=40, buffer_regions='Tumor 1', verbose=True)
    
        # Apply buffering only to a list of ROIs
        adata = sm.hl.addROI_omero(adata, roi=roi_df, label='Sample_ROI', buffer_roi=40, buffer_regions=['Tumor 1', 'TLS'], verbose=True)
        ```
    """
    # Create a DataFrame with spatial coordinates and image identifiers.
    data = pd.DataFrame(adata.obs)[[x_coordinate, y_coordinate, imageid]]
    if verbose:
        print(f"[INFO] Total number of cells in dataset: {data.shape[0]}")

    # Subset data if a specific image subset is provided.
    if subset is not None:
        if isinstance(subset, str): 
            subset = [subset]
        sub_data = data[data['imageid'].isin(subset)]
        if verbose:
            print(f"[INFO] Using subset {subset}: {sub_data.shape[0]} cells selected.")
    else:
        sub_data = data

    # Helper function to parse coordinate strings from the ROI CSV.
    def parse_roi_points(all_points):
        return np.array(
            re.findall(r'\d+\.?\d+', all_points), dtype=float
        ).reshape(-1, 2)

    # Function to create an ellipse patch from four vertex points.
    def ellipse_points_to_patch(vertex_1, vertex_2, co_vertex_1, co_vertex_2):
        v_and_co_v = np.array([vertex_1, vertex_2, co_vertex_1, co_vertex_2])
        centers = v_and_co_v.mean(axis=0)
        d = sdistance.cdist(v_and_co_v, v_and_co_v, metric='euclidean')
        width = d[0, 1]
        height = d[2, 3]
        vector_2 = v_and_co_v[1] - v_and_co_v[0]
        vector_2 /= np.linalg.norm(vector_2)
        angle = np.degrees(np.arccos([1, 0] @ vector_2))
        ellipse_patch = mpatches.Ellipse(centers, width=width, height=height, angle=angle)
        return ellipse_patch

    # Function to retrieve a matplotlib patch corresponding to the ROI type.
    def get_mpatch(roi_row):
        points = parse_roi_points(roi_row['all_points'])
        roi_type = roi_row['type']
        if roi_type in ['Point', 'Line']:
            roi_mpatch = mpatches.Polygon(points, closed=False)
        elif roi_type in ['Rectangle', 'Polygon', 'Polyline']:
            roi_mpatch = mpatches.Polygon(points, closed=True)
        elif roi_type == 'Ellipse':
            roi_mpatch = ellipse_points_to_patch(*points)
        else:
            raise ValueError("Unrecognized ROI type: {}".format(roi_type))
        return roi_mpatch

    # Build a dictionary mapping each ROI Id to its corresponding shapely Polygon.
    main_polys = {}
    unique_ids = roi['Id'].unique()
    if verbose:
        print(f"[INFO] Found {len(unique_ids)} unique ROIs to process.")
    for roi_id in unique_ids:
        roi_row = roi[roi['Id'] == roi_id].iloc[0]
        points = parse_roi_points(roi_row['all_points'])
        try:
            poly = Polygon(points)
        except Exception as e:
            raise ValueError("Error creating polygon for ROI {}: {}".format(roi_row[naming_column], e))
        main_polys[roi_id] = poly

    # Main ROI assignment: assign each spatial coordinate to its corresponding ROI.
    def add_roi_internal(roi_id):
        roi_subset = roi[roi['Id'] == roi_id].iloc[0]
        roi_mpatch = get_mpatch(roi_subset)
        inside = sub_data[roi_mpatch.contains_points(sub_data[[x_coordinate, y_coordinate]])].copy()
        inside['ROI_internal'] = roi_subset[naming_column]
        return inside

    if verbose:
        print("[INFO] Assigning cells to main ROIs...")
    main_roi_list = Parallel(n_jobs=n_jobs, verbose=verbose)(
        delayed(add_roi_internal)(roi_id=i) for i in unique_ids
    )
    final_roi = pd.concat(main_roi_list)[['ROI_internal']]
    # Aggregate overlapping main ROI assignments into a composite label.
    final_roi = final_roi.groupby(final_roi.index).agg(
        {'ROI_internal': lambda x: '_'.join(sorted(set(x)))}
    )
    if verbose:
        print(f"[INFO] Main ROI assignment complete: {final_roi.shape[0]} cells assigned.")

    # If buffering is requested, compute border ROIs for eligible regions.
    if buffer_roi > 0:
        if verbose:
            print("[INFO] Buffering enabled. Computing border ROIs...")
        half_buffer = buffer_roi / 2.0

        def add_border_roi(roi_row):
            # Determine whether this ROI should be buffered.
            if buffer_regions is not None:
                allowed = buffer_regions if isinstance(buffer_regions, list) else [buffer_regions]
                if roi_row[naming_column] not in allowed:
                    return pd.DataFrame(columns=['ROI_internal'])
            roi_id = roi_row['Id']
            poly = main_polys[roi_id]
            try:
                inner = poly.buffer(-half_buffer)
            except Exception:
                inner = poly  # If inward buffering fails, revert to original polygon.
            outer = poly.buffer(half_buffer)
            border_poly = outer.difference(inner)
            # Enforce a minimum 1-pixel gap from other ROIs.
            other_polys = [main_polys[other_id].buffer(1) for other_id in main_polys if other_id != roi_id]
            if other_polys:
                union_other = unary_union(other_polys)
                border_poly = border_poly.difference(union_other)
            border_indices = []
            for idx, row in sub_data.iterrows():
                pt = Point(row[x_coordinate], row[y_coordinate])
                if border_poly.contains(pt):
                    border_indices.append(idx)
            df_border = pd.DataFrame(index=border_indices)
            df_border['ROI_internal'] = str(roi_row[naming_column]) + '_border'
            return df_border

        border_roi_list = Parallel(n_jobs=n_jobs, verbose=verbose)(
            delayed(add_border_roi)(roi_row) for idx, roi_row in roi.drop_duplicates('Id').iterrows()
        )
        border_df = pd.concat(border_roi_list)
        if verbose:
            print(f"[INFO] Border ROI assignment complete: {border_df.shape[0]} cells assigned to border regions.")
        # Exclude cells already assigned a main ROI.
        border_df = border_df[~border_df.index.isin(final_roi.index)]
        # Aggregate overlapping border ROI assignments.
        border_df = border_df.groupby(border_df.index).agg(
            {'ROI_internal': lambda x: '_'.join(sorted(set(x)))}
        )
        result = pd.merge(data, final_roi, left_index=True, right_index=True, how='outer')
        result = result.merge(border_df, left_index=True, right_index=True, how='outer', suffixes=('', '_border'))
        # Prefer main ROI; if missing, fill with border ROI.
        result['ROI_internal'] = result['ROI_internal'].fillna(result['ROI_internal_border'])
        result = result.reindex(adata.obs.index)
    else:
        result = pd.merge(data, final_roi, left_index=True, right_index=True, how='outer')
        result = result.reindex(adata.obs.index)
        result['ROI_internal'] = result['ROI_internal'].fillna('Other')

    # Handle the case where the label already exists in adata.obs.
    if label in adata.obs.columns:
        if not overwrite:
            old_roi = adata.obs[label]
            combined_roi = pd.merge(result, old_roi, left_index=True, right_index=True, how='outer')
            combined_roi['ROI_internal'] = combined_roi['ROI_internal'].fillna(combined_roi[label])
        else:
            combined_roi = result.copy()
            combined_roi['ROI_internal'] = combined_roi['ROI_internal'].fillna('Other')
    else:
        combined_roi = result.copy()
        combined_roi['ROI_internal'] = combined_roi['ROI_internal'].fillna('Other')

    adata.obs[label] = combined_roi['ROI_internal']

    # Generate summary statistics on ROI assignments.
    roi_counts = adata.obs[label].value_counts()
    if verbose:
        print("[INFO] ROI assignment summary:")
        print(roi_counts.to_string())
        print(f"[INFO] Total cells annotated: {adata.obs[label].notnull().sum()}")

    return adata
