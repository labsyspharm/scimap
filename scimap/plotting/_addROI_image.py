# -*- coding: utf-8 -*-
#Created on Mon Apr  4 11:12:59 2022
#@author: Ajit Johnson Nirmal
#Adding ROI to data via Napari.

"""
!!! abstract "Short Description"
    `sm.pl.addROI_image`: The function allows users to add ROIs to the data. 
    The function opens the image in a `napari` viewer and allows users to use the 
    shapes layer to add ROIs.   
    
    If the user is interesed in adding different ROI's (e.g. Tumor, Stroma, Tumor-Stromal-interface), 
    each of these should be drawn in a seperate `shape layer`. 
    
    Please note that a single `shape layer` can contain multiple seperate annotated regions. 
    All annotated regions within a single shape layer will be pooled into a single ROI. 
    
    The `shape layers` can be **renamed** to users choice of ROI name.  
    
    Please note: All ROI's have to be unique and cannot overlap.  
    
    *This function could also be used as a QC step to annotate regions of poor/ good quality 
    and later subsetted (keep or remove) for further analysis.*  
    
## Function
"""


# Library
#%gui qt
try:
    import napari
    from napari.layers import Shapes
except:
    pass
import pandas as pd
import random
import tifffile as tiff

import dask.array as da
from dask.cache import Cache
import zarr
import os
import matplotlib.patches as mpatches
import numpy as np
import scipy.spatial.distance as sdistance
from joblib import Parallel, delayed
from napari.utils.notifications import show_info


cache = Cache(2e9)  # Leverage two gigabytes of memory
cache.register()


# Function
def addROI_image (image_path, adata, subset=None,imageid='imageid', overlay=None, flip_y=True,
                    overlay_category=None,markers=None,channel_names='default',
                    x_coordinate='X_centroid',y_coordinate='Y_centroid',point_size=10,
                    point_color=None,seg_mask=None,
                    n_jobs=-1, verbose=False, 
                    overwrite=True, label='ROI', **kwargs):
    """
Parameters:
    image_path : string  
        Location to the image file (TIFF, OME.TIFF, ZARR supported)  

    adata : AnnData Object  
    
    subset : list, optional  
        Name of the image to which the ROI is to be added. Note if you have multiple images in the 
        adata object, you will need to add ROI's to each image one after the other independently.  
        
    imageid : string, optional  
        In the event that the adata object contains multiple images, it is
        important that ROIs are added to each image seperately. Pass the name 
        of the column that contains the `imageid` and use it in conjunction with
        the `subset` parameter to add ROI's to a specific image.

    seg_mask: string  
        Location to the segmentation mask file.  

    overlay : string, optional  
        Name of the column with any categorical data such as phenotypes or clusters.
        
    flip_y : bool, optional  
        Flip the Y-axis if needed. Some algorithms output the XY with the Y-coordinates flipped.
        If the image overlays do not align to the cells, try again by setting this to `False`.
        
    overlay_category : list, optional  
        If only specfic categories within the overlay column is needed, pass their names as a list.
        If None, all categories will be used.

    markers : list, optional  
        Markers to be included. If none, all markers will be displayed.

    channel_names : list, optional  
        List of channels in the image in the exact order as image. The default is `adata.uns['all_markers']`

    x_coordinate : string, optional  
        X axis coordinate column name in AnnData object.

    y_coordinate : string, optional  
        Y axis coordinate column name in AnnData object.

    point_size : int, optional  
        point size in the napari plot.
    
    overwrite : bool, optional  
        In the event you have multiple images in the adata object, ROI can be added to each image
        independently using the `imageid` and `subset` parameter. If you wish the results to be
        all saved with in the same column set this parameter to `False`. By default, the 
        function will overwrite the `label` column. 
        
    n_jobs : int, optional  
        Number of cores to use. Default is to use all available cores.  
    
    label : string, optional  
        Key for the returned data, stored in `adata.obs`.  

    **kwargs  
        Other arguments that can be passed to napari viewer
    
Returns:
    Modified Anndata object.

Example:
```python
    image_path = '/Users/aj/Desktop/ptcl_tma/image.ome.tif'
    adata = sm.pl.addROI_image (image_path, adata, label="Tumor Regions")
```
    """
    
    #TODO
    # - ADD Subset markers for ZARR ssection
    # - Ability to use ZARR metadata if available
    
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
        
    
    # adding option to load just the image without an adata object
    if adata is None:
        channel_names = None
    else: 
        # All operations on the AnnData object is performed first
        # Plot only the Image that is requested
        if subset is not None:
            adata = adata[adata.obs[imageid] == subset]
    
        # Recover the channel names from adata
        if channel_names == 'default':
            channel_names = adata.uns['all_markers']
        else:
            channel_names = channel_names
    
        # Index of the marker of interest and corresponding names
        if markers is None:
            idx = list(range(len(channel_names)))
            channel_names = channel_names
        else:
            idx = []
            for i in markers:
                idx.append(list(channel_names).index(i))
            channel_names = markers
        
        # Load the segmentation mask
        if seg_mask is not None:
            seg_m = tiff.imread(seg_mask)
    
 
    # Operations on the OME TIFF image is performed next
    # check the format of image
    if os.path.isfile(image_path) is True:  
        image = tiff.TiffFile(image_path, is_ome=False) #is_ome=False
        z = zarr.open(image.aszarr(), mode='r') # convert image to Zarr array
        # Identify the number of pyramids
        n_levels = len(image.series[0].levels) # pyramid
        
        # If and if not pyramids are available
        if n_levels > 1:
            pyramid = [da.from_zarr(z[i]) for i in range(n_levels)]
            multiscale = True
        else:
            pyramid = da.from_zarr(z)
            multiscale = False  
   
        # subset channels of interest
        if markers is not None:
            if n_levels > 1:
                for i in range(n_levels-1):
                    pyramid[i] = pyramid[i][idx, :, :]
                n_channels = pyramid[0].shape[0] # identify the number of channels
            else:
                pyramid = pyramid[idx, :, :]
                n_channels = pyramid.shape[0] # identify the number of channels
        else:
            if n_levels > 1:
                n_channels = pyramid[0].shape[0]
            else:
                n_channels = pyramid.shape[0]
            
        # check if channel names have been passed to all channels
        if channel_names is not None:
            assert n_channels == len(channel_names), (
                f'number of channel names ({len(channel_names)}) must '
                f'match number of channels ({n_channels})'
            )

        # Load the viewer
        viewer = napari.view_image(
            pyramid, multiscale=multiscale, channel_axis=0,
            visible=False, 
            name = None if channel_names is None else channel_names,
            **kwargs
        )
            
    # Operations on the ZARR image
    # check the format of image
    if os.path.isfile(image_path) is False: 
        #print(image_path)
        viewer = napari.Viewer()
        viewer.open(image_path, multiscale=True,
                    visible=False,
                    name = None if channel_names is None else channel_names)

    # Add the seg mask
    if seg_mask is not None:
        viewer.add_labels(seg_m, name='segmentation mask')

    # Add phenotype layer function
    def add_phenotype_layer (adata, overlay, phenotype_layer,x,y,viewer,point_size,point_color):
        coordinates = adata[adata.obs[overlay] == phenotype_layer]
        # Flip Y AXIS if needed
        if flip_y is True:
            coordinates = pd.DataFrame({'y': coordinates.obs[y],'x': coordinates.obs[x]})
        else:  
            coordinates = pd.DataFrame({'x': coordinates.obs[x],'y': coordinates.obs[y]})

        #points = coordinates.values.tolist()
        points = coordinates.values
        if point_color is None:
            r = lambda: random.randint(0,255) # random color generator
            point_color = '#%02X%02X%02X' % (r(),r(),r()) # random color generator
        viewer.add_points(points, size=point_size,face_color=point_color,visible=False,name=phenotype_layer)

    if overlay is not None:
        # categories under investigation
        if overlay_category is None:
            available_phenotypes = list(adata.obs[overlay].unique())
        else:
            available_phenotypes = overlay_category

        # Run the function on all phenotypes
        for i in available_phenotypes:
            add_phenotype_layer (adata=adata, overlay=overlay,
                                    phenotype_layer=i, x=x_coordinate, y=y_coordinate, viewer=viewer,
                                    point_size=point_size,point_color=point_color)
     
    # Intiate an ROI layer
    shape_layer = viewer.add_shapes(name=label)
    shape_layer.mode = 'add_polygon'
    _ = show_info('Draw ROIs')
    
    
    # helper functions
    def ellipse_points_to_patch(vertex_1, vertex_2,co_vertex_1, co_vertex_2):
        
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
    
    # block the viewer until ROI is added
    a = """
        Opening Napari;
        Add shape layers (on left) to draw ROI's. 
        Rename the shape layer to give a name to your ROI
        Multiple shape layers are supported
        ROI's should not overlap
        Close Napari to save ROI's.
        """
    print(a)
    viewer.show(block=True)    
    

    # Find all the shape layers
    my_shapes = [layer for layer in viewer.layers if isinstance(layer, Shapes)]
    # loop through the layers to find their names
    shape_names = []
    added_rois = []
    for i in my_shapes:
        shape_names.append(i.name)
        added_rois.append(len(viewer.layers[i.name].data))
        
    if any(y > 0 for y in added_rois):
        # Loop through all the Shape layers and extract the vertices and shape type
        all_rois = pd.DataFrame()
        for i in shape_names:
            # return shape vertices
            ver = viewer.layers[i].data
            # return shape shape
            structure = viewer.layers[i].shape_type
            # Each layer may contain multiple ROI's with different shapes (handle that)
            napari_roi_table = pd.DataFrame(dict(vertices=[np.fliplr(v) for v in ver], type=[str(s) for s in structure], ROI=i))
                
            # convert gathered ROIs into mpatches
            for i in range(len(napari_roi_table.index)):
                # ellipse
                if napari_roi_table['type'][i] == 'ellipse':
                    napari_roi_table.loc[i:, 'mpatch'] = ellipse_points_to_patch(napari_roi_table['vertices'][i][0],
                                                                                napari_roi_table['vertices'][i][1],
                                                                                napari_roi_table['vertices'][i][2],
                                                                                napari_roi_table['vertices'][i][3])
                # polygon, rectangle, line
                elif napari_roi_table['type'][i] in ['rectangle', 'polygon', 'path']:
                    napari_roi_table.loc[i:, 'mpatch'] = mpatches.Polygon(napari_roi_table['vertices'][i], closed=True)
                else:
                    raise ValueError
    
            # merge the final ROI's across all shape layers
            all_rois = pd.concat ([all_rois, napari_roi_table])
            all_rois['id'] = range(all_rois.shape[0])
            all_rois.index = all_rois['id']
            
        # Find cells within ROI and add it to the dataframe
        def add_roi_internal (roi_id):
            roi_mpatch = all_rois[all_rois['id'] == roi_id]['mpatch'][roi_id]
            inside = sub_data[roi_mpatch.contains_points(sub_data[[x_coordinate, y_coordinate]])]
            inside['ROI_internal'] = all_rois[all_rois['id'] == roi_id]['ROI'][roi_id]
            # return
            return inside
        
        print("Identifying cells within selected ROI's")
        # all ROI cells
        roi_list = all_rois['id'].unique()
        #final_roi = list()
        #for i in roi_list: 
        #    roi_mpatch = all_rois[all_rois['id'] == i]['mpatch'][i]
        #    inside = sub_data[roi_mpatch.contains_points(sub_data[[x_coordinate, y_coordinate]])]
        #    inside['ROI_internal'] = all_rois[all_rois['id'] == i]['ROI'][i]
        #    final_roi.append(inside)
        
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



