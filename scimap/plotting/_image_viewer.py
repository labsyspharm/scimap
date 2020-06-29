#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 21:57:54 2020
@author: Ajit Johnson Nirmal
Using Napari to Visualize images overlayed with phenotypes or any categorical column
"""

#%gui qt
import napari
import pandas as pd
import random
import tifffile as tiff

def image_viewer (image_path, adata, overlay=None,
                    overlay_category=None,markers=None,channel_names='default',
                    x='X_centroid',y='Y_centroid',point_size=10,
                    point_color=None,image_id=None,seg_mask=None):
    """
    Parameters
    ----------
    image_path : string
        Location to the image file.
    seg_mask: string (The default is None)
        Location to the segmentation mask file.
    adata : AnnData Object
    overlay : string, optional (The default is None)
        Name of the column with any categorical data such as phenotypes or clusters.
    overlay_category : list, optional (The default is None)
        If only specfic categories within the overlay column is needed, pass their names as a list.
        If None, all categories will be used.
    markers : list, optional (The default is None)
        Markers to be included. If none, all markers will be displayed.
    channel_names : list, optional (The default is `adata.uns['all_markers']`)
        List of channels in the image in the exact order as image.
    x : string, optional (The default is 'X_centroid')
        X axis coordinate column name in AnnData object.
    y : string, optional (The default is 'Y_centroid')
        Y axis coordinate column name in AnnData object.
    point_size : int, optional (The default is 10)
        point size in the napari plot.
    image_id: string, optional (The default is None)
        The ID under 'imageid' to load the image of interest. Only useful when multiple images are being analyzed together.

    Returns
    -------
    None.

    Example
    -------
    image_path = '/Users/aj/Desktop/ptcl_tma/image.tif'
    sm.pl.image_viewer (image_path, adata, overlay='phenotype',overlay_category=None,
                markers=['CD31', "CD3D","DNA11",'CD19','CD45','CD163','FOXP3'],
                point_size=7,point_color='white')

    """
    # Plot only the Image that is requested
    if image_id is not None:
        adata = adata[adata.obs['imageid'] == image_id]

    # Recover the channel names from adata
    if channel_names is 'default':
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


    # Load the image
    image = tiff.imread(image_path, key = idx)

    # Load the segmentation mask
    if seg_mask is not None:
        seg_m = tiff.imread(seg_mask)

    # Load the viewer
    viewer = napari.view_image(
    image,
    #is_pyramid=False,
    channel_axis=0,
    name = None if channel_names is None else channel_names,
    visible = False)

    # Add the seg mask
    if seg_mask is not None:
        viewer.add_labels(seg_m, name='segmentation mask')

    # Add phenotype layer function
    def add_phenotype_layer (adata, overlay, phenotype_layer,x,y,viewer,point_size,point_color):
        coordinates = adata[adata.obs[overlay] == phenotype_layer]
        coordinates = pd.DataFrame({'y': coordinates.obs[y],'x': coordinates.obs[x]})
        points = coordinates.values.tolist()
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
                                    phenotype_layer=i, x=x, y=y, viewer=viewer,
                                    point_size=point_size,point_color=point_color)
