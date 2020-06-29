#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 23:47:52 2020
@author: Ajit Johnson Nirmal
View the gating in Napari viewer
"""

import napari
import pandas as pd
import tifffile as tiff
import numpy as np

def gate_finder (image_path, adata, marker_of_interest, from_gate = 6, to_gate = 8, increment = 0.1,
                 markers=None, channel_names = 'default',
                 x='X_centroid',y='Y_centroid',point_size=10,image_id=None,seg_mask=None):
    """


    Parameters
    ----------
    image_path : string
        Location to the image file..
    adata : Ann Data Object
    marker_of_interest : string
        Marker for which gate is to be defined e.g. 'CD45'.
    from_gate : int, optional (The default is 6)
        Start value gate of interest.
    to_gate : int, optional (The default is 8)
        End value of the gate of interest.
    increment : float, optional (The default is 0.1)
        Increments between the start and end values.
    markers : string, optional (The default is None)
        Additional markers to be included in the plot for evaluation.
    channel_names : list, optional (The default is `adata.uns['all_markers']`)
        List of channels in the image in the exact order as image.
    x : string, optional (The default is 'X_centroid')
        X axis coordinate column name in AnnData object.
    y : string, optional (The default is 'Y_centroid')
        Y axis coordinate column name in AnnData object.
    point_size : int, optional (The default is 10)
        point size in the napari plot.
    image_id : string, optional (The default is None)
        The ID under 'ImageId' for the image of interest.
    seg_mask : string, optional (The default is None)
        Location to the segmentation mask file.

    Example
    -------
    image_path = '/Users/aj/Desktop/ptcl_tma/image.tif'
    sm.pl.gate_finder (image_path, adata, marker_of_interest='CD45',
                 from_gate = 6, to_gate = 8, increment = 0.1,
                 markers=['DNA10'], channel_names = 'default',
                 x='X_position',y='Y_position',point_size=10,
                 image_id= '77', seg_mask=None)

    """



    # If no raw data is available make a copy
    if adata.raw is None:
        adata.raw = adata

    # Copy of the raw data if it exisits
    if adata.raw is not None:
        adata.X = adata.raw.X

    # Make a copy of the data with the marker of interest
    data = pd.DataFrame(np.log1p(adata.X), columns = adata.var.index, index= adata.obs.index)[[marker_of_interest]]

    # Generate a dataframe with various gates
    def gate (g, d):
        dd = d.values
        dd = np.where(dd < g, np.nan, dd)
        np.warnings.filterwarnings('ignore')
        dd = np.where(dd > g, 1, dd)
        dd = pd.DataFrame(dd, index = d.index, columns = ['gate-' + str(g)])
        return dd

    # Identify the list of increments
    inc = list(np.arange (from_gate, to_gate, increment))
    inc = [round(num,3) for num in inc]

    # Apply the function
    r_gate = lambda x: gate(g=x, d=data) # Create lamda function
    gated_data = list(map(r_gate, inc)) # Apply function
    # Concat all the results into a single dataframe
    gates = pd.concat(gated_data, axis=1)

    ##########################################################################
    # Visulaisation using Napari
    # Plot only the Image that is requested
    if image_id is not None:
        adata = adata[adata.obs['ImageId'] == image_id]

    # Recover the channel names from adata
    if channel_names is 'default':
        channel_names = adata.uns['all_markers']
    else:
        channel_names = channel_names


    # Index of the marker of interest and corresponding names
    if markers is None:
        idx = [0, list(channel_names).index(marker_of_interest)]
        channel_names = [adata.uns['all_markers'][0], marker_of_interest]
    elif markers == 'all':
        idx = list(range(len(channel_names)))
        channel_names = channel_names
    else:
        markers = [marker_of_interest] + list(markers)
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
        viewer.add_labels(seg_m,
                          name='segmentation mask',
                          visible = False)

    # subset the gates to include only the image of interest
    gates = gates.loc[adata.obs.index,]

    # Add phenotype layer function
    def add_phenotype_layer (adata, gates, phenotype_layer,x,y,viewer,point_size):
        cells = gates[gates[phenotype_layer] == 1].index
        coordinates = adata[cells]
        coordinates = pd.DataFrame({'y': coordinates.obs[y],'x': coordinates.obs[x]})
        points = coordinates.values.tolist()
        viewer.add_points(points, size=point_size, face_color='white',visible=False,name=phenotype_layer)

    # Run the function on all phenotypes
    for i in gates.columns:
        add_phenotype_layer (adata=adata, gates=gates, phenotype_layer=i, x=x, y=y, viewer=viewer,
                                                           point_size=point_size)
