#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Tue May 12 23:47:52 2020
# @author: Ajit Johnson Nirmal
"""
!!! abstract "Short Description"
    `sm.pl.gate_finder`: This function leverages Napari to display OME-TIFF images, 
    overlaying points that assist in manually determining gating thresholds for specific markers. 
    By visualizing marker expression spatially, users can more accurately define gates. 
    Subsequently, the identified gating parameters can be applied to the dataset using `sm.pp.rescale`, 
    enabling precise control over data segmentation and analysis based on marker expression levels.

## Function
"""

try:
    import napari
except:
    pass

import pandas as pd
import tifffile as tiff
import numpy as np

import dask.array as da
from dask.cache import Cache
import zarr
import os

cache = Cache(2e9)  # Leverage two gigabytes of memory
cache.register()


def gate_finder(
    image_path,
    adata,
    marker_of_interest,
    layer='raw',
    log=True,
    from_gate=6,
    to_gate=8,
    increment=0.1,
    markers=None,
    channel_names='default',
    flip_y=True,
    x_coordinate='X_centroid',
    y_coordinate='Y_centroid',
    point_size=10,
    imageid='imageid',
    subset=None,
    seg_mask=None,
    **kwargs,
):
    """
Parameters:
        image_path (str):
            Path to the high-resolution image file (supports formats like TIFF, OME.TIFF).

        adata (anndata.AnnData):
            The annotated data matrix.

        marker_of_interest (str):
            The target marker for which the gating threshold is to be determined.

        layer (str, optional):
            Specifies the layer in `adata` containing expression data. Defaults to 'raw' for `adata.raw.X`.

        log (bool, optional):
            Applies log transformation to expression data if set to True. 

        from_gate (int, optional):
            Starting gate threshold value for the marker of interest. 

        to_gate (int, optional):
            Ending gate threshold value for the marker of interest. 

        increment (float, optional):
            Incremental step size between `from_gate` and `to_gate`. 

        markers (list, optional):
            A list of additional markers to include in visualization for context.

        channel_names (list or str, optional):
            Names of the channels in the image, in order. Defaults to 'default', using `adata.uns['all_markers']`.

        flip_y (bool, optional):
            Inverts the Y-axis to match image coordinates if set to True. Defaults to True.

        x_coordinate, y_coordinate (str, optional):
            Columns in `adata.obs` specifying cell coordinates. Defaults are 'X_centroid' and 'Y_centroid'.

        point_size (int, optional):
            Size of points in the visualization. 

        imageid (str, optional):
            Column in `adata.obs` identifying images for datasets with multiple images. 

        subset (str, optional):
            Specific image identifier for targeted analysis, typically an image ID. 

        seg_mask (str, optional):
            Path to a segmentation mask file to overlay.

        **kwargs:
            Additional arguments passed to the visualization tool.

Returns:
        Image (napari): 
            Displays the visualization using napari viewer.

Example:
    ```python
    
    # Visualize gating thresholds for CD45 on a specific image
    sm.pl.gate_finder(
        image_path='/path/to/image.ome.tif', adata=adata, marker_of_interest='CD45',
        from_gate=4, to_gate=10, increment=0.2, flip_y=False, point_size=12,
        subset='Sample1', seg_mask='/path/to/seg_mask.tif')

    # Log-transformed gating for a marker with additional markers and custom channel names
    sm.pl.gate_finder(
        image_path='/path/to/image.ome.tif', adata=adata, marker_of_interest='CD3',
        log=True, from_gate=3, to_gate=7, increment=0.1, markers=['CD19', 'CD4'],
        channel_names=['DAPI', 'CD3', 'CD19', 'CD4'], point_size=15)

    # Explore gating for multiple markers across different segments
    sm.pl.gate_finder(
        image_path='/path/to/image.ome.tif', adata=adata, marker_of_interest='CD8',
        layer='expression', from_gate=5, to_gate=9, increment=0.05, markers=['CD8', 'PD1'],
        subset='TumorRegion', seg_mask='/path/to/tumor_seg_mask.tif')
    
    ```
    """

    # If no raw data is available make a copy
    if adata.raw is None:
        adata.raw = adata

    # subset data if neede
    if subset is not None:
        if isinstance(subset, str):
            subset = [subset]
        if layer == 'raw':
            bdata = adata.copy()
            bdata.X = adata.raw.X
            bdata = bdata[bdata.obs[imageid].isin(subset)]
        else:
            bdata = adata.copy()
            bdata = bdata[bdata.obs[imageid].isin(subset)]
    else:
        bdata = adata.copy()

    # isolate the data
    if layer is None:
        data = pd.DataFrame(bdata.X, index=bdata.obs.index, columns=bdata.var.index)[
            [marker_of_interest]
        ]
    elif layer == 'raw':
        data = pd.DataFrame(
            bdata.raw.X, index=bdata.obs.index, columns=bdata.var.index
        )[[marker_of_interest]]
    else:
        data = pd.DataFrame(
            bdata.layers[layer], index=bdata.obs.index, columns=bdata.var.index
        )[[marker_of_interest]]

    if log is True:
        data = np.log1p(data)

    # Copy of the raw data if it exisits
    # if adata.raw is not None:
    #    adata.X = adata.raw.X

    # Plot only the Image that is requested
    # if subset is not None:
    #    adata = adata[adata.obs[imageid] == subset]

    # Make a copy of the data with the marker of interest
    # data = pd.DataFrame(np.log1p(adata.X), columns = adata.var.index, index= adata.obs.index)[[marker_of_interest]]

    # Generate a dataframe with various gates
    def gate(g, d):
        dd = d.values
        dd = np.where(dd < g, np.nan, dd)
        # np.warnings.filterwarnings('ignore')
        np.seterr('ignore')
        dd = np.where(dd > g, 1, dd)
        dd = pd.DataFrame(dd, index=d.index, columns=['gate-' + str(g)])
        return dd

    # Identify the list of increments
    inc = list(np.arange(from_gate, to_gate, increment))
    inc = [round(num, 3) for num in inc]

    # Apply the function
    r_gate = lambda x: gate(g=x, d=data)  # Create lamda function
    gated_data = list(map(r_gate, inc))  # Apply function
    # Concat all the results into a single dataframe
    gates = pd.concat(gated_data, axis=1)

    # Recover the channel names from adata
    if channel_names == 'default':
        channel_names = adata.uns['all_markers']
    else:
        channel_names = channel_names

    # if markers is a string convert to list
    if isinstance(markers, str):
        markers = [markers]

    # Index of the marker of interest and corresponding names
    if markers is not None:
        markers.extend([marker_of_interest])
        idx = np.where(np.isin(channel_names, markers))[0]
        channel_names = [channel_names[i] for i in idx]
    else:
        idx = list(range(len(channel_names)))
        channel_names = channel_names

    # Load the segmentation mask
    if seg_mask is not None:
        seg_m = tiff.imread(seg_mask)
        if (len(seg_m.shape) > 2) and (seg_m.shape[0] > 1):
            seg_m = seg_m[0]

    ##########################################################################
    # Visulaisation using Napari

    # load OME TIFF
    if os.path.isfile(image_path) is True:
        # Load the image
        image = tiff.TiffFile(image_path, is_ome=False)
        z = zarr.open(image.aszarr(), mode='r')  # convert image to Zarr array
        # Identify the number of pyramids and number of channels
        n_levels = len(image.series[0].levels)  # pyramid
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
                for i in range(n_levels - 1):
                    pyramid[i] = pyramid[i][idx, :, :]
                n_channels = pyramid[0].shape[0]  # identify the number of channels
            else:
                pyramid = pyramid[idx, :, :]
                n_channels = pyramid.shape[0]  # identify the number of channels
        else:
            if n_levels > 1:
                n_channels = pyramid[0].shape[0]
            else:
                n_channels = pyramid.shape[
                    0
                ]  # check if channel names have been passed to all channels
        if channel_names is not None:
            assert n_channels == len(channel_names), (
                f'number of channel names ({len(channel_names)}) must '
                f'match number of channels ({n_channels})'
            )

        # Load the viewer
        viewer = napari.view_image(
            pyramid,
            channel_axis=0,
            multiscale=multiscale,
            name=None if channel_names is None else channel_names,
            visible=False,
            **kwargs,
        )

    # Operations on the ZARR image
    # check the format of image
    if os.path.isfile(image_path) is False:
        # print(image_path)
        viewer = napari.Viewer()
        viewer.open(
            image_path,
            multiscale=True,
            visible=False,
            name=None if channel_names is None else channel_names,
        )

    # Add the seg mask
    if seg_mask is not None:
        viewer.add_labels(seg_m, name='segmentation mask', visible=False)

    # subset the gates to include only the image of interest
    gates = gates.loc[bdata.obs.index,]

    # Add gating layer
    def add_phenotype_layer(adata, gates, phenotype_layer, x, y, viewer, point_size):
        cells = gates[gates[phenotype_layer] == 1].index
        coordinates = adata[cells]
        # Flip Y axis if needed
        if flip_y is True:
            coordinates = pd.DataFrame(
                {'y': coordinates.obs[y], 'x': coordinates.obs[x]}
            )
        else:
            coordinates = pd.DataFrame(
                {'x': coordinates.obs[x], 'y': coordinates.obs[y]}
            )
        # points = coordinates.values.tolist()
        points = coordinates.values
        # import time
        # start = time.time()
        viewer.add_points(
            points,
            size=point_size,
            face_color='white',
            visible=False,
            name=phenotype_layer,
        )
        # stop = time.time()
        # print(stop-start)

    # Run the function on all gating layer
    for i in gates.columns:
        add_phenotype_layer(
            adata=bdata,
            gates=gates,
            phenotype_layer=i,
            x=x_coordinate,
            y=y_coordinate,
            viewer=viewer,
            point_size=point_size,
        )
