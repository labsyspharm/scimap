#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Wed Apr  1 21:57:54 2020
# @author: Ajit Johnson Nirmal, Yuan Chen
"""
!!! abstract "Short Description"
    `sm.pl.image_viewer`: This function enables users to view OME-TIFF images within 
    the Napari viewer, offering the capability to overlay points based on any categorical 
    column, such as cluster annotations or phenotypes. It provides a dynamic and interactive 
    way to visually explore spatial distributions and relationships of cells directly 
    on the source images, enriching the analysis with spatial context and insights.

## Function
"""

# %gui qt
try:
    import napari
except:
    pass
import pandas as pd
import random
import tifffile as tiff

import dask.array as da
from dask.cache import Cache
import zarr
import os

# cache = Cache(2e9)  # Leverage two gigabytes of memory
# cache.register()


def image_viewer(
    image_path,
    adata,
    overlay=None,
    flip_y=True,
    overlay_category=None,
    markers=None,
    channel_names='default',
    x_coordinate='X_centroid',
    y_coordinate='Y_centroid',
    point_size=10,
    point_color=None,
    subset=None,
    imageid='imageid',
    seg_mask=None,
    **kwargs,
):
    """
    Parameters:
            image_path (str):
                Path to the image file. Supports TIFF, OME.TIFF, and ZARR formats.

            adata (anndata.AnnData):
                The annotated data matrix.

            overlay (str, optional):
                Column in `adata.obs` containing categorical data for visualization, such as phenotypes or clusters.

            flip_y (bool, optional):
                If True, inverts the Y-axis to match image coordinates.

            overlay_category (list, optional):
                Specific categories within `overlay` to display. If None, all categories are shown.

            markers (list, optional):
                List of markers to include in the visualization. If None, all available markers are displayed.

            channel_names (list or str, optional):
                Specifies the order of channels in the image. Default uses all markers in `adata.uns['all_markers']`.

            x_coordinate, y_coordinate (str, optional):
                Columns in `adata.obs` specifying cell coordinates. Default to 'X_centroid' and 'Y_centroid'.

            point_size (int, optional):
                Size of the points in the visualization.

            point_color (str or dict, optional):
                Color(s) for the points. Can specify a single color or a dictionary mapping categories to colors.

            imageid (str, optional):
                Column in `adata.obs` identifying different images in the dataset. Useful for datasets with multiple images.

            subset (str, optional):
                Identifier for a specific image to analyze, used in conjunction with `imageid`.

            seg_mask (str, optional):
                Path to a segmentation mask file to overlay.

            **kwargs:
                Additional arguments passed to the napari viewer.

    Returns:
            Image (napari):
                Displays the visualization using napari viewer.

    Example:
        ```python

        # Basic visualization with phenotype overlay
        sm.pl.image_viewer(image_path='/path/to/image.ome.tif', adata=adata, overlay='phenotype', point_size=5)

        # Visualization with segmentation mask and custom point colors
        sm.pl.image_viewer(image_path='/path/to/image.ome.tif', adata=adata, seg_mask='/path/to/mask.tif',
                     overlay='phenotype', point_color={'T cell': 'green', 'B cell': 'blue'}, point_size=7)

        ```
    """

    # TODO
    # - ADD Subset markers for ZARR ssection
    # - Ability to use ZARR metadata if available

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
            if (len(seg_m.shape) > 2) and (seg_m.shape[0] > 1):
                seg_m = seg_m[0]

    # Operations on the OME TIFF image is performed next
    # check the format of image
    if os.path.isfile(image_path) is True:
        image = tiff.TiffFile(image_path, is_ome=False)  # is_ome=False
        z = zarr.open(image.aszarr(), mode='r')  # convert image to Zarr array
        # Identify the number of pyramids
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
                n_channels = pyramid.shape[0]

        # check if channel names have been passed to all channels
        if channel_names is not None:
            assert n_channels == len(channel_names), (
                f'number of channel names ({len(channel_names)}) must '
                f'match number of channels ({n_channels})'
            )

        # Load the viewer
        viewer = napari.view_image(
            pyramid,
            multiscale=multiscale,
            channel_axis=0,
            visible=False,
            name=None if channel_names is None else channel_names,
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

    # Add phenotype layer function
    def add_phenotype_layer(
        adata,
        overlay,
        phenotype_layer,
        x,
        y,
        viewer,
        point_size,
        point_color,
        available_phenotypes,
    ):
        coordinates = adata[adata.obs[overlay] == phenotype_layer]
        # Flip Y AXIS if needed
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
        if point_color is None:
            r = lambda: random.randint(0, 255)  # random color generator
            point_color = '#%02X%02X%02X' % (r(), r(), r())  # random color generator
        elif isinstance(point_color, dict):
            # if dict identify the color for the given phenotype
            # also if a color is not provided in the dict assign it to white
            try:
                point_color = point_color[available_phenotypes]
            except KeyError:
                point_color = 'white'
                # if the dict has list, we need to account for it and so the following two lines
                if isinstance(point_color, list):
                    point_color = point_color[0]

        # check if point_color is a dict and if so isolate the color to the specific categoty
        viewer.add_points(
            points,
            size=point_size,
            face_color=point_color,
            visible=False,
            name=phenotype_layer,
        )

    if overlay is not None:
        # categories under investigation
        if overlay_category is None:
            available_phenotypes = list(adata.obs[overlay].unique())
        else:
            available_phenotypes = overlay_category

        # Run the function on all phenotypes
        for i in available_phenotypes:
            add_phenotype_layer(
                adata=adata,
                overlay=overlay,
                phenotype_layer=i,
                x=x_coordinate,
                y=y_coordinate,
                viewer=viewer,
                point_size=point_size,
                point_color=point_color,
                available_phenotypes=i,
            )
