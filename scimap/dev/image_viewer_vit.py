#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Wed Apr  1 21:57:54 2020
# @author: Ajit Johnson Nirmal, Yuan Chen
"""
!!! abstract "Short Description"
    `sm.pl.image_viewer`: This function enables users to view OME-TIFF images using either 
    Napari or Vitessce viewers, offering the capability to overlay points based on any categorical 
    column, such as cluster annotations or phenotypes. It provides a dynamic and interactive 
    way to visually explore spatial distributions and relationships of cells directly 
    on the source images, enriching the analysis with spatial context and insights.

## Function
"""

# Imports
try:
    import napari
except:
    pass

try:
    import vitessce
    from vitessce import (
        VitessceConfig,
        CoordinationType,
        Component as cm,
        FileType,
        DataType,
        ViewType,
        MultiImageWrapper,
        OmeTiffWrapper,
        ObsLocationsWrapper,
        ImageWrapper,
    )
except:
    pass

import pandas as pd
import random
import tifffile as tiff
import dask.array as da
from dask.cache import Cache
import zarr
import os


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
    flip_y=True,
):
    """Helper function to add phenotype layers for napari"""
    coordinates = adata[adata.obs[overlay] == phenotype_layer]
    if flip_y:
        coordinates = pd.DataFrame({'y': coordinates.obs[y], 'x': coordinates.obs[x]})
    else:
        coordinates = pd.DataFrame({'x': coordinates.obs[x], 'y': coordinates.obs[y]})

    points = coordinates.values
    if point_color is None:
        r = lambda: random.randint(0, 255)
        point_color = '#%02X%02X%02X' % (r(), r(), r())
    elif isinstance(point_color, dict):
        try:
            point_color = point_color[available_phenotypes]
        except KeyError:
            point_color = 'white'
            if isinstance(point_color, list):
                point_color = point_color[0]

    viewer.add_points(
        points,
        size=point_size,
        face_color=point_color,
        visible=False,
        name=phenotype_layer,
    )


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
    backend='napari',
    embed_notebook=True,
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

        backend (str, optional):
            Visualization backend to use. Options: 'napari' (default) or 'vitessce'

        embed_notebook (bool, optional):
            If True, embeds the visualization in the notebook (default). For napari,
            this will show a static screenshot with a message about interactivity options.

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

    # Data preprocessing
    if adata is None:
        channel_names = None
    else:
        if subset is not None:
            adata = adata[adata.obs[imageid] == subset]

        if channel_names == 'default':
            channel_names = adata.uns['all_markers']

        if markers is None:
            idx = list(range(len(channel_names)))
        else:
            idx = [list(channel_names).index(i) for i in markers]
            channel_names = markers

        if seg_mask is not None:
            seg_m = tiff.imread(seg_mask)
            if (len(seg_m.shape) > 2) and (seg_m.shape[0] > 1):
                seg_m = seg_m[0]

    # Image loading
    if os.path.isfile(image_path):
        image = tiff.TiffFile(image_path, is_ome=False)
        z = zarr.open(image.aszarr(), mode='r')
        n_levels = len(image.series[0].levels)

        if n_levels > 1:
            pyramid = [da.from_zarr(z[i]) for i in range(n_levels)]
            multiscale = True
        else:
            pyramid = da.from_zarr(z)
            multiscale = False

        if markers is not None:
            if n_levels > 1:
                for i in range(n_levels - 1):
                    pyramid[i] = pyramid[i][idx, :, :]
                n_channels = pyramid[0].shape[0]
            else:
                pyramid = pyramid[idx, :, :]
                n_channels = pyramid.shape[0]
        else:
            if n_levels > 1:
                n_channels = pyramid[0].shape[0]
            else:
                n_channels = pyramid.shape[0]

        if channel_names is not None:
            assert n_channels == len(channel_names), (
                f'number of channel names ({len(channel_names)}) must '
                f'match number of channels ({n_channels})'
            )

        # Viewer creation based on backend choice
        if backend == 'napari':
            if embed_notebook:
                print(
                    "Note: Napari doesn't support interactive embedding in notebooks. For interactive visualization, either:\n"
                    "  - Use embed_notebook=False to open in a separate interactive window\n"
                    "  - Use backend='vitessce' for interactive visualization within the notebook"
                )

                viewer = napari.view_image(
                    pyramid,
                    multiscale=multiscale,
                    channel_axis=0,
                    visible=False,
                    name=None if channel_names is None else channel_names,
                    **kwargs,
                )
                return napari.utils.nbscreenshot(viewer)

            else:
                viewer = napari.view_image(
                    pyramid,
                    multiscale=multiscale,
                    channel_axis=0,
                    visible=False,
                    name=None if channel_names is None else channel_names,
                    **kwargs,
                )

        elif backend == 'vitessce':
            # Initialize Vitessce config
            vc = VitessceConfig(schema_version="1.0.15", name="Scimap Viewer")

            # Create dataset
            dataset = vc.add_dataset(name="scimap")

            # Add image data using zarr store
            dataset.add_object(
                DataType.IMAGE,
                FileType.IMAGE_OME_ZARR,
                path=z.store.path if hasattr(z.store, 'path') else str(image_path),
                name="Image",
                metadata={
                    "channel_names": (
                        channel_names if channel_names != 'default' else None
                    ),
                    "is_pyramid": n_levels > 1,
                },
            )

            # Add views
            spatial_view = vc.add_view(ViewType.SPATIAL, dataset=dataset)
            layer_controller = vc.add_view(ViewType.LAYER_CONTROLLER, dataset=dataset)
            status = vc.add_view(ViewType.STATUS, dataset=dataset)

            # Add points if overlay is specified
            if overlay is not None and adata is not None:
                coordinates = pd.DataFrame(
                    {
                        'x': adata.obs[x_coordinate],
                        'y': (
                            adata.obs[y_coordinate]
                            if not flip_y
                            else max(adata.obs[y_coordinate]) - adata.obs[y_coordinate]
                        ),
                        'category': adata.obs[overlay],
                    }
                )

                dataset.add_object(
                    DataType.CELLS,
                    FileType.OBS_LOCATIONS_CELLS_JSON,
                    cells=coordinates[['x', 'y']].values.tolist(),
                    cell_metadata={'category': coordinates['category'].values.tolist()},
                    name="Cells",
                )

            # Add coordination
            zoom_scope, target_x, target_y = vc.add_coordination(
                CoordinationType.SPATIAL_ZOOM,
                CoordinationType.SPATIAL_TARGET_X,
                CoordinationType.SPATIAL_TARGET_Y,
            )

            # Set initial view state
            zoom_scope.set_value(1.0)
            target_x.set_value(0)
            target_y.set_value(0)

            # Link views to coordination
            spatial_view.use_coordination(zoom_scope, target_x, target_y)
            layer_controller.use_coordination(zoom_scope, target_x, target_y)

            # Configure layout
            vc.layout(spatial_view | (layer_controller / status))

            # Set view properties
            spatial_view.set_props(
                disable3d=True,
                imageLayer={
                    "visible": True,
                    "opacity": 1.0,
                },
            )

            # Return the widget or standalone viewer
            return (
                vc.widget(proxy=kwargs.get('proxy', False), theme="light")
                if embed_notebook
                else vc.standalone()
            )

        else:
            raise ValueError(
                f"Unsupported backend: {backend}. Use 'napari' or 'vitessce'."
            )

        # Add phenotype layers for napari backend
        if backend == 'napari' and not embed_notebook:
            if overlay is not None and adata is not None:
                available_phenotypes = (
                    list(adata.obs[overlay].unique())
                    if overlay_category is None
                    else overlay_category
                )

                for phenotype in available_phenotypes:
                    add_phenotype_layer(
                        adata=adata,
                        overlay=overlay,
                        phenotype_layer=phenotype,
                        x=x_coordinate,
                        y=y_coordinate,
                        viewer=viewer,
                        point_size=point_size,
                        point_color=point_color,
                        available_phenotypes=phenotype,
                        flip_y=flip_y,
                    )

            if seg_mask is not None:
                viewer.add_labels(seg_m, name='segmentation mask', visible=False)

    else:
        # Handle ZARR format
        if backend == 'napari':
            viewer = napari.Viewer()
            viewer.open(
                image_path,
                multiscale=True,
                visible=False,
                name=None if channel_names is None else channel_names,
            )

            if embed_notebook:
                return napari.utils.nbscreenshot(viewer)

        elif backend == 'vitessce':
            # Initialize Vitessce config for ZARR
            vc = VitessceConfig(schema_version="1.0.15", name="Spatial View")
            dataset = vc.add_dataset(name="scimap")
            spatial_view = vc.add_view(cm.SPATIAL, dataset=dataset)

            # Add ZARR image data
            image_layer = dataset.add_object(
                DataType.IMAGE,
                FileType.ZARR,
                url=image_path,
                name="Image",
                metadata={"channel_names": channel_names} if channel_names else {},
            )

            return vc.widget() if embed_notebook else vc.standalone()
