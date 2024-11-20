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

import warnings

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
from tqdm.auto import tqdm

# cache = Cache(2e9)  # Leverage two gigabytes of memory
# cache.register()


def get_marker_data(marker, adata, layer, log, verbose=False):
    # """Helper function to get consistent data for a marker"""
    if layer == 'raw':
        data = pd.DataFrame(
            adata.raw.X, index=adata.obs.index, columns=adata.var.index
        )[[marker]]
        if verbose:
            print(
                f"Raw data range - min: {float(data.min().iloc[0])}, max: {float(data.max().iloc[0])}"
            )
    else:
        data = pd.DataFrame(
            adata.layers[layer], index=adata.obs.index, columns=adata.var.index
        )[[marker]]
        if verbose:
            print(
                f"Layer data range - min: {float(data.min().iloc[0])}, max: {float(data.max().iloc[0])}"
            )

    if log:
        data = np.log1p(data)
        if verbose:
            print(
                f"After log transform - min: {float(data.min().iloc[0])}, max: {float(data.max().iloc[0])}"
            )

    return data


def initialize_gates(adata, imageid):
    # """Initialize gates DataFrame if it doesn't exist"""
    from sklearn.mixture import GaussianMixture
    from sklearn.preprocessing import StandardScaler

    # Create gates DataFrame if it doesn't exist
    if 'gates' not in adata.uns:
        print("Initializing gates with GMM...")
        adata.uns['gates'] = pd.DataFrame(
            index=adata.var.index, columns=adata.obs[imageid].unique(), dtype=float
        )
        adata.uns['gates'].iloc[:, :] = np.nan

        # Convert to list and use tqdm properly
        markers = list(adata.var.index)
        with tqdm(total=len(markers), desc="Computing gates", leave=False) as pbar:
            for marker in markers:
                # Get log-transformed data
                data = get_marker_data(marker, adata, 'raw', log=True, verbose=False)

                # Preprocess for GMM
                values = data.values.flatten()
                values = values[~np.isnan(values)]

                # Cut outliers
                p01, p99 = np.percentile(values, [0.1, 99.9])
                values = values[(values >= p01) & (values <= p99)]

                # Scale data
                scaler = StandardScaler()
                values_scaled = scaler.fit_transform(values.reshape(-1, 1))

                # Fit GMM
                gmm = GaussianMixture(n_components=3, random_state=42)
                gmm.fit(values_scaled)

                # Sort components by their means
                means = scaler.inverse_transform(gmm.means_)
                sorted_idx = np.argsort(means.flatten())
                sorted_means = means[sorted_idx]

                # Calculate gate as midpoint between middle and high components
                gate_value = np.mean([sorted_means[1], sorted_means[2]])

                # Ensure gate value is within data range
                min_val = float(data.min().iloc[0])
                max_val = float(data.max().iloc[0])
                gate_value = np.clip(gate_value, min_val, max_val)

                # Store gate value for all images
                adata.uns['gates'].loc[marker, :] = gate_value
                pbar.update(1)

    return adata


def calculate_auto_contrast(img, percentile_low=1, percentile_high=99, padding=0.1):
    # """Calculate contrast limits using histogram analysis with padding"""
    # If image is dask or zarr array, compute on smallest pyramid if available
    if isinstance(img, (da.Array, zarr.Array)):
        # Get smallest pyramid level if available
        if hasattr(img, 'shape') and len(img.shape) > 2:
            img = img[-1]  # Use smallest pyramid level
        # Compute statistics on a subset of data
        sample = img[::10, ::10]  # Sample every 10th pixel
        if hasattr(sample, 'compute'):
            sample = sample.compute()
    else:
        sample = img

    # Calculate percentiles for contrast
    low = np.percentile(sample, percentile_low)
    high = np.percentile(sample, percentile_high)

    # Add padding
    range_val = high - low
    low = max(0, low - (range_val * padding))  # Ensure we don't go below 0
    high = high + (range_val * padding)

    return low, high


def initialize_contrast_settings(
    adata, img, channel_names, imageid='imageid', subset=None
):
    #  """Initialize contrast settings if they don't exist"""
    if 'image_contrast_settings' not in adata.uns:
        print("Initializing contrast settings...")
        adata.uns['image_contrast_settings'] = {}

    current_image = adata.obs[imageid].iloc[0] if subset is None else subset

    if current_image not in adata.uns['image_contrast_settings']:
        tiff_file = img._store._source
        contrast_settings = {}

        # Use tqdm for contrast calculation progress
        for i, channel in enumerate(
            tqdm(channel_names, desc="Calculating contrast", leave=False)
        ):
            try:
                channel_data = tiff_file.series[0].pages[i].asarray()
                low, high = calculate_auto_contrast(channel_data)
                contrast_settings[channel] = {
                    'low': float(low),
                    'high': float(high),
                }
            except Exception as e:
                # Set default contrast values if calculation fails
                contrast_settings[channel] = {
                    'low': 0.0,
                    'high': 1.0,
                }

        adata.uns['image_contrast_settings'][current_image] = contrast_settings

    return adata


def napariGater(
    image_path,
    adata,
    layer='raw',
    log=True,
    x_coordinate='X_centroid',
    y_coordinate='Y_centroid',
    imageid='imageid',
    subset=None,
    flip_y=True,
    channel_names='default',
    point_size=10,
):
    """
    Parameters:
        image_path (str):
            Path to the high-resolution image file (supports formats like TIFF, OME.TIFF, ZARR).

        adata (anndata.AnnData):
            The annotated data matrix.

        layer (str, optional):
            Specifies the layer in `adata` containing expression data. Defaults to 'raw'.

        log (bool, optional):
            Applies log transformation to expression data if True. Defaults to True.

        x_coordinate, y_coordinate (str, optional):
            Columns in `adata.obs` specifying cell coordinates. Defaults are 'X_centroid' and 'Y_centroid'.

        imageid (str, optional):
            Column in `adata.obs` identifying images for datasets with multiple images. Defaults to 'imageid'.

        subset (str, optional):
            Specific image identifier for targeted analysis. Defaults to None.

        flip_y (bool, optional):
            Inverts the Y-axis to match image coordinates if True. Defaults to True.

        channel_names (list or str, optional):
            Names of the channels in the image. Defaults to 'default', using `adata.uns['all_markers']`.

        point_size (int, optional):
            Size of points in the visualization. Defaults to 10.

    Returns:
        None:
            Updates `adata.uns['gates']` with the gating thresholds.

    Example:
        ```python
        # Basic usage with default mcmicro parameters
        sm.pl.napariGater(
            image_path='path/to/image.ome.tif',
            adata=adata
        )

        # Custom settings with specific channels and coordinate columns
        sm.pl.napariGater(
            image_path='path/to/image.ome.tif',
            adata=adata,
            x_coordinate='X_position',
            y_coordinate='Y_position',
            channel_names=['DAPI', 'CD45', 'CD3', 'CD8'],
            point_size=15
        )

        # Working with specific image from a multi-image dataset
        sm.pl.napariGater(
            image_path='path/to/image.ome.tif',
            adata=adata,
            subset='sample1',
            imageid='imageid'
        )
        ```
    """
    import napari
    from magicgui import magicgui
    import time
    import warnings
    import os

    # Suppress macOS-specific warnings
    os.environ['QT_MAC_WANTS_LAYER'] = '1'

    # Show warning when function is called
    warnings.warn(
        "NOTE: napariGater() is currently in beta testing. "
        "If you encounter any issues, please report them at: "
        "https://github.com/labsyspharm/scimap/issues",
        UserWarning,
        stacklevel=2,
    )

    start_time = time.time()

    # Initialize gates with GMM if needed
    adata = initialize_gates(adata, imageid)

    print(f"Opening napari viewer...")

    # Recover the channel names from adata
    if channel_names == 'default':
        channel_names = adata.uns['all_markers']
    else:
        channel_names = channel_names

    # Load the image
    if isinstance(image_path, str):
        if image_path.endswith(('.tiff', '.tif')):
            image = tiff.TiffFile(image_path, is_ome=False)
            store = image.aszarr()
            img = zarr.open(store, mode='r')
            # Store the TiffFile object for later use
            img._store._source = image

            # Get shape from the TiffFile object
            shape = image.series[0].shape
            ndim = len(shape)
            is_multichannel = ndim > 2
            num_channels = shape[0] if is_multichannel else 1

            # print(f"Image shape: {shape}")
            # print(f"Number of channels: {num_channels}")

        elif image_path.endswith(('.zarr', '.zr')):
            img = zarr.open(image_path, mode='r')
            shape = img.shape
            ndim = len(shape)
            is_multichannel = ndim > 2
            num_channels = shape[0] if is_multichannel else 1
    else:
        img = image_path
        shape = img.shape
        ndim = len(shape)
        is_multichannel = ndim > 2
        num_channels = shape[0] if is_multichannel else 1

    # Initialize contrast settings if needed
    adata = initialize_contrast_settings(
        adata,
        img,
        channel_names,
        imageid=imageid,
        subset=subset,
    )

    # Create the viewer and add the image
    viewer = napari.Viewer()

    # Define a list of colormaps to cycle through
    colormaps = ['magenta', 'cyan', 'yellow', 'red', 'green', 'blue']

    # Add each channel as a separate layer with saved contrast settings
    current_image = adata.obs[imageid].iloc[0] if subset is None else subset

    # Get the TiffFile object
    tiff_file = img._store._source

    # Suppress progress bar output
    with tqdm(total=len(channel_names), desc="Loading channels", leave=False) as pbar:
        for i, channel_name in enumerate(channel_names):
            try:
                contrast_limits = (
                    adata.uns['image_contrast_settings'][current_image][channel_name][
                        'low'
                    ],
                    adata.uns['image_contrast_settings'][current_image][channel_name][
                        'high'
                    ],
                )

                try:
                    # Try direct page access first
                    channel_data = tiff_file.series[0].pages[i].asarray()
                except:
                    # Fallback to zarr array if direct access fails
                    channel_data = img[i]
                    if isinstance(channel_data, zarr.core.Array):
                        channel_data = channel_data[:]

                viewer.add_image(
                    channel_data,
                    name=channel_name,
                    visible=False,
                    colormap=colormaps[i % len(colormaps)],
                    blending='additive',
                    contrast_limits=contrast_limits,
                )
                pbar.update(1)
            except Exception as e:
                print(f"Failed to load channel {channel_name}: {type(e).__name__}")
                pbar.update(1)
                continue

    # Verify loaded channels
    loaded_channels = [
        layer.name for layer in viewer.layers if isinstance(layer, napari.layers.Image)
    ]
    if len(loaded_channels) != len(channel_names):
        print(
            f"\nWarning: Only loaded {len(loaded_channels)}/{len(channel_names)} channels"
        )
        missing = set(channel_names) - set(loaded_channels)
        if missing:
            print(f"Missing channels: {', '.join(missing)}")

    # Create points layer
    points_layer = viewer.add_points(
        np.zeros((0, 2)),
        size=point_size,
        face_color='white',
        name='gated_points',
        visible=True,
    )

    # Create initial marker data before creating GUI
    initial_marker = list(adata.var.index)[0]
    initial_data = get_marker_data(initial_marker, adata, 'raw', log, verbose=False)

    # Calculate initial min/max from expression values
    marker_data = pd.DataFrame(adata.raw.X, columns=adata.var.index)[initial_marker]
    if log:
        marker_data = np.log1p(marker_data)
    min_val = float(marker_data.min())
    max_val = float(marker_data.max())

    # Get initial gate value
    current_image = adata.obs[imageid].iloc[0] if subset is None else subset
    initial_gate = adata.uns['gates'].loc[initial_marker, current_image]
    if pd.isna(initial_gate) or initial_gate < min_val or initial_gate > max_val:
        initial_gate = min_val

    @magicgui(
        auto_call=True,
        marker={'choices': list(adata.var.index), 'value': initial_marker},
        gate={
            'widget_type': 'FloatSpinBox',
            'min': min_val,
            'max': max_val,
            'value': initial_gate,
            'step': 0.1,
        },
        confirm_gate={'widget_type': 'PushButton', 'text': 'Confirm Gate'},
        finish={'widget_type': 'PushButton', 'text': 'Finish Gating'},
    )
    def gate_controls(
        marker: str,
        gate: float = initial_gate,
        confirm_gate=False,
        finish=False,
    ):
        # Get data using helper function
        data = get_marker_data(marker, adata, layer, log)

        # Apply gate
        mask = data.values >= gate
        cells = data.index[mask.flatten()]

        # Update points
        coordinates = adata[cells]
        if flip_y:
            coordinates = pd.DataFrame(
                {'y': coordinates.obs[y_coordinate], 'x': coordinates.obs[x_coordinate]}
            )
        else:
            coordinates = pd.DataFrame(
                {'x': coordinates.obs[x_coordinate], 'y': coordinates.obs[y_coordinate]}
            )
        points_layer.data = coordinates.values

    # Add a separate handler for marker changes
    @gate_controls.marker.changed.connect
    def _on_marker_change(marker: str):
        # Calculate min/max from expression values
        marker_data = pd.DataFrame(adata.raw.X, columns=adata.var.index)[marker]
        if log:
            marker_data = np.log1p(marker_data)
        min_val = float(marker_data.min())
        max_val = float(marker_data.max())

        # Get existing gate value
        current_image = adata.obs[imageid].iloc[0] if subset is None else subset
        existing_gate = adata.uns['gates'].loc[marker, current_image]
        if pd.isna(existing_gate) or existing_gate < min_val or existing_gate > max_val:
            value = min_val
        else:
            value = existing_gate

        # Update the spinbox properties
        gate_controls.gate.min = min_val
        gate_controls.gate.max = max_val
        gate_controls.gate.value = value

        # Update layer visibility
        for layer in viewer.layers:
            if isinstance(layer, napari.layers.Image):
                if layer.name == marker:
                    layer.visible = True
                else:
                    layer.visible = False

        # Force viewer update
        viewer.reset_view()

    @gate_controls.confirm_gate.clicked.connect
    def _on_confirm():
        marker = gate_controls.marker.value
        gate = gate_controls.gate.value
        current_image = adata.obs[imageid].iloc[0] if subset is None else subset
        adata.uns['gates'].loc[marker, current_image] = float(gate)
        # print(f"Gate for {marker} in image {current_image} set to {gate}")

    # Add handler for finish button
    @gate_controls.finish.clicked.connect
    def _on_finish():
        viewer.close()

    # Initialize with empty points
    points_layer.data = np.zeros((0, 2))

    # Add the GUI to the viewer
    viewer.window.add_dock_widget(gate_controls)

    # Start the viewer
    napari.run()

    print(f"Napari viewer initialized in {time.time() - start_time:.2f} seconds")

    # return adata
