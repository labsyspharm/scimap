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

    Repacement for `sm.pl.gate_finder()`

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
from dask import delayed
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

    # Initialize provenance tracking
    if 'napariGaterProvenance' not in adata.uns:
        adata.uns['napariGaterProvenance'] = {
            'manually_adjusted': {},  # Track adjusted markers per image
            'timestamp': {},          # Track when adjustments were made
            'original_values': {},    # Track original GMM values
        }
    
    # Initialize for current image if needed
    current_image = adata.obs[imageid].iloc[0]
    if current_image not in adata.uns['napariGaterProvenance']['manually_adjusted']:
        adata.uns['napariGaterProvenance']['manually_adjusted'][current_image] = {}
        adata.uns['napariGaterProvenance']['timestamp'][current_image] = {}
        adata.uns['napariGaterProvenance']['original_values'][current_image] = {}
        
        # Store initial GMM values
        for marker in adata.var.index:
            adata.uns['napariGaterProvenance']['original_values'][current_image][marker] = \
                float(adata.uns['gates'].loc[marker, current_image])
    
    return adata


def calculate_auto_contrast(img, percentile_low=1, percentile_high=99, padding=0.4, sample_size=100000):
    # """Calculate contrast limits using efficient sampling strategies.
    
    # Args:
    #     img: Input image (numpy array, dask array, or zarr array)
    #     percentile_low: Lower percentile for contrast (default: 1)
    #     percentile_high: Upper percentile for contrast (default: 99)
    #     padding: Padding factor for contrast range (default: 0.1)
    #     sample_size: Target number of pixels to sample (default: 100000)
    # """
    # Handle different array types
    if isinstance(img, (da.Array, zarr.Array)):
        # Get smallest pyramid if available
        if hasattr(img, 'shape') and len(img.shape) > 2:
            img = img[-1]  # Use smallest pyramid
        
        # Calculate stride to achieve target sample size
        total_pixels = np.prod(img.shape[-2:])  # Only consider spatial dimensions
        stride = max(1, int(np.sqrt(total_pixels / sample_size)))
        
        # Sample the data efficiently
        sample = img[..., ::stride, ::stride]  # Use ellipsis for any leading dimensions
        
        # Compute only the sampled data
        if hasattr(sample, 'compute'):
            sample = sample.compute()
    else:
        # For numpy arrays, use similar strided sampling
        total_pixels = np.prod(img.shape[-2:])
        stride = max(1, int(np.sqrt(total_pixels / sample_size)))
        sample = img[..., ::stride, ::stride]

    # Flatten the sample while preserving only finite values
    sample_flat = sample[np.isfinite(sample)].ravel()
    
    # Calculate robust statistics
    low = np.percentile(sample_flat, percentile_low)
    high = np.percentile(sample_flat, percentile_high)
    
    # Add padding with bounds checking
    range_val = high - low
    low = max(0, low - (range_val * padding))
    high = high + (range_val * padding)
    
    if high == low:  # Handle edge case of uniform values
        high = low + 1

    return float(low), float(high)




def initialize_contrast_settings(adata, img_data, channel_names, imageid='imageid', subset=None):
    #"""Initialize contrast settings for all channels"""
    if 'image_contrast_settings' not in adata.uns:
        adata.uns['image_contrast_settings'] = {}

    current_image = adata.obs[imageid].iloc[0] if subset is None else subset
    
    # Create a new entry for this image if it doesn't exist
    # or if we're loading it with a different set of channels
    should_initialize = False
    if current_image not in adata.uns['image_contrast_settings']:
        should_initialize = True
    else:
        # Check if the channels match the existing ones
        existing_channels = set(adata.uns['image_contrast_settings'][current_image].keys())
        new_channels = set(channel_names)
        if existing_channels != new_channels:
            print(f"New channel configuration detected for {current_image}. Updating contrast settings.")
            should_initialize = True

    if should_initialize:
        contrast_settings = {}
        with tqdm(total=len(channel_names), desc=f"Calculating contrast for {current_image}", leave=False) as pbar:
            for i, channel in enumerate(channel_names):
                try:
                    # Handle pyramidal vs non-pyramidal data
                    if isinstance(img_data, list):  # Pyramidal
                        channel_data = img_data[-1][i]  # Last level is smallest
                    else:  # Non-pyramidal
                        channel_data = img_data[i]
                    
                    # Calculate contrast limits
                    low, high = calculate_auto_contrast(
                        channel_data,
                        percentile_low=1,
                        percentile_high=99,
                        sample_size=100000
                    )
                    
                    contrast_settings[channel] = {
                        'low': low,
                        'high': high,
                    }
                except Exception as e:
                    print(f"Warning: Failed to calculate contrast for {channel} in {current_image}: {str(e)}")
                    contrast_settings[channel] = {'low': 0.0, 'high': 1.0}
                finally:
                    pbar.update(1)

        # Save settings for this specific image
        adata.uns['image_contrast_settings'][current_image] = contrast_settings
        print(f"Saved contrast settings for {current_image} with {len(channel_names)} channels")

    return adata


def load_image_efficiently(image_path):
    #"""Efficiently load image using zarr conversion"""
    if isinstance(image_path, str):
        if image_path.endswith(('.tiff', '.tif')):
            image = tiff.TiffFile(image_path, is_ome=False)
            z = zarr.open(image.aszarr(), mode='r')
            
            # Check for pyramids
            n_levels = len(image.series[0].levels)
            
            if n_levels > 1:
                data = [da.from_zarr(z[i]) for i in range(n_levels)]
                multiscale = True
            else:
                data = da.from_zarr(z)
                multiscale = False
                
            return data, image, multiscale
    
    return None, None, False

def add_channels_to_viewer(viewer, img_data, channel_names, contrast_settings, colormaps):
    #"""Add channels maintaining pyramid structure if available"""
    n_channels = len(channel_names)
    extended_colormaps = [colormaps[i % len(colormaps)] for i in range(n_channels)]
    
    # Get contrast limits for each channel
    contrast_limits = []
    for channel in channel_names:
        if channel in contrast_settings:
            contrast_limits.append(
                (contrast_settings[channel]['low'], contrast_settings[channel]['high'])
            )
        else:
            print(f"Warning: Missing contrast settings for channel {channel}. Using defaults.")
            contrast_limits.append((0.0, 1.0))
    
    viewer.add_image(
        img_data,
        channel_axis=0,
        name=channel_names,
        visible=False,
        colormap=extended_colormaps,
        blending='additive',
        contrast_limits=contrast_limits,
        multiscale=isinstance(img_data, list),
        rendering='mip',
        interpolation2d='nearest'
    )


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
    calculate_contrast=True,
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

        calculate_contrast (bool, optional):
            Whether to calculate contrast limits automatically. If False, uses full data range.
            Defaults to True.

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
            channel_names=['DAPI', 'CD45', 'CD3', 'CD8'], # note this much include all channels in the image but also match the names in `adata.var.index`
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

    print("Initializing...")
    start_time = time.time()

    # Initialize gates with GMM if needed
    adata = initialize_gates(adata, imageid)

    # Handle channel names#
    if channel_names == 'default':
        channel_names = adata.uns['all_markers']
    else:
        channel_names = channel_names 


    # Load image efficiently
    print("Loading image data...")
    img_data, tiff_file, multiscale = load_image_efficiently(image_path)
    if img_data is None:
        raise ValueError("Failed to load image data")
    
    
    # Initialize contrast settings if needed
    current_image = adata.obs[imageid].iloc[0] if subset is None else subset
    
    if calculate_contrast:
        print("Calculating contrast settings...")
        adata = initialize_contrast_settings(
            adata,
            img_data,
            channel_names,
            imageid=imageid,
            subset=subset,
        )
    else:
        # Initialize contrast settings structure if it doesn't exist
        if 'image_contrast_settings' not in adata.uns:
            adata.uns['image_contrast_settings'] = {}
        
        current_image = adata.obs[imageid].iloc[0] if subset is None else subset
        
        # Check if we need to initialize or update settings for this image
        should_initialize = False
        if current_image not in adata.uns['image_contrast_settings']:
            should_initialize = True
        else:
            # Check if we have settings for all current channels
            existing_channels = set(adata.uns['image_contrast_settings'][current_image].keys())
            new_channels = set(channel_names)
            missing_channels = new_channels - existing_channels
            if missing_channels:
                print(f"Adding default contrast settings for new channels: {missing_channels}")
                should_initialize = True
        
        if should_initialize:
            # Keep existing settings if any
            contrast_settings = adata.uns['image_contrast_settings'].get(current_image, {})
            
            # Add default settings for new channels
            for channel in channel_names:
                if channel not in contrast_settings:
                    # Try to get data range from the image if possible
                    try:
                        if isinstance(img_data, list):  # Pyramidal
                            channel_idx = channel_names.index(channel)
                            channel_data = img_data[-1][channel_idx]
                            if hasattr(channel_data, 'compute'):
                                min_val = float(channel_data.min().compute())
                                max_val = float(channel_data.max().compute())
                            else:
                                min_val = float(channel_data.min())
                                max_val = float(channel_data.max())
                        else:  # Non-pyramidal
                            channel_idx = channel_names.index(channel)
                            min_val = float(img_data[channel_idx].min())
                            max_val = float(img_data[channel_idx].max())
                    except Exception as e:
                        print(f"Warning: Could not determine data range for {channel}, using defaults. Error: {str(e)}")
                        min_val, max_val = 0.0, 1.0
                    
                    contrast_settings[channel] = {'low': min_val, 'high': max_val}
            
            # Save settings for this image
            adata.uns['image_contrast_settings'][current_image] = contrast_settings
            print(f"Initialized contrast settings for {current_image} with {len(channel_names)} channels")

    print(f"Initialization completed in {time.time() - start_time:.2f} seconds")
    print("Opening napari viewer...")
    
    # Create the viewer and add all channels efficiently
    viewer = napari.Viewer()
    
    default_colormaps = [
        'magenta', 'cyan', 'yellow', 'red', 'green', 'blue',
        'magenta', 'cyan', 'yellow', 'red', 'green', 'blue'
    ]  # Basic colors that will be cycled

    add_channels_to_viewer(
        viewer,
        img_data,
        channel_names,
        adata.uns['image_contrast_settings'][current_image],
        colormaps=default_colormaps
    )
    
    # Verify loaded channels
    loaded_channels = [layer.name for layer in viewer.layers if isinstance(layer, napari.layers.Image)]
    if len(loaded_channels) != len(channel_names):
        print(f"\nWarning: Only loaded {len(loaded_channels)}/{len(channel_names)} channels")
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
        layout='vertical',
        marker={
            'choices': list(adata.var.index), 
            'value': initial_marker,
            'label': 'Select Marker:'
        },
        gate={
            'widget_type': 'FloatSpinBox',
            'min': min_val,
            'max': max_val,
            'value': initial_gate,
            'step': 0.01,
            'label': 'Gate Threshold:'
        },
        marker_status={
            'widget_type': 'Label',
            'value': '⚪ Not adjusted'  # Initial value
        },
        confirm_gate={
            'widget_type': 'PushButton', 
            'text': 'Confirm Gate'
        },
        finish={
            'widget_type': 'PushButton', 
            'text': 'Finish Gating'
        },
    )
    def gate_controls(
        marker: str,
        gate: float = initial_gate,
        marker_status: str = '⚪ Not adjusted',
        confirm_gate=False,
        finish=False,
    ):
        # Get data using helper function
        data = get_marker_data(marker, adata, layer, log)
        
        # Filter data for subset if specified
        if subset is not None:
            mask = adata.obs[imageid] == subset
            data = data[mask]

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
        # Store current view state
        current_state = {
            'zoom': viewer.camera.zoom,
            'center': viewer.camera.center
        }
        
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

        # Update layer visibility and selection
        for layer in viewer.layers:
            if isinstance(layer, napari.layers.Image):
                if layer.name == marker:
                    layer.visible = True
                    viewer.layers.selection.active = layer  # Select the layer
                    viewer.layers.selection.clear()
                    viewer.layers.selection.add(layer)
                else:
                    layer.visible = False

        # Restore view state
        viewer.camera.zoom = current_state['zoom']
        viewer.camera.center = current_state['center']

        # Update status with more visible formatting and shorter timestamp
        current_image = adata.obs[imageid].iloc[0] if subset is None else subset
        is_adjusted = marker in adata.uns['napariGaterProvenance']['manually_adjusted'].get(current_image, {})
        if is_adjusted:
            status_text = "✓ ADJUSTED"
            # Get and format timestamp
            timestamp = adata.uns['napariGaterProvenance']['timestamp'][current_image][marker]
            # Convert stored timestamp to shorter format
            from datetime import datetime
            try:
                dt = datetime.strptime(timestamp, "%Y-%m-%d %H:%M:%S")
                short_timestamp = dt.strftime("%y-%m-%d %H:%M")
                status_text += f" ({short_timestamp})"
            except:
                status_text += f" ({timestamp})"
        else:
            status_text = "⚪ NOT ADJUSTED"
        
        gate_controls.marker_status.value = status_text

    @gate_controls.confirm_gate.clicked.connect
    def _on_confirm():
        marker = gate_controls.marker.value
        gate = gate_controls.gate.value
        current_image = adata.obs[imageid].iloc[0] if subset is None else subset
        
        # Update gate value
        adata.uns['gates'].loc[marker, current_image] = float(gate)
        
        # Update provenance tracking
        from datetime import datetime
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        # Ensure all necessary dictionaries exist
        if current_image not in adata.uns['napariGaterProvenance']['manually_adjusted']:
            adata.uns['napariGaterProvenance']['manually_adjusted'][current_image] = {}
        if current_image not in adata.uns['napariGaterProvenance']['timestamp']:
            adata.uns['napariGaterProvenance']['timestamp'][current_image] = {}
        if current_image not in adata.uns['napariGaterProvenance']['original_values']:
            adata.uns['napariGaterProvenance']['original_values'][current_image] = {}
        
        # Store the adjustment
        adata.uns['napariGaterProvenance']['manually_adjusted'][current_image][marker] = float(gate)
        adata.uns['napariGaterProvenance']['timestamp'][current_image][marker] = timestamp
        
        # Store original value if not already stored
        if marker not in adata.uns['napariGaterProvenance']['original_values'][current_image]:
            original_value = float(adata.uns['gates'].loc[marker, current_image])
            adata.uns['napariGaterProvenance']['original_values'][current_image][marker] = original_value
        
        # Update status with confirmation message
        short_timestamp = datetime.strptime(timestamp, "%Y-%m-%d %H:%M:%S").strftime("%y-%m-%d %H:%M")
        gate_controls.marker_status.value = f"✓ ADJUSTED ({short_timestamp})"
        
        print(f"Gate confirmed for {marker} at {gate:.2f}")

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
