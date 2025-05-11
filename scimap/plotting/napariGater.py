#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Tue May 12 23:47:52 2020
# @author: Ajit Johnson Nirmal
"""
!!! abstract "Short Description"
    `sm.pl.napariGater()`: This function leverages Napari to display OME-TIFF images, 
    overlaying points that assist in manually determining gating thresholds for specific markers. 
    By visualizing marker expression spatially, users can more accurately define gates. 
    Subsequently, the identified gating parameters can be applied to the dataset using `sm.pp.rescale`, 
    enabling precise control over data segmentation and analysis based on marker expression levels.

    Replacement for `sm.pl.gate_finder()`
"""

import warnings

try:
    import napari
    from magicgui import magicgui
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from qtpy.QtWidgets import QWidget, QVBoxLayout
except ImportError:
    pass

import pandas as pd
import tifffile as tiff
import numpy as np

import dask.array as da
from dask import delayed
#from dask.cache import Cache
import zarr
import os
from tqdm.auto import tqdm

from matplotlib.figure import Figure

import time

# cache = Cache(2e9)  # Leverage two gigabytes of memory
# cache.register()


def get_marker_data(marker, adata, layer, log, verbose):
    # Always use the specified layer for data retrieval.
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
                f"Layer '{layer}' data range - min: {float(data.min().iloc[0])}, max: {float(data.max().iloc[0])}"
            )

    # Apply log transformation if requested.
    if log:
        data = np.log1p(data)
        if verbose:
            print(
                f"After log transform - min: {float(data.min().iloc[0])}, max: {float(data.max().iloc[0])}"
            )

    return data


def initialize_gates(adata, imageid, layer, log, verbose):
    from sklearn.mixture import GaussianMixture
    from sklearn.preprocessing import StandardScaler

    if 'gates' not in adata.uns:
        print("Initializing gates with GMM...")
        adata.uns['gates'] = pd.DataFrame(
            index=adata.var.index, columns=adata.obs[imageid].unique(), dtype=float
        )
        adata.uns['gates'].iloc[:, :] = np.nan

        markers = list(adata.var.index)
        with tqdm(total=len(markers), desc="Computing gates", leave=False) as pbar:
            for marker in markers:
                data = get_marker_data(marker, adata, layer, log, verbose)
                values = data.values.flatten()
                values = values[~np.isnan(values)]
                p01, p99 = np.percentile(values, [0.1, 99.9])
                values = values[(values >= p01) & (values <= p99)]
                scaler = StandardScaler()
                values_scaled = scaler.fit_transform(values.reshape(-1, 1))
                gmm = GaussianMixture(n_components=3, random_state=42)
                gmm.fit(values_scaled)
                means = scaler.inverse_transform(gmm.means_)
                sorted_idx = np.argsort(means.flatten())
                sorted_means = means[sorted_idx]
                gate_value = np.mean([sorted_means[1], sorted_means[2]])
                min_val = float(data.min().iloc[0])
                max_val = float(data.max().iloc[0])
                gate_value = np.clip(gate_value, min_val, max_val)
                adata.uns['gates'].loc[marker, :] = gate_value
                pbar.update(1)

    if 'napariGaterProvenance' not in adata.uns:
        adata.uns['napariGaterProvenance'] = {
            'manually_adjusted': {},
            'timestamp': {},
            'original_values': {},
        }
    
    current_image = adata.obs[imageid].iloc[0]
    if current_image not in adata.uns['napariGaterProvenance']['manually_adjusted']:
        adata.uns['napariGaterProvenance']['manually_adjusted'][current_image] = {}
        adata.uns['napariGaterProvenance']['timestamp'][current_image] = {}
        adata.uns['napariGaterProvenance']['original_values'][current_image] = {}
        
        for marker in adata.var.index:
            adata.uns['napariGaterProvenance']['original_values'][current_image][marker] = \
                float(adata.uns['gates'].loc[marker, current_image])
    
    return adata


def calculate_auto_contrast(img, percentile_low=1, percentile_high=99, padding=0.4, sample_size=100000):
    if isinstance(img, (da.Array, zarr.Array)):
        if hasattr(img, 'shape') and len(img.shape) > 2:
            img = img[-1]
        total_pixels = np.prod(img.shape[-2:])
        stride = max(1, int(np.sqrt(total_pixels / sample_size)))
        sample = img[..., ::stride, ::stride]
        if hasattr(sample, 'compute'):
            sample = sample.compute()
    else:
        total_pixels = np.prod(img.shape[-2:])
        stride = max(1, int(np.sqrt(total_pixels / sample_size)))
        sample = img[..., ::stride, ::stride]

    sample_flat = sample[np.isfinite(sample)].ravel()
    low = np.percentile(sample_flat, percentile_low)
    high = np.percentile(sample_flat, percentile_high)
    range_val = high - low
    low = max(0, low - (range_val * padding))
    high = high + (range_val * padding)
    if high == low:
        high = low + 1

    return float(low), float(high)


def initialize_contrast_settings(adata, img_data, channel_names, imageid='imageid', subset=None):
    if 'image_contrast_settings' not in adata.uns:
        adata.uns['image_contrast_settings'] = {}

    current_image = adata.obs[imageid].iloc[0] if subset is None else subset

    should_initialize = False
    if current_image not in adata.uns['image_contrast_settings']:
        should_initialize = True
    else:
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
                    if isinstance(img_data, list):
                        channel_data = img_data[-1][i]
                    else:
                        channel_data = img_data[i]
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
                    contrast_settings[channel] = {'low': 0.0, 'high': 100}
                finally:
                    pbar.update(1)

        adata.uns['image_contrast_settings'][current_image] = contrast_settings
        print(f"Saved contrast settings for {current_image} with {len(channel_names)} channels")

    return adata


def load_image_efficiently(image_path):
    if isinstance(image_path, str):
        if image_path.endswith(('.tiff', '.tif')):
            image = tiff.TiffFile(image_path, is_ome=False)
            z = zarr.open(image.aszarr(), mode='r')
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
    n_channels = len(channel_names)
    extended_colormaps = [colormaps[i % len(colormaps)] for i in range(n_channels)]
    contrast_limits = []
    for channel in channel_names:
        if channel in contrast_settings:
            contrast_limits.append(
                (contrast_settings[channel]['low'], contrast_settings[channel]['high'])
            )
        else:
            print(f"Warning: Missing contrast settings for channel {channel}. Using defaults.")
            contrast_limits.append((0.0, 100))
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
    verbose=False
):
    
    """    
    Parameters:
        image_path (str):
            Path to the high-resolution multi-channel image (TIFF or OME-TIFF supported).
    
        adata (anndata.AnnData):
            Annotated data matrix containing single-cell expression and spatial metadata.
    
        layer (str, optional):
            Specifies which layer in `adata` to use for expression data (e.g., 'raw' or a named layer). Defaults to 'raw'.
    
        log (bool, optional):
            If True, applies a log1p transformation to expression values before visualization. Defaults to True.
    
        x_coordinate, y_coordinate (str, optional):
            Keys in `adata.obs` specifying X and Y spatial coordinates of cells. Defaults are 'X_centroid' and 'Y_centroid'.
    
        imageid (str, optional):
            Column name in `adata.obs` indicating the image source (used for filtering and metadata grouping). Defaults to 'imageid'.
    
        subset (str, optional):
            Specific image ID or sample name to filter and visualize. If None, uses the first available entry in `adata.obs[imageid]`.
    
        flip_y (bool, optional):
            If True, inverts the Y-axis to match image coordinate system. Defaults to True.
    
        channel_names (list or str, optional):
            List of marker/channel names corresponding to the order in the image. 
            Defaults to 'default', which uses `adata.uns['all_markers']`.
    
        point_size (int, optional):
            Size of the points representing gated cells. Defaults to 10.
    
        calculate_contrast (bool, optional):
            If True, contrast settings are estimated automatically for each channel. If False, existing settings in `adata.uns` are reused or defaulted. Defaults to True.
    
        verbose (bool, optional):
            If True, prints detailed information about data ranges and transformations.
    
    Returns:
        None:
            Launches an interactive napari viewer for manual gate threshold adjustment.
            The gating values are stored in `adata.uns['gates']`, and user edits are tracked under `adata.uns['napariGaterProvenance']`.
    
    Example:
        ```python
        # Launch napariGater with default settings
        sm.pl.napariGater('/path/to/image.ome.tif', adata=adata)
    
        # Use expression layer, flip Y-axis off, and focus on a specific sample
        sm.pl.napariGater('/path/to/image.tif', adata=adata, layer='expression',
                          flip_y=False, subset='Sample_A')
    
        # Specify custom channels and disable contrast calculation
        sm.pl.napariGater('/path/to/image.tif', adata=adata,
                          channel_names=['DAPI', 'CD45', 'CD3'], calculate_contrast=False)
        ```
"""

    os.environ['QT_MAC_WANTS_LAYER'] = '1'

    warnings.warn(
        "NOTE: napariGater() is currently in beta testing. "
        "If you encounter any issues, please report them at: "
        "https://github.com/labsyspharm/scimap/issues",
        UserWarning,
        stacklevel=2,
    )

    print("Initializing...")
    start_time = time.time()

    # Initialize gates with the chosen layer and log settings.
    adata = initialize_gates(adata, imageid, layer, log, verbose)

    if channel_names == 'default':
        channel_names = adata.uns['all_markers']
    else:
        channel_names = channel_names 

    print("Loading image data...")
    img_data, tiff_file, multiscale = load_image_efficiently(image_path)
    if img_data is None:
        raise ValueError("Failed to load image data")
    
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
        if 'image_contrast_settings' not in adata.uns:
            adata.uns['image_contrast_settings'] = {}
        current_image = adata.obs[imageid].iloc[0] if subset is None else subset
        should_initialize = False
        if current_image not in adata.uns['image_contrast_settings']:
            should_initialize = True
        else:
            existing_channels = set(adata.uns['image_contrast_settings'][current_image].keys())
            new_channels = set(channel_names)
            missing_channels = new_channels - existing_channels
            if missing_channels:
                print(f"Adding default contrast settings for new channels: {missing_channels}")
                should_initialize = True
        if should_initialize:
            contrast_settings = adata.uns['image_contrast_settings'].get(current_image, {})
            for channel in channel_names:
                if channel not in contrast_settings:
                    try:
                        if isinstance(img_data, list):
                            channel_idx = channel_names.index(channel)
                            channel_data = img_data[-1][channel_idx]
                            if hasattr(channel_data, 'compute'):
                                min_val = float(channel_data.min().compute())
                                max_val = float(channel_data.max().compute())
                            else:
                                min_val = float(channel_data.min())
                                max_val = float(channel_data.max())
                        else:
                            channel_idx = channel_names.index(channel)
                            min_val = float(img_data[channel_idx].min())
                            max_val = float(img_data[channel_idx].max())
                    except Exception as e:
                        print(f"Warning: Could not determine data range for {channel}, using defaults. Error: {str(e)}")
                        min_val, max_val = 0.0, 1.0
                    contrast_settings[channel] = {'low': min_val, 'high': max_val}
            adata.uns['image_contrast_settings'][current_image] = contrast_settings
            print(f"Initialized contrast settings for {current_image} with {len(channel_names)} channels")

    print(f"Initialization completed in {time.time() - start_time:.2f} seconds")
    print("Opening napari viewer...")
    
    viewer = napari.Viewer()
    
    default_colormaps = [
        'magenta', 'cyan', 'yellow', 'red', 'green', 'blue',
        'magenta', 'cyan', 'yellow', 'red', 'green', 'blue'
    ]

    add_channels_to_viewer(
        viewer,
        img_data,
        channel_names,
        adata.uns['image_contrast_settings'][current_image],
        colormaps=default_colormaps
    )
    
    loaded_channels = [lyr.name for lyr in viewer.layers if isinstance(lyr, napari.layers.Image)]
    if len(loaded_channels) != len(channel_names):
        print(f"\nWarning: Only loaded {len(loaded_channels)}/{len(channel_names)} channels")
        missing = set(channel_names) - set(loaded_channels)
        if missing:
            print(f"Missing channels: {', '.join(missing)}")

    points_layer = viewer.add_points(
        np.zeros((0, 2)),
        size=point_size,
        face_color='white',
        name='gated_points',
        visible=True,
    )

    # Use the chosen layer when creating the initial marker data.
    initial_marker = list(adata.var.index)[0]
    initial_data = get_marker_data(initial_marker, adata, layer, log, verbose)

    # Calculate initial min/max from the specified layer using .iloc[0] for explicit float conversion
    marker_data = get_marker_data(initial_marker, adata, layer, log, verbose=False)
    min_val = round(float(marker_data.min().iloc[0]), 2)
    max_val = round(float(marker_data.max().iloc[0]), 2)
    
    current_image = adata.obs[imageid].iloc[0] if subset is None else subset
    initial_gate = adata.uns['gates'].loc[initial_marker, current_image]
    if pd.isna(initial_gate) or initial_gate < min_val or initial_gate > max_val:
        initial_gate = min_val
    else:
        initial_gate = round(float(initial_gate), 2)

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
            'value': '⚪ Not adjusted'
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
        data = get_marker_data(marker, adata, layer, log, verbose)
        if subset is not None:
            mask = adata.obs[imageid] == subset
            data = data[mask]
        mask = data.values >= gate
        cells = data.index[mask.flatten()]
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

    # Create histogram widget
    hist_fig = Figure(figsize=(3, 3))
    hist_canvas = FigureCanvas(hist_fig)
    hist_ax = hist_fig.add_subplot(111)

    # Add to layout manually using QWidget wrapper
    hist_widget = QWidget()
    hist_layout = QVBoxLayout()
    hist_layout.addWidget(hist_canvas)
    hist_widget.setLayout(hist_layout)
    hist_widget.setMinimumHeight(250)

    def update_histogram(marker: str, gate_value: float = None):
        hist_ax.clear()
        data = get_marker_data(marker, adata, layer, log, verbose)
        if subset is not None:
            data = data[adata.obs[imageid] == subset]
        flat_data = data.values.flatten()
        hist_ax.hist(flat_data, bins=80, color='gray', edgecolor='black')
    
        if gate_value is None:
            gate_value = adata.uns['gates'].loc[marker, current_image]
        if pd.notna(gate_value):
            hist_ax.axvline(gate_value, color='red', linestyle='--', label=f'Gate: {gate_value:.2f}')
            hist_ax.legend(loc='upper right', fontsize=8)
   
        title = f'{marker}'
        if log: 
            title = f'{marker} (log scale)'
        hist_ax.set_title(title)
        hist_ax.set_xlabel('Expression')
        hist_ax.set_ylabel('Cell Count')
        hist_canvas.draw()

    # Initial histogram
    update_histogram(initial_marker, initial_gate)

    @gate_controls.marker.changed.connect
    def _on_marker_change(marker: str):
        current_state = {
            'zoom': viewer.camera.zoom,
            'center': viewer.camera.center
        }
        marker_data = get_marker_data(marker, adata, layer, log, verbose=False)
        min_val = round(float(marker_data.min().iloc[0]), 2)
        max_val = round(float(marker_data.max().iloc[0]), 2)
        current_image = adata.obs[imageid].iloc[0] if subset is None else subset
        existing_gate = adata.uns['gates'].loc[marker, current_image]
        value = existing_gate if (not pd.isna(existing_gate) and min_val <= existing_gate <= max_val) else min_val
        gate_controls.gate.min = min_val
        gate_controls.gate.max = max_val
        gate_controls.gate.value = round(float(value), 2)
        for lyr in viewer.layers:
            if isinstance(lyr, napari.layers.Image):
                if lyr.name == marker:
                    lyr.visible = True
                    viewer.layers.selection.active = lyr
                    viewer.layers.selection.clear()
                    viewer.layers.selection.add(lyr)
                else:
                    lyr.visible = False
        viewer.camera.zoom = current_state['zoom']
        viewer.camera.center = current_state['center']
        current_image = adata.obs[imageid].iloc[0] if subset is None else subset
        is_adjusted = marker in adata.uns['napariGaterProvenance']['manually_adjusted'].get(current_image, {})
        if is_adjusted:
            status_text = "✓ ADJUSTED"
            timestamp = adata.uns['napariGaterProvenance']['timestamp'][current_image][marker]
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
        update_histogram(marker, gate_controls.gate.value)

    @gate_controls.gate.changed.connect
    def _on_gate_change(gate: float):
        marker = gate_controls.marker.value
        update_histogram(marker, gate)

    @gate_controls.confirm_gate.clicked.connect
    def _on_confirm():
        marker = gate_controls.marker.value
        gate = gate_controls.gate.value
        current_image = adata.obs[imageid].iloc[0] if subset is None else subset
        adata.uns['gates'].loc[marker, current_image] = float(gate)
        from datetime import datetime
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        if current_image not in adata.uns['napariGaterProvenance']['manually_adjusted']:
            adata.uns['napariGaterProvenance']['manually_adjusted'][current_image] = {}
        if current_image not in adata.uns['napariGaterProvenance']['timestamp']:
            adata.uns['napariGaterProvenance']['timestamp'][current_image] = {}
        if current_image not in adata.uns['napariGaterProvenance']['original_values']:
            adata.uns['napariGaterProvenance']['original_values'][current_image] = {}
        adata.uns['napariGaterProvenance']['manually_adjusted'][current_image][marker] = float(gate)
        adata.uns['napariGaterProvenance']['timestamp'][current_image][marker] = timestamp
        if marker not in adata.uns['napariGaterProvenance']['original_values'][current_image]:
            original_value = float(adata.uns['gates'].loc[marker, current_image])
            adata.uns['napariGaterProvenance']['original_values'][current_image][marker] = original_value
        short_timestamp = (datetime.strptime(timestamp, "%Y-%m-%d %H:%M:%S")
                           .strftime("%y-%m-%d %H:%M"))
        gate_controls.marker_status.value = f"✓ ADJUSTED ({short_timestamp})"
        print(f"Gate confirmed for {marker} at {gate:.2f}")

    @gate_controls.finish.clicked.connect
    def _on_finish():
        viewer.close()

    points_layer.data = np.zeros((0, 2))
    viewer.window.add_dock_widget(gate_controls, name='Gate Controls')
    viewer.window.add_dock_widget(hist_widget, name='Marker Histogram')
    napari.run()

    print(f"Napari viewer initialized in {time.time() - start_time:.2f} seconds")
    # Optionally, you can return adata if desired.
