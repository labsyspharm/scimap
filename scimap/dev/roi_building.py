#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 20:23:12 2025

@author: aj
"""


import scimap as sm
import anndata as ad


adata = sm.pp.mcmicro_to_scimap(feature_table_path=['/Users/aj/Downloads/exampleSpatialTable.csv','/Users/aj/Downloads/exampleSpatialTable2.csv'])

adata = sm.tl.cluster(adata, k=4)


import anndata as ad
adata = ad.read_h5ad('/Users/aj/Partners HealthCare Dropbox/Ajit Nirmal/nirmal lab/resources/exemplarData/mcmicroExemplar001/example_data.h5ad')


adata = ad.read_h5ad('/Users/aj/Partners HealthCare Dropbox/Ajit Nirmal/nirmal lab/resources/exemplarData/cspotExampleData/CSPOT/csPhenotype/exampleImage_cspotPredict.ome.h5ad')


adata = ad.read_h5ad('/Users/aj/Partners HealthCare Dropbox/Ajit Nirmal/nirmal lab/ztempData/CRC45.h5ad')



adata = find_sample_region(adata, 
                           x_coordinate='X_centroid', y_coordinate='Y_centroid',
                           phenotype='phenotype')




visualize_sample_regions(adata, 
                             x_coordinate='X_centroid', 
                             y_coordinate='Y_centroid',
                             imageid='imageid',
                             #subset='exampleSpatialTable',
                             label='sample_regions',
                             downsample=10,
                             simplification_tolerance=0.0)




x = adata.uns['sample_regions']

sm.pl.spatial_scatterPlot (adata, colorBy ='phenotype',s=0.2)




def find_sample_region(adata, 
                       x_coordinate='X_position', 
                       y_coordinate='Y_position',
                       imageid='imageid',
                       subset=None,
                       phenotype=None,
                       label='sample_regions',
                       eps=2000, 
                       min_samples=100, 
                       alpha_value=0.0005, 
                       area_threshold_mm2=0.5,
                       pixel_size_x=0.3,
                       grid_size=100,
                       cleaning_z=1.3):
    """
    Computes tissue sample region(s) from cell centroids using DBSCAN clustering and alphashape.
    Only clusters with an area above area_threshold_mm2 are retained.
    
    A grid-based cleaning step is applied on the final subset of coordinates before clustering.
    In this cleaning, the coordinates are binned into a grid (with resolution defined by grid_size)
    and points in bins with low density (below [mean - cleaning_z * std] of bin counts) are removed.
    
    If a phenotype column is provided, the function computes:
      - The shape for the whole sample (ignoring phenotype differences), stored under 'whole_sample'
      - The shape for each unique phenotype category.
      
    When subset is None, the function processes each image independently (based on the
    'imageid' column). The results are saved in a nested dictionary in adata.uns[label] with the structure:
    
        {
            image_id_1: {
                <if phenotype is None>:
                    "sample_region": [<GeoJSON dict>, <GeoJSON dict>, ...],
                <if phenotype is provided>:
                    "whole_sample": [<GeoJSON dict>, <GeoJSON dict>, ...],
                    phenotype_category_1: [<GeoJSON dict>, <GeoJSON dict>, ...],
                    phenotype_category_2: [<GeoJSON dict>, <GeoJSON dict>, ...],
                    ...
            },
            image_id_2: { ... },
            ...
        }
    
    Parameters:
        adata : AnnData
            Annotated data object containing cell information. Each cell's centroid is based on pixels.
        x_coordinate : str
            Column name in adata.obs containing the x-coordinate.
        y_coordinate : str
            Column name in adata.obs containing the y-coordinate.
        imageid : str
            Column name in adata.obs containing the image identifier.
        subset : str or None, optional
            If provided, only cells with adata.obs[imageid] equal to this value are processed.
        phenotype : str or None, optional
            Column name in adata.obs for phenotype labels. If provided, shapes are computed
            for each unique phenotype category in addition to the whole sample.
        label : str, optional
            Key under which the nested dictionary will be stored in adata.uns.
        eps : float, optional
            Maximum distance (in pixels) between two samples for them to be considered as in the same neighborhood.
        min_samples : int, optional
            Minimum number of samples (cells) required to form a cluster.
        alpha_value : float, optional
            Parameter controlling the detail of the alphashape boundary.
        area_threshold_mm2 : float, optional
            Minimum area (in mm^2) for a cluster to be considered valid.
        pixel_size_x : float, optional
            The size of one pixel in mm (assumes square pixels).
        grid_size : int, optional
            Number of bins along each axis for grid-based cleaning.
        cleaning_z : float, optional
            Z-score factor used to determine the minimum bin count threshold.
    
    Returns:
        The modified AnnData object with computed sample region(s) stored in adata.uns[label].
        For each image (and for each phenotype, if provided), the stored value is a list of GeoJSON
        dictionaries representing the individual ROIs.
    """
    import numpy as np
    from sklearn.cluster import DBSCAN
    import alphashape
    from shapely.geometry import mapping

    # Helper function: compute sample region(s) for a given AnnData subset.
    def compute_sample_region(adata_subset):
        # Downsample to speed up clustering (e.g., every 5th cell)
        idx = np.arange(0, adata_subset.n_obs, 5)
        adata_subsampled = adata_subset[idx].copy()
        
        # Extract coordinates.
        x_coords = adata_subsampled.obs[x_coordinate].values
        y_coords = adata_subsampled.obs[y_coordinate].values
        coordinates = np.column_stack((x_coords, y_coords))
        
        # --- Grid-based cleaning ---
        if coordinates.shape[0] == 0:
            return None
        x = coordinates[:, 0]
        y = coordinates[:, 1]
        # Create grid edges.
        x_edges = np.linspace(x.min(), x.max(), grid_size + 1)
        y_edges = np.linspace(y.min(), y.max(), grid_size + 1)
        # Assign each point to a grid cell.
        x_idx = np.searchsorted(x_edges, x, side='right') - 1
        y_idx = np.searchsorted(y_edges, y, side='right') - 1
        # Flatten the 2D grid to 1D index.
        flat_idx = x_idx + y_idx * grid_size
        # Count number of points per bin.
        unique_bins, bin_counts = np.unique(flat_idx, return_counts=True)
        bin_count_dict = dict(zip(unique_bins, bin_counts))
        # For each point, look up its bin count.
        point_bin_counts = np.array([bin_count_dict.get(idx, 0) for idx in flat_idx])
        # Z-score based thresholding on bin counts.
        mean_count = np.mean(bin_counts)
        std_count = np.std(bin_counts)
        threshold = mean_count - cleaning_z * std_count
        # Keep only points in bins with count >= threshold.
        mask = point_bin_counts >= threshold
        coordinates = coordinates[mask]
        if coordinates.shape[0] == 0:
            return None
        
        # --- DBSCAN Clustering ---
        db = DBSCAN(eps=eps, min_samples=min_samples).fit(coordinates)
        labels = db.labels_
        unique_labels = set(labels)
        valid_labels = [lbl for lbl in unique_labels if lbl != -1]
        if not valid_labels:
            return None
        
        clusters_polygons = []
        for lbl in valid_labels:
            cluster_coords = coordinates[labels == lbl]
            shp = alphashape.alphashape(cluster_coords, alpha_value)
            if shp.is_empty:
                continue
            
            # Ensure we have a list of polygons.
            if shp.geom_type == "Polygon":
                polys = [shp]
            elif shp.geom_type == "MultiPolygon":
                polys = list(shp.geoms)
            else:
                polys = [g for g in shp.geoms if g.geom_type == "Polygon"]
            
            for poly in polys:
                area_px = poly.area
                area_mm2 = area_px * (pixel_size_x ** 2)
                if area_mm2 >= area_threshold_mm2:
                    clusters_polygons.append(poly)
        
        if not clusters_polygons:
            return None
        
        # Return a list of individual cluster geometries.
        return clusters_polygons

    # Determine which image ids to process.
    if imageid not in adata.obs.columns:
        raise ValueError(f"Column '{imageid}' not found in adata.obs.")
    if subset is not None:
        image_ids = [subset]
        print(f"Subsetting adata to cells with {imageid} == '{subset}'.")
    else:
        image_ids = adata.obs[imageid].unique()
    
    result = {}  # To store results per image.
    
    for im in image_ids:
        # Filter the AnnData object for the current image.
        adata_im = adata[adata.obs[imageid] == im].copy()
        print(f"\n... Processing image '{im}' with {adata_im.n_obs} cells ...")
        
        if phenotype is None:
            sample_regions = compute_sample_region(adata_im)
            if sample_regions is None:
                result[im] = None
                print(f"No valid sample region found for image '{im}'.")
            else:
                result[im] = [mapping(poly) for poly in sample_regions]
                print(f"Sample regions computed for image '{im}'.")
        else:
            if phenotype not in adata_im.obs.columns:
                raise ValueError(f"Phenotype column '{phenotype}' not found in adata.obs for image '{im}'.")
            pheno_dict = {}
            # Compute whole sample region (grid cleaning is applied here as well).
            whole_sample_regions = compute_sample_region(adata_im)
            if whole_sample_regions is None:
                pheno_dict["whole_sample"] = None
                print(f"No valid whole sample region found for image '{im}'.")
            else:
                pheno_dict["whole_sample"] = [mapping(poly) for poly in whole_sample_regions]
                print(f"Whole sample regions computed for image '{im}'.")
            
            unique_phenotypes = adata_im.obs[phenotype].unique()
            for pheno in unique_phenotypes:
                adata_pheno = adata_im[adata_im.obs[phenotype] == pheno].copy()
                print(f"Processing phenotype '{pheno}' in image '{im}' with {adata_pheno.n_obs} cells.")
                sample_regions = compute_sample_region(adata_pheno)
                if sample_regions is None:
                    pheno_dict[pheno] = None
                    print(f"No valid sample region found for phenotype '{pheno}' in image '{im}'.")
                else:
                    pheno_dict[pheno] = [mapping(poly) for poly in sample_regions]
                    print(f"Sample regions computed for phenotype '{pheno}' in image '{im}'.")
            result[im] = pheno_dict

    # Store the nested dictionary in adata.uns under the provided label.
    adata.uns[label] = result
    print(f"\nAll sample regions stored in adata.uns under key '{label}'.")
    
    return adata












def find_sample_region(adata, 
                       x_coordinate='X_centroid', 
                       y_coordinate='Y_centroid',
                       imageid='imageid',
                       subset=None,
                       phenotype=None,
                       label='sample_regions',
                       eps=1000, 
                       min_samples=100, 
                       alpha_value=0.001, 
                       area_threshold_mm2=0.5,
                       pixel_size_x=0.3):
    """
    Computes tissue sample region(s) from cell centroids using DBSCAN clustering and alphashape.
    Only clusters with an area above area_threshold_mm2 are retained.
    
    If a phenotype column is provided (i.e., phenotype is not None), the function computes:
      - The shape for the whole sample (ignoring phenotype differences), stored under 'whole_sample'
      - The shape for each unique phenotype category.
      
    When subset is None, the function processes each image independently (based on the
    'imageid' column). The results are saved in a nested dictionary in adata.uns[label] with the structure:
    
        {
            image_id_1: {
                <if phenotype is None>:
                    "sample_region": [<GeoJSON dict>, <GeoJSON dict>, ...],
                <if phenotype is provided>:
                    "whole_sample": [<GeoJSON dict>, <GeoJSON dict>, ...],
                    phenotype_category_1: [<GeoJSON dict>, <GeoJSON dict>, ...],
                    phenotype_category_2: [<GeoJSON dict>, <GeoJSON dict>, ...],
                    ...
            },
            image_id_2: { ... },
            ...
        }
    
    Parameters:
        adata : AnnData
            Annotated data object containing cell information.
        x_coordinate : str
            Column name in adata.obs containing the x-coordinate.
        y_coordinate : str
            Column name in adata.obs containing the y-coordinate.
        imageid : str
            Column name in adata.obs containing the image identifier.
        subset : str or None, optional
            If provided, only cells with adata.obs[imageid] equal to this value are processed.
        phenotype : str or None, optional
            Column name in adata.obs for phenotype labels. If provided, shapes are computed
            for each unique phenotype category in addition to the whole sample.
        label : str, optional
            Key under which the nested dictionary will be stored in adata.uns.
        eps : float, optional
            Maximum distance between two samples for them to be considered as in the same neighborhood.
        min_samples : int, optional
            Minimum number of samples required to form a cluster.
        alpha_value : float, optional
            Parameter controlling the detail of the alphashape boundary.
        area_threshold_mm2 : float, optional
            Minimum area (in mm^2) for a cluster to be considered valid.
        pixel_size_x : float, optional
            The size of one pixel in mm (assumes square pixels).
    
    Returns:
        The modified AnnData object with computed sample region(s) stored in adata.uns[label].
        For each image (and for each phenotype, if provided), the stored value is a list of GeoJSON
        dictionaries representing the individual ROIs.
    """
    import numpy as np
    from sklearn.cluster import DBSCAN
    import alphashape
    from shapely.geometry import mapping

    # Helper function: compute sample region(s) for a given AnnData subset.
    def compute_sample_region(adata_subset):
        # Downsample to speed up clustering (e.g., every 5th cell)
        idx = np.arange(0, adata_subset.n_obs, 5)
        adata_subsampled = adata_subset.copy() #[idx]
        
        # Extract coordinates.
        x_coords = adata_subsampled.obs[x_coordinate].values
        y_coords = adata_subsampled.obs[y_coordinate].values
        coordinates = np.column_stack((x_coords, y_coords))
        
        # Run DBSCAN.
        db = DBSCAN(eps=eps, min_samples=min_samples).fit(coordinates)
        labels = db.labels_
        unique_labels = set(labels)
        valid_labels = [lbl for lbl in unique_labels if lbl != -1]
        if not valid_labels:
            return None
        
        clusters_polygons = []
        for lbl in valid_labels:
            cluster_coords = coordinates[labels == lbl]
            shp = alphashape.alphashape(cluster_coords, alpha_value)
            if shp.is_empty:
                continue
            
            # Ensure we have a list of polygons.
            if shp.geom_type == "Polygon":
                polys = [shp]
            elif shp.geom_type == "MultiPolygon":
                polys = list(shp.geoms)
            else:
                polys = [g for g in shp.geoms if g.geom_type == "Polygon"]
            
            for poly in polys:
                area_px = poly.area
                area_mm2 = area_px * (pixel_size_x ** 2)
                if area_mm2 >= area_threshold_mm2:
                    clusters_polygons.append(poly)
        
        if not clusters_polygons:
            return None
        
        # Always return the list of individual cluster geometries.
        return clusters_polygons

    # Determine which image ids to process.
    if imageid not in adata.obs.columns:
        raise ValueError(f"Column '{imageid}' not found in adata.obs.")
    
    if subset is not None:
        image_ids = [subset]
        print(f"Subsetting adata to cells with {imageid} == '{subset}'.")
    else:
        image_ids = adata.obs[imageid].unique()
    
    result = {}  # To store results per image.
    
    for im in image_ids:
        # Filter the AnnData object for the current image.
        adata_im = adata[adata.obs[imageid] == im].copy()
        print(f"\n... Processing image '{im}' with {adata_im.n_obs} cells ...")
        
        if phenotype is None:
            sample_regions = compute_sample_region(adata_im)
            if sample_regions is None:
                result[im] = None
                print(f"No valid sample region found for image '{im}'.")
            else:
                result[im] = [mapping(poly) for poly in sample_regions]
                print(f"Sample regions computed for image '{im}'.")
        else:
            if phenotype not in adata_im.obs.columns:
                raise ValueError(f"Phenotype column '{phenotype}' not found in adata.obs for image '{im}'.")
            pheno_dict = {}
            # Whole sample region.
            whole_sample_regions = compute_sample_region(adata_im)
            if whole_sample_regions is None:
                pheno_dict["whole_sample"] = None
                print(f"No valid whole sample region found for image '{im}'.")
            else:
                pheno_dict["whole_sample"] = [mapping(poly) for poly in whole_sample_regions]
                print(f"Whole sample regions computed for image '{im}'.")
            
            unique_phenotypes = adata_im.obs[phenotype].unique()
            for pheno in unique_phenotypes:
                adata_pheno = adata_im[adata_im.obs[phenotype] == pheno].copy()
                print(f"Processing phenotype '{pheno}' in image '{im}' with {adata_pheno.n_obs} cells.")
                sample_regions = compute_sample_region(adata_pheno)
                if sample_regions is None:
                    pheno_dict[pheno] = None
                    print(f"No valid sample region found for phenotype '{pheno}' in image '{im}'.")
                else:
                    pheno_dict[pheno] = [mapping(poly) for poly in sample_regions]
                    print(f"Sample regions computed for phenotype '{pheno}' in image '{im}'.")
            result[im] = pheno_dict

    # Store the nested dictionary in adata.uns under the provided label.
    adata.uns[label] = result
    print(f"\nAll sample regions stored in adata.uns under key '{label}'.")
    
    return adata





def visualize_sample_regions(adata, 
                             x_coordinate='X_position', 
                             y_coordinate='Y_position',
                             imageid='imageid',
                             subset=None,
                             label='sample_regions',
                             downsample=10,
                             simplification_tolerance=0.0):
    """
    Visualizes the cell coordinates with overlays of the computed sample region(s) using Plotly.
    
    The function creates an interactive Plotly scatter plot of the cell centroids (using the
    WebGL-optimized Scattergl for improved performance on large datasets) and overlays the shape(s)
    stored in adata.uns[label]. If multiple shapes exist (e.g., when phenotype-specific regions are 
    present), each is added as a separate trace with a legend entry that can be toggled on/off.
    
    The function also accepts a subset parameter so that only one image is visualized at a time.
    Additionally, if simplification_tolerance > 0, the polygon geometries are simplified before plotting,
    reducing the number of vertices.
    
    Parameters:
        adata : AnnData
            Annotated data object containing cell information.
        x_coordinate : str
            Column name in adata.obs containing the x-coordinate.
        y_coordinate : str
            Column name in adata.obs containing the y-coordinate.
        imageid : str
            Column name in adata.obs containing the image identifier.
        subset : str or None, optional
            If provided, only cells with adata.obs[imageid] equal to this value are visualized.
            If not provided and multiple image ids are present, an error is raised.
        label : str, optional
            Key under which the nested dictionary with sample region(s) is stored in adata.uns.
        downsample : int, optional
            Downsampling factor for the scatter plot points (default is 10, to reduce load).
        simplification_tolerance : float, optional
            Tolerance for simplifying the polygon geometries. A value > 0 will reduce the number 
            of vertices in the overlay shapes.
    
    Opens an interactive Plotly plot in your browser.
    """
    import plotly.graph_objects as go
    import plotly.offline as pyo
    import numpy as np
    from shapely.geometry import shape, mapping
    
    # Ensure that we are visualizing one image at a time.
    if subset is not None:
        image_ids = [subset]
        adata_subset = adata[adata.obs[imageid] == subset].copy()
        image_id = subset
    else:
        unique_ids = adata.obs[imageid].unique()
        if len(unique_ids) != 1:
            raise ValueError("Multiple image ids found. Please provide a subset to visualize one image at a time.")
        image_id = unique_ids[0]
        adata_subset = adata.copy()

    # Check that sample region data exists.
    if label not in adata.uns:
        raise ValueError(f"No sample region data found in adata.uns under key '{label}'.")
    regions_dict = adata.uns[label]
    if image_id not in regions_dict:
        raise ValueError(f"Image id '{image_id}' not found in adata.uns['{label}'].")
    region_data = regions_dict[image_id]

    # Create the base scatter plot using Scattergl for efficient rendering.
    x_points = adata_subset.obs[x_coordinate].values[::downsample]
    y_points = adata_subset.obs[y_coordinate].values[::downsample]
    
    fig = go.Figure()
    fig.add_trace(go.Scattergl(
        x=x_points,
        y=y_points,
        mode='markers',
        marker=dict(size=3, color='gray'),
        name='Cells'
    ))
    
    # Updated helper function to add a polygon trace from a GeoJSON geometry.
    def add_polygon_trace(geojson, trace_name):
        if geojson is None:
            return
        
        # If geojson is a list, iterate over its elements.
        if isinstance(geojson, list):
            for i, single_geo in enumerate(geojson):
                add_polygon_trace(single_geo, f"{trace_name} (ROI {i+1})")
            return
        
        # If simplification is requested, simplify the geometry.
        if simplification_tolerance > 0:
            try:
                poly_obj = shape(geojson)
                simplified = poly_obj.simplify(simplification_tolerance, preserve_topology=True)
                geojson = mapping(simplified)
            except Exception as e:
                print(f"Error simplifying geometry for {trace_name}: {e}")
        
        geom_type = geojson.get("type", None)
        if geom_type == "Polygon":
            # A polygon's "coordinates" is a list of rings; we use the exterior ring (first one).
            for ring in geojson["coordinates"]:
                xs = [pt[0] for pt in ring]
                ys = [pt[1] for pt in ring]
                # Close the polygon by repeating the first point.
                xs.append(xs[0])
                ys.append(ys[0])
                fig.add_trace(go.Scattergl(
                    x=xs,
                    y=ys,
                    mode='lines',
                    fill='toself',
                    name=trace_name
                ))
        elif geom_type == "MultiPolygon":
            # For MultiPolygon, iterate over each polygon.
            for i, poly in enumerate(geojson["coordinates"]):
                for ring in poly:
                    xs = [pt[0] for pt in ring]
                    ys = [pt[1] for pt in ring]
                    xs.append(xs[0])
                    ys.append(ys[0])
                    fig.add_trace(go.Scattergl(
                        x=xs,
                        y=ys,
                        mode='lines',
                        fill='toself',
                        name=f"{trace_name} (part {i+1})"
                    ))
        else:
            print(f"Unsupported geometry type: {geom_type}")

    # Add polygon overlays.
    # If region_data is a plain GeoJSON dict (phenotype was not used), add it directly.
    if isinstance(region_data, dict) and ("type" in region_data):
        add_polygon_trace(region_data, "Sample Region")
    elif isinstance(region_data, dict):
        # region_data is assumed to be a nested dict with keys like "whole_sample" and phenotype categories.
        for key, geojson in region_data.items():
            add_polygon_trace(geojson, key)
    else:
        print("No valid sample region geometry found.")

    # Update layout settings.
    fig.update_layout(
        title=f"Sample Regions for Image {image_id}",
        xaxis_title=x_coordinate,
        yaxis_title=y_coordinate,
        showlegend=True
    )
    
    # Open the interactive plot in a web browser.
    pyo.plot(fig, auto_open=True)




import numpy as np

x= adata.obs['X_centroid'].values
y= adata.obs['Y_centroid'].values

# Assume x, y are 1D NumPy arrays with millions of points
coords = np.column_stack((x, y))

# Define grid resolution (adjust based on scale)
grid_size = 100  # e.g., 100x100 bins

x_edges = np.linspace(x.min(), x.max(), grid_size + 1)
y_edges = np.linspace(y.min(), y.max(), grid_size + 1)

# Assign each point to a grid cell
x_idx = np.searchsorted(x_edges, x, side='right') - 1
y_idx = np.searchsorted(y_edges, y, side='right') - 1

# Combine indices for counting
flat_idx = x_idx + y_idx * grid_size
unique, counts = np.unique(flat_idx, return_counts=True)

# Create a lookup of counts per grid cell
count_dict = dict(zip(unique, counts))
point_counts = np.array([count_dict.get(idx, 0) for idx in flat_idx])

# Filter points in sparse bins (e.g., less than 3 neighbors)
mask = point_counts >= 30
filtered_coords = coords[mask]


import matplotlib.pyplot as plt

# Assume x and y are 1D arrays or lists of equal length
plt.figure(figsize=(10,10))
plt.scatter(x, y, s=0.01, alpha=0.5)

plt.figure(figsize=(10,10))
plt.scatter(filtered_coords[:, 0], filtered_coords[:, 1], s=0.01, alpha=0.5)






import numpy as np
import matplotlib.pyplot as plt

x= adata.obs['X_centroid'].values
y= adata.obs['Y_centroid'].values
coords = np.column_stack((x, y))

# ---- Define the grid ----
grid_size = 100  # Adjust resolution as needed
x_edges = np.linspace(x.min(), x.max(), grid_size + 1)
y_edges = np.linspace(y.min(), y.max(), grid_size + 1)

# Assign each point to a grid cell
x_idx = np.searchsorted(x_edges, x, side='right') - 1
y_idx = np.searchsorted(y_edges, y, side='right') - 1

# Flatten the 2D grid to 1D index for counting
flat_idx = x_idx + y_idx * grid_size

# Count number of points per bin
unique_bins, bin_counts = np.unique(flat_idx, return_counts=True)

# Build a fast lookup dictionary for bin densities
bin_count_dict = dict(zip(unique_bins, bin_counts))

# For each point, look up how many points are in its bin
point_bin_counts = np.array([bin_count_dict.get(idx, 0) for idx in flat_idx])

# ---- Z-score based thresholding ----
mean_count = np.mean(bin_counts)
std_count = np.std(bin_counts)

# Define how strict to be: keep bins above (mean - z * std)
z = 1.
threshold = mean_count - z * std_count
print(f"Mean bin count: {mean_count:.2f}, Std: {std_count:.2f}, Threshold: {threshold:.2f}")

# Filter points based on bin density
mask = point_bin_counts >= threshold
filtered_coords = coords[mask]
print(f"Filtered points: {len(filtered_coords)} / {len(coords)}")

# ---- Plot before and after filtering (optional) ----
fig, axes = plt.subplots(1, 2, figsize=(12, 6))
axes[0].scatter(coords[:, 0], coords[:, 1], s=0.1, alpha=0.3)
axes[0].set_title("Original")
axes[1].scatter(filtered_coords[:, 0], filtered_coords[:, 1], s=0.1, alpha=0.3, color='green')
axes[1].set_title("After Z-score Filtering")
for ax in axes:
    ax.axis("equal")
    ax.grid(True)
plt.tight_layout()
plt.show()




