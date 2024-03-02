#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Thu Feb 11 20:19:50 2021
# @author: Ajit Johnson Nirmal
# Adapted from https://stackoverflow.com/questions/20515554/colorize-voronoi-diagram/20678647#20678647

"""
!!! abstract "Short Description"
    `sm.pl.voronoi`: This function enables the visualization of spatial data through the creation of Voronoi diagrams, 
    offering a distinctive method to explore the spatial distribution of cells or features within a defined area. 
    Users can color these diagrams according to values from any categorical column in their dataset, 
    such as cell type or tissue compartment, particularly good looking for manuscripts.

    Key considerations for optimal use of this function include:
    
    - **Application Scope**: Voronoi diagrams are particularly effective for analyzing small to moderately sized spatial regions, 
    supporting regions with up to 5,000 cells. This constraint ensures both interpretability and performance are maintained, 
    as larger regions can result in complex visualizations that are difficult to interpret and may require extended processing times.
    
    - **Performance and Interpretability**: While Voronoi diagrams offer insightful visualizations for spatial data, their 
    utility diminishes with increasing dataset size. For regions containing more than 5,000 cells, the generated plots may become 
    cluttered and challenging to interpret, alongside experiencing significant delays in generation time.
    Users are encouraged to segment larger datasets into smaller, manageable regions or utilize alternative visualization methods 
    suitable for high-density spatial data.
        

## Function
"""

# import lib
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
try:
  from shapely.geometry import MultiPoint, Point, Polygon
except:
  print("If you are using Conda in Windows please install Shapely uisng the following command to use the Voronoi function\n\nconda install -c conda-forge shapely")
from scipy.spatial import Voronoi
import matplotlib.patches as mpatches

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

# Function
def voronoi_finite_polygons_2d(vor, 
                               radius=None):

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)

# Actual function

def voronoi (adata, 
             color_by=None, 
             colors=None, 
             x_coordinate='X_centroid', 
             y_coordinate='Y_centroid',
             imageid='imageid',
             subset=None, 
             x_lim=None, 
             y_lim=None, 
             flip_y=True,
             voronoi_edge_color='black', 
             voronoi_line_width=0.1, 
             voronoi_alpha=0.5, 
             size_max=np.inf,
             overlay_points=None, 
             overlay_points_categories=None, 
             overlay_drop_categories=None, 
             overlay_points_colors=None,
             overlay_point_size = 5, 
             overlay_point_alpha= 1, 
             overlay_point_shape=".", 
             plot_legend=True, 
             legend_size = 6, **kwargs):
    """
Parameters:
        adata (anndata.AnnData): 
            An AnnData object containing the spatial data to be visualized.
        
        color_by (str, optional):  
            The name of the column used to color the Voronoi diagram. Typically, this column represents categorical 
            variables such as cell types or tissue compartments.
        
        colors (str or Dict, optional):  
            Custom color mapping for the Voronoi diagram. Can be specified as a seaborn color palette name or a dictionary 
            mapping categories to colors.
        
        x_coordinate (str, optional):  
            The column name containing the x-coordinates.
        
        y_coordinate (str, optional):  
            The column name containing the y-coordinates.
        
        imageid (str, optional):  
            The column name containing identifiers for different images or spatial contexts.
        
        subset (str, optional):  
            Specifies the identifier of a single image to focus the visualization on.
        
        x_lim (list, optional):  
            The x-axis limits for the plot as a list of two elements: [xmin, xmax].
        
        y_lim (list, optional):  
            The y-axis limits for the plot as a list of two elements: [ymin, ymax].
        
        flip_y (bool, optional):  
            If set to True, the y-axis will be flipped. This may be necessary for some datasets where y-coordinates are 
            inverted.
        
        voronoi_edge_color (str, optional):  
            The color of the edges of the Voronoi cells.
        
        voronoi_line_width (float, optional):  
            The line width of the Voronoi cell edges.
        
        voronoi_alpha (float, optional):  
            The opacity of the Voronoi cells, ranging from 0 (completely transparent) to 1 (completely opaque).
        
        size_max (float, optional): 
            The maximum size for the Voronoi cells. Can be used to limit the cell size in the visualization.
        
        overlay_points (str, optional):  
            The name of the column to use for overlaying points on the Voronoi diagram.
        
        overlay_points_categories (list, optional):  
            Specific categories within the `overlay_points` column to include in the overlay.
        
        overlay_drop_categories (list, optional):  
            Specific categories within the `overlay_points` column to exclude from the overlay.
        
        overlay_points_colors (str or Dict, optional):  
            Custom color mapping for the overlay points. Can be specified as a seaborn color palette name or a dictionary 
            mapping categories to colors.
        
        overlay_point_size (float, optional):  
            The size of the overlay points.
        
        overlay_point_alpha (float, optional):  
            The opacity of the overlay points, ranging from 0 (completely transparent) to 1 (completely opaque).
        
        overlay_point_shape (str, optional):  
            The shape of the overlay points.
        
        plot_legend (bool, optional):  
            Whether to display a legend for the plot.
        
        legend_size (float, optional):  
            The font size of the legend text.

Returns:
        Plot (matplotlib):
                Returns a plot.

Example:
        ```python
        
    
        #E xample 1: Basic Voronoi plot with default settings
    
        sm.pl.voronoi(adata)
    
        
        # Example 2: Voronoi plot colored by cell type, with customized colors and overlay points for a specific phenotype
    
        sm.pl.voronoi(adata, color_by='cell_type', colors='Set2', overlay_points='phenotype', overlay_points_colors={'phenotype1': 'red', 'phenotype2': 'blue'}, plot_legend=True)
    
        
        # Example 3: Voronoi plot for a specific image subset, with adjusted alpha and line width, and without flipping the y-axis
    
        sm.pl.voronoi(adata, subset='image_01', voronoi_alpha=0.7, voronoi_line_width=0.5)
        
        ```

    """
    
    
    # create the data frame needed
    data = adata.obs
        
    # Subset the image of interest
    if subset is not None:
        data = data[data[imageid] == subset]
        

    # subset coordinates if needed
    if x_lim is not None:
        x1 = x_lim [0]
        if len(x_lim) < 2:
            x2 = max(data[x_coordinate])
        else:
            x2 = x_lim [1]
    if y_lim is not None:
        y1 = y_lim [0]
        if len(y_lim) < 2:
            y2 = min(data[y_coordinate])
        else:
            y2 = y_lim [1]
    # do the actuall subsetting
    if x_lim is not None:
        data = data[data[x_coordinate] >= x1]
        data = data[data[x_coordinate] <= x2]
    if y_lim is not None:
        data = data[data[y_coordinate] <= y1]
        data = data[data[y_coordinate] >= y2]

                    
    # create an extra column with index information
    data['index_info'] = np.arange(data.shape[0])
        
    # generate the x and y coordinates
    points = data[[x_coordinate,y_coordinate]].values
    
    # invert the Y-axis
    if flip_y is True:
        points[:,1] = max(points[:,1])-points[:,1]
    
    # Generate colors
    if color_by is None:
        colors = np.repeat('#e5e5e5', len(data))
#    elif color_by is None and colors is not None:
#        if isinstance(colors,str):
#            colors = np.repeat(colors, len(data))
    elif color_by is not None and colors is None:
        # auto color the samples
        if len(np.unique(data[color_by])) <= 9:
            c = sns.color_palette('Set1')[0:len(np.unique(data[color_by]))]
        if len(np.unique(data[color_by])) > 9 and len(np.unique(data[color_by])) <= 20:
            c = sns.color_palette('tab20')[0:len(np.unique(data[color_by]))]
        if len(np.unique(data[color_by])) > 20:
            # For large categories generate random colors 
            np.random.seed(0) ; c = np.random.rand(len(np.unique(data[color_by])),3).tolist()
        # merge colors with phenotypes/ categories of interest
        p = np.unique(data[color_by])
        c_p = dict(zip(p, c))
        # map to colors
        colors = list(map(c_p.get, list(data[color_by].values)))
    elif color_by is not None and colors is not None:
        # check if colors is a dictionary or a sns color scale
        if isinstance(colors,str): 
            if len(sns.color_palette(colors)) < len(np.unique(data[color_by])):
                raise ValueError(str(colors) + ' includes a maximun of ' + str(len(sns.color_palette(colors))) + ' colors, while your data need '  + str(len(np.unique(data[color_by]))) + ' colors')
            else:
                c = sns.color_palette(colors)[0:len(np.unique(data[color_by]))]    
                # merge colors with phenotypes/ categories of interest
                p = np.unique(data[color_by])
                c_p = dict(zip(p, c))
        if isinstance(colors,dict):
            if len(colors) < len(np.unique(data[color_by])):
                raise ValueError('Color mapping is not provided for all categories. Please check')
            else:
                c_p = colors
        # map to colors
        colors = list(map(c_p.get, list(data[color_by].values)))
        
        
    
    # create the voronoi object
    vor = Voronoi(points)
    
    # trim the object
    regions, vertices = voronoi_finite_polygons_2d(vor)
    
    # plotting
    pts = MultiPoint([Point(i) for i in points])
    mask = pts.convex_hull
    new_vertices = []
    if type(voronoi_alpha)!=list:
        voronoi_alpha = [voronoi_alpha]*len(points)
    areas = []
    for i,(region,alph) in enumerate(zip(regions,voronoi_alpha)):
        polygon = vertices[region]
        shape = list(polygon.shape)
        shape[0] += 1
        p = Polygon(np.append(polygon, polygon[0]).reshape(*shape)).intersection(mask)
        areas+=[p.area]
        if p.area <size_max:
            poly = np.array(list(zip(p.boundary.coords.xy[0][:-1], p.boundary.coords.xy[1][:-1])))
            new_vertices.append(poly)
            if voronoi_edge_color == 'facecolor':
                plt.fill(*zip(*poly), alpha=alph, edgecolor=colors[i], linewidth = voronoi_line_width , facecolor = colors[i])
                plt.xticks([]) ; plt.yticks([]);
            else:
                plt.fill(*zip(*poly), alpha=alph, edgecolor=voronoi_edge_color, linewidth = voronoi_line_width, facecolor = colors[i])
                plt.xticks([]) ; plt.yticks([]);
                #plt.xlim([1097.5,1414.5])
                #plt.ylim([167.3,464.1])

    
    # Add scatter on top of the voronoi if user requests
    if overlay_points is not None:
        if overlay_points_categories is None:
            d = data
        if overlay_points_categories is not None:
            # convert to list if needed (cells to keep)
            if isinstance(overlay_points_categories,str): 
                overlay_points_categories = [overlay_points_categories]
            # subset cells needed
            d = data[data[overlay_points].isin(overlay_points_categories)]    
        if overlay_drop_categories is not None:
            # conver to list if needed (cells to drop)
            if isinstance(overlay_drop_categories,str): 
                overlay_drop_categories = [overlay_drop_categories]
            # subset cells needed
            d = d[-d[overlay_points].isin(overlay_drop_categories)]
        
        # Find the x and y coordinates for the overlay category
        #points_scatter = d[[x_coordinate,y_coordinate]].values
        points_scatter = points[d.index_info.values]
    
        # invert the Y-axis
        #points_scatter[:,1] = max(points_scatter[:,1])-points_scatter[:,1]
        
        # Generate colors for the scatter plot
        if overlay_points_colors is None and color_by == overlay_points:
            # Borrow color from vornoi
            wanted_keys = np.unique(d[overlay_points]) # The keys to extract
            c_p_scatter = dict((k, c_p[k]) for k in wanted_keys if k in c_p)
        elif overlay_points_colors is None and color_by != overlay_points:
            # Randomly generate colors for all the categories in scatter plot
            # auto color the samples
            if len(np.unique(d[overlay_points])) <= 9:
                c_scatter = sns.color_palette('Set1')[0:len(np.unique(d[overlay_points]))]
            if len(np.unique(d[overlay_points])) > 9 and len(np.unique(d[overlay_points])) <= 20:
                c_scatter = sns.color_palette('tab20')[0:len(np.unique(d[overlay_points]))]
            if len(np.unique(d[overlay_points])) > 20:
                # For large categories generate random colors 
                np.random.seed(1) ; c_scatter = np.random.rand(len(np.unique(d[overlay_points])),3).tolist()
            # merge colors with phenotypes/ categories of interest
            p_scatter = np.unique(d[overlay_points])
            c_p_scatter = dict(zip(p_scatter, c_scatter))
        elif  overlay_points_colors is not None:
            # check if the overlay_points_colors is a pallete
            if isinstance(overlay_points_colors,str):
                try:
                    c_scatter = sns.color_palette(overlay_points_colors)[0:len(np.unique(d[overlay_points]))]
                    if len(sns.color_palette(overlay_points_colors)) < len(np.unique(d[overlay_points])):
                        raise ValueError(str(overlay_points_colors) + ' pallete includes a maximun of ' + str(len(sns.color_palette(overlay_points_colors))) + ' colors, while your data (overlay_points_colors) need '  + str(len(np.unique(d[overlay_points]))) + ' colors') 
                except:
                    c_scatter = np.repeat(overlay_points_colors,len(np.unique(d[overlay_points])))   #[overlay_points_colors]
                # create a dict
                p_scatter = np.unique(d[overlay_points])
                c_p_scatter = dict(zip(p_scatter, c_scatter))
            if isinstance(overlay_points_colors,dict):
                if len(overlay_points_colors) < len(np.unique(d[overlay_points])):
                    raise ValueError('Color mapping is not provided for all categories. Please check overlay_points_colors')
                else:
                    c_p_scatter = overlay_points_colors
        # map to colors
        colors_scatter = list(map(c_p_scatter.get, list(d[overlay_points].values)))
            
        #plt.scatter(x = points_scatter[:,0], y = points_scatter[:,1], s= overlay_point_size, alpha= overlay_point_alpha, c= colors_scatter, marker=overlay_point_shape)
        plt.scatter(x = points_scatter[:,0], y = points_scatter[:,1], s= overlay_point_size, alpha= overlay_point_alpha, c= colors_scatter, marker=overlay_point_shape,**kwargs)
        plt.xticks([]) ; plt.yticks([]);


    if plot_legend is True:
        # Add legend to voronoi
        patchList = []
        for key in c_p:
                data_key = mpatches.Patch(color=c_p[key], label=key)
                patchList.append(data_key)
    
        first_legend = plt.legend(handles=patchList, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., prop={'size': legend_size})
        plt.tight_layout()
        # Add the legend manually to the current Axes.
        ax = plt.gca().add_artist(first_legend)
        
        if overlay_points is not None:
            # Add legend to scatter
            patchList_scatter = []
            for key in c_p_scatter:
                    data_key_scatter = mpatches.Patch(color=c_p_scatter[key], label=key)
                    patchList_scatter.append(data_key_scatter)
        
            plt.legend(handles=patchList_scatter, bbox_to_anchor=(-0.05, 1), loc=1, borderaxespad=0., prop={'size': legend_size})
            #plt.tight_layout()
    
    
