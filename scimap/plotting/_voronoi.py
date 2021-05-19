#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 20:19:50 2021
@author: Ajit Johnson Nirmal
Generating a voronoi diagram
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

# Function
def voronoi_finite_polygons_2d(vor, radius=None):
    
    """
    Adapted from https://stackoverflow.com/questions/20515554/colorize-voronoi-diagram/20678647#20678647
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.

    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.

    """

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

def voronoi (adata, color_by=None, colors=None, x_coordinate='X_centroid', y_coordinate='Y_centroid',
             imageid='imageid',subset=None, x_lim=None,y_lim=None,
             voronoi_edge_color = 'black',voronoi_line_width = 0.1, voronoi_alpha = 0.5, size_max=np.inf,
             overlay_points=None, overlay_points_categories=None, overlay_drop_categories=None, overlay_points_colors=None,
             overlay_point_size = 5, overlay_point_alpha= 1, overlay_point_shape=".", plot_legend=True, legend_size = 6,**kwargs):
    """
    

    Parameters
    ----------
    adata : Anndata object

    color_by : string, optional
        Color the voronoi diagram based on categorical variable (e.g. cell types or neighbourhoods).
        Pass the name of the column which contains the categorical variable.
        The default is None.
    colors : string or Dict, optional
        Custom coloring the voronoi diagram. The parameter accepts `sns color palettes` or a python dictionary
        mapping the categorical variable with the required color. The default is None.
    x_coordinate : float, required
        Column name containing the x-coordinates values. The default is 'X_centroid'.
    y_coordinate : float, required
        Column name containing the y-coordinates values. The default is 'Y_centroid'.
    imageid : string, optional
        Column name of the column containing the image id. The default is 'imageid'.
    subset : string, optional
        imageid of a single image to be subsetted for plotting. The default is None.
    voronoi_edge_color : string, optional
        A Matplotlib color for marking the edges of the voronoi. 
        If `facecolor` is passed, the edge color will always be the same as the face color.
        The default is 'black'.
    voronoi_line_width : float, optional
        The linewidth of the marker edges. Note: The default edgecolors is 'face'. You may want to change this as well. 
        The default is 0.1.
    voronoi_alpha : float, optional
        The alpha blending value, between 0 (transparent) and 1 (opaque). The default is 0.5.
    x_lim : list, optional
        Pass the x-coordinates range [x1,x2]. The default is None.
    y_lim : list, optional
        Pass the y-coordinates range [y1,y2]. The default is None.
    overlay_points : string, optional
        It is possible to overlay a scatter plot on top of the voronoi diagram.
        Pass the name of the column which contains categorical variable to be overlayed.
        The default is None.
    overlay_points_categories : list, optional
        If the passed column in `overlay_points` contains multiple categories, however the user is only
        interested in a subset of categories, those specific names can be passed as a list. By default all 
        categories will be overlayed on the voronoi diagram. The default is None.
    overlay_drop_categories : list, optional
        Similar to `overlay_points_categories`. Here for ease of use, especially if large number of categories are present.
        The user can drop a set of categories. The default is None.
    overlay_points_colors : string or dict, optional
        Similar to `colors`. 
        User can pass in a 
        a) solid color (like `black`)  
        b) sns palettes name (like `Set1`)
        c) python dictionary mapping the categories with custom colors
    The default is None.
    overlay_point_size : float, optional
        Overlay scatter plot point size. The default is 5.
    overlay_point_alpha : float, optional
        The alpha blending value for the overlay, between 0 (transparent) and 1 (opaque). The default is 1.
    overlay_point_shape : string, optional
        The marker style. marker can be either an instance of the class or the text shorthand for a particular marker.
        The default is ".".
    plot_legend : bool, optional
        Define if the figure legend should be plotted.
        Please note the figure legend may be out of view and you may need to resize the image to see it, especially 
        the legend for the scatter plot which will be on the left side of the plot.
        The default is True.
    legend_size : float, optional
        Resize the legend if needed. The default is 6.

    Example
    -------
    
    ```
    sm.pl.voronoi(adata, color_by='phenotype', colors=None, x_coordinate='X_position', y_coordinate='Y_position',
             imageid='ImageId',subset=None,
             voronoi_edge_color = 'black',voronoi_line_width = 0.2, voronoi_alpha = 0.5, size_max=np.inf,
             overlay_points='phenotype', overlay_points_categories=None, overlay_drop_categories=None,
             overlay_point_size = 5, overlay_point_alpha= 1, overlay_point_shape=".", plot_legend=False, legend_size=6)
    
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
    
    
