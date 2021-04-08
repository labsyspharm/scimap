#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 21:43:29 2021
@author: Ajit Johnson Nirmal
Lasso selector to annotate regions of interest to the anndata object
"""

# Library
import matplotlib.pyplot as plt
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
import numpy as np
import pandas as pd


# Function

def add_roi_scatter (adata, method=None, marker=None, threshold=0.5, raw=False, log=False, subset=None, imageid='imageid',
                     x_coordinate='X_centroid', y_coordinate='Y_centroid',size=5,cmap='gist_heat',
                     lasso_alpha=0.8, lasso_linewidth=3, lasso_color='red',
                     roi_name='selected_roi', label='roi_scatter'):
    
    
    # create a copy of the anndata object
    adata = combined_excluded_RJP.copy()
    bdata = adata.copy()
    
    # Subset the image of interest
    if subset is not None:
        bdata = bdata[bdata.obs[imageid] == subset]
        
    # indentify the index of the marker of interest
    if marker is not None:
        marker_index = bdata.var.index.get_loc(marker)
        
    # generate the dataframe needed with 
    if method is None:
        data = pd.DataFrame({'x': bdata.obs[x_coordinate], 'y': bdata.obs[y_coordinate]})
    else:
        if marker is None:
            raise ValueError('Please choose a marker to be able to use the selected method')
    
    # If a method is passed first create the data frame with expression values
    if method is not None and raw is False:
        data = pd.DataFrame({'x': bdata.obs[x_coordinate], 'y': bdata.obs[y_coordinate], marker: bdata.X[:,marker_index]})
        if log is True:
            data[marker] = np.log1p(data[marker])
    elif method is not None and raw is True:
        data = pd.DataFrame({'x': bdata.obs[x_coordinate], 'y': bdata.obs[y_coordinate], marker: bdata.raw.X[:,marker_index]})
        if log is True:
            data[marker] = np.log1p(data[marker])
            
    # If threshold method is passed convert the marker expression to binary
    if method == 'threshold':
        hue = np.array(data[marker])
        hue = ['red' if x >= threshold else 'grey' for x in hue]
    

    # Scatter plot start
    fig, ax = plt.subplots()

    # plot
    if method == 'threshold':
        collection = ax.scatter(x=data['x'], y=data['y'], c=hue, s=size)
    elif method == 'expression':
        collection = ax.scatter(x=data['x'], y=data['y'], c=data[marker], cmap=cmap, s=size)
    else:
        collection = ax.scatter(x=data['x'], y=data['y'], s=size)
    ax.invert_yaxis()
    plt.xticks([]) ; plt.yticks([]);
    
    
    xys = collection.get_offsets()
    ind = []
    
    def onSelect(verts):
        path = Path(verts)
        ind = np.nonzero(path.contains_points(xys))[0]
        adata.obs.loc[data.iloc[ind].index, label] = roi_name
        print('Added Selected ROI to adata.obs.' + str(label))
    
    lineprops = {'color': lasso_color, 'linewidth': lasso_linewidth, 'alpha': lasso_alpha}
    lsso = LassoSelector(ax=ax, onselect=onSelect, lineprops=lineprops)
    
    
    plt.show()
        
        




