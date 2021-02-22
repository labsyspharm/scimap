#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 13:24:11 2021
@author: Ajit Johnson Nirmal
Scatter plot with real coordinates
"""

# load Library
import numpy as np
import mpl_scatter_density
import matplotlib.pyplot as plt


# Function

def hotspot (adata, method='matplotlib'):
    # Start
    
    x = adata.obs['X_centroid']
    y = adata.obs['Y_centroid']
    c= adata.X[:,adata.var.index.get_loc('PD1')]
    c = [1 if x >= 0.5 else 0 for x in c]
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
    ax.scatter_density(x, y, dpi=80, c=c, cmap='magma')
    ax.invert_yaxis()
    plt.xticks([]) ; plt.yticks([]);