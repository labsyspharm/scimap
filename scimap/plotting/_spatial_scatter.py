#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Thu Feb 18 13:24:11 2021
# @author: Ajit Johnson Nirmal

"""
!!! abstract "Short Description"
    `sm.pl.hotspot`:  The function allows users to generate a scatter plot with X and Y coordinates.

## Function
"""

# load Library
import numpy as np
import mpl_scatter_density
import matplotlib.pyplot as plt


# Function

def hotspot (adata, method='matplotlib'):
    # Start
    bdata = adata.copy()
    
    x = bdata.obs['X_centroid']
    y = bdata.obs['Y_centroid']
    c= bdata.X[:,bdata.var.index.get_loc('KI67')]
    c = [1 if x >= 0.5 else 0 for x in c]
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
    ax.scatter_density(x, y, dpi=80, c=c, cmap='magma')
    ax.invert_yaxis()
    plt.xticks([]) ; plt.yticks([]);
    
    
    
    # identify cells based on marker postivity
    bdata = sm.hl.classify(bdata, pos=['KI67','MITF'], collapse_failed=True, label='pcna_mitf')
    bdata.obs['pcna_mitf'].value_counts()
    
    c= bdata.obs['pcna_mitf'].values
    c = [1 if x == 'passed_classify' else 0 for x in c]
    
        
    c= bdata.obs['phenotype_final'].values
    c = [1 if x == 'T5b' else 0 for x in c]
