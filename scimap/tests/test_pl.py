#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 21:22:21 2024
@author: aj
Test plotting functions
"""

import pytest
import os
import anndata as ad
import unittest


# load data
@pytest.fixture
def adata():
    image_path = os.getcwd() + '/scimap/tests/scimapExampleData/scimapExampleData.h5ad'
    adata = ad.read_h5ad(image_path)
    return adata


# heatmap
def test_heatmap (adata):
    from scimap.plotting.heatmap import heatmap
    saveDir = os.getcwd() + '/testFigures'
    fileName = 'heatmap.png'
    heatmap(adata, groupBy='phenotype', standardScale='column', saveDir=saveDir, fileName=fileName)
    # check the file exist
    full_path = os.path.join(saveDir, fileName)
    assert os.path.exists(full_path), f"File was not created: {full_path}"
    




# markerCorrelation
# groupCorrelation
# distPlot
# densityPlot2D
# cluster_plots
# umap
# foldchange
# spatial_scatterPlot
# spatial_distance
# spatial_interaction
# spatialInteractionNetwork
# spatial_pscore
# stacked_barplot
# pie
# voronoi


# image_viewer
# addROI_image
# gate_finder