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
def test_markerCorrelation (adata):
    from scimap.plotting.markerCorrelation import markerCorrelation
    saveDir = os.getcwd() + '/testFigures'
    fileName = 'markerCorrelation.png'
    markerCorrelation(adata, saveDir=saveDir, fileName=fileName)
    # check the file exist
    full_path = os.path.join(saveDir, fileName)
    assert os.path.exists(full_path), f"File was not created: {full_path}"


# groupCorrelation
def test_groupCorrelation (adata):
    from scimap.plotting.groupCorrelation import groupCorrelation
    saveDir = os.getcwd() + '/testFigures'
    fileName = 'groupCorrelation.png'
    groupCorrelation(adata, groupBy='ROI', condition='phenotype', saveDir=saveDir, fileName=fileName)
    # check the file exist
    full_path = os.path.join(saveDir, fileName)
    assert os.path.exists(full_path), f"File was not created: {full_path}"


# distPlot
def test_distPlot (adata):
    from scimap.plotting.distPlot import distPlot
    saveDir = os.getcwd() + '/testFigures'
    fileName = 'distPlot.png'
    distPlot(adata, markers=['ECAD','FOXP3'], plotGrid=True, ncols=2, saveDir=saveDir, fileName=fileName)
    # check the file exist
    full_path = os.path.join(saveDir, fileName)
    assert os.path.exists(full_path), f"File was not created: {full_path}"


# densityPlot2D
def test_densityPlot2D (adata):
    from scimap.plotting.densityPlot2D import densityPlot2D
    saveDir = os.getcwd() + '/testFigures'
    fileName = 'densityPlot2D.png'
    densityPlot2D(adata, markerA='ECAD', markerB='FOXP3',saveDir=saveDir, fileName=fileName)
    # check the file exist
    full_path = os.path.join(saveDir, fileName)
    assert os.path.exists(full_path), f"File was not created: {full_path}"


# =============================================================================
# # cluster_plots
# def test_cluster_plots (adata):
#     from scimap.plotting.cluster_plots import cluster_plots
#     output_dir = os.getcwd() + '/testFigures'
#     fileName = '_matrixplot.pdf'
#     cluster_plots(adata, group_by='phenotype', output_dir=output_dir)
#     # check the file exist
#     full_path = os.path.join(output_dir, fileName)
#     assert os.path.exists(full_path), f"File was not created: {full_path}"
# =============================================================================


# umap
def test_umap (adata):
    from scimap.plotting.umap import umap
    saveDir = os.getcwd() + '/testFigures'
    fileName = 'umap.png'
    umap(adata, color='phenotype', saveDir=saveDir, fileName=fileName)
    # check the file exist
    full_path = os.path.join(saveDir, fileName)
    assert os.path.exists(full_path), f"File was not created: {full_path}"


# foldchange
def test_foldchange (adata):
    from scimap.tools.foldchange import foldchange
    adata = foldchange(adata, from_group = 'ROI1', imageid='ROI')
    # plotting
    from scimap.plotting.foldchange import foldchange
    saveDir = os.getcwd() + '/testFigures'
    fileName = 'foldchange.png'
    foldchange(adata, label='foldchange', method='heatmap', saveDir=saveDir, fileName=fileName)
    # check the file exist
    full_path = os.path.join(saveDir, fileName)
    assert os.path.exists(full_path), f"File was not created: {full_path}"





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