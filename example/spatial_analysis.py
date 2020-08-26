#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 09:39:41 2020
@author: Ajit Johnson Nirmal
Spatial Analysis testing
"""

# Library
import scimap as sm
import anndata as ad
import os

# WD
os.chdir ("/Users/aj/Desktop/scimap_tutorial/")


# Data
adata = ad.read('tutorial_data.h5ad')


# Functions
# 1. spatial count
adata = sm.tl.spatial_count (adata,x_coordinate='X_centroid',y_coordinate='Y_centroid',
                           phenotype='phenotype',method='radius',radius=30,
                           imageid='imageid',subset=None,label='spatial_count_radius')

adata = sm.tl.spatial_count (adata,x_coordinate='X_centroid',y_coordinate='Y_centroid',
                           phenotype='phenotype',method='knn',radius=30,
                           imageid='imageid',subset=None,label='spatial_count_knn')

# 2. spatial expression

adata = sm.tl.spatial_expression (adata, x_coordinate='X_centroid',y_coordinate='Y_centroid',
                                method='radius', radius=30, imageid='imageid', 
                                use_raw=True,subset=None,label='spatial_expression_radius')

adata = sm.tl.spatial_expression (adata, x_coordinate='X_centroid',y_coordinate='Y_centroid',
                                method='knn', radius=30, imageid='imageid', 
                                use_raw=True,subset=None,label='spatial_expression_knn')

# 3. spatial aggregates


adata = sm.tl.spatial_aggregate (adata, x_coordinate='X_centroid',y_coordinate='Y_centroid',
                        phenotype='phenotype', method='radius', radius=30, purity = 60,
                        imageid='imageid',subset=None,label='spatial_aggregate_radius')