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
import scanpy as sc

from sklearn.cluster import MiniBatchKMeans

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
                        phenotype='phenotype', method='radius', radius=30, purity = 95,
                        imageid='imageid',subset=None,label='spatial_aggregate_95')

adata = sm.tl.spatial_aggregate (adata, x_coordinate='X_centroid',y_coordinate='Y_centroid',
                        phenotype='phenotype', method='knn', radius=30, purity = 60,
                        imageid='imageid',subset=None,label='spatial_aggregate_knn')


# Clustering the spatial count and expression
kmeans = MiniBatchKMeans(n_clusters=6, random_state=0).fit(adata.uns['spatial_expression_radius'])

# Rename the labels
cluster_labels = list(map(str,kmeans.labels_))
cluster_labels = list(map(lambda orig_string: 'kmeans' + '-' + orig_string, cluster_labels))
adata.obs['kmeans'] = cluster_labels

# Percent plot
percent_plot (adata,x_axis='kmeans',y_axis='phenotype',method='percent',figsize=(10, 10))

# viz clusters
import plotly.express as px
import plotly.io as pio
pio.renderers.default = 'browser'

data = pd.DataFrame({'x':adata.obs['X_centroid'], 'y':adata.obs['Y_centroid'],'col': adata.obs['spatial_aggregate_80']})
fig = px.scatter(data, x="x", y="y", color="col")
fig.update_traces(marker=dict(size=8),selector=dict(mode='markers'))
fig.update_yaxes(autorange="reversed")

# Genes overexpressed in clusters
bdata = adata.copy()
bdata.X = np.log1p(bdata.raw.X)
sc.tl.rank_genes_groups(bdata, 'kmeans', method='t-test')
sc.pl.rank_genes_groups(bdata, n_genes=5, sharey=False, fontsize=12, ncols=6)





# Save data
adata.write('tutorial_data.h5ad')
