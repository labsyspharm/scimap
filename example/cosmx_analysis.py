# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 13:18:24 2021
@author: Ajit Johnson Nirmal
COSMX analysis
"""

import scimap as sm
import scanpy as sc
import numpy as np
import pandas as pd


# %% Impot data
# import cosmx data
data_path = {
"D:/cosmx/Lung5_Rep1/Lung5_Rep1-Flat_files_and_images/Lung5_Rep1_exprMat_file.csv": "D:/cosmx/Lung5_Rep1/Lung5_Rep1-Flat_files_and_images/Lung5_Rep1_metadata_file.csv",
"D:/cosmx/Lung5_Rep2/Lung5_Rep2-Flat_files_and_images/Lung5_Rep2_exprMat_file.csv": "D:/cosmx/Lung5_Rep2/Lung5_Rep2-Flat_files_and_images/Lung5_Rep2_metadata_file.csv",
"D:/cosmx/Lung5_Rep3/Lung5_Rep3-Flat_files_and_images/Lung5_Rep3_exprMat_file.csv": "D:/cosmx/Lung5_Rep3/Lung5_Rep3-Flat_files_and_images/Lung5_Rep3_metadata_file.csv",
"D:/cosmx/Lung6/Lung6-Flat_files_and_images/Lung6_exprMat_file.csv": "D:/cosmx/Lung6/Lung6-Flat_files_and_images/Lung6_metadata_file.csv",
"D:/cosmx/Lung9_Rep1/Lung9_Rep1-Flat_files_and_images/Lung9_Rep1_exprMat_file.csv": "D:/cosmx/Lung9_Rep1/Lung9_Rep1-Flat_files_and_images/Lung9_Rep1_metadata_file.csv",
"D:/cosmx/Lung9_Rep2/Lung9_Rep2-Flat_files_and_images/Lung9_Rep2_exprMat_file.csv": "D:/cosmx/Lung9_Rep2/Lung9_Rep2-Flat_files_and_images/Lung9_Rep2_metadata_file.csv",
"D:/cosmx/Lung12/Lung12-Flat_files_and_images/Lung12_exprMat_file.csv": "D:/cosmx/Lung12/Lung12-Flat_files_and_images/Lung12_metadata_file.csv",
"D:/cosmx/Lung13/Lung13-Flat_files_and_images/Lung13_exprMat_file.csv": "D:/cosmx/Lung13/Lung13-Flat_files_and_images/Lung13_metadata_file.csv",
            }

data_path = {
"D:/cosmx/Lung9_Rep1/Lung9_Rep1-Flat_files_and_images/Lung9_Rep1_exprMat_file.csv": "D:/cosmx/Lung9_Rep1/Lung9_Rep1-Flat_files_and_images/Lung9_Rep1_metadata_file.csv",
  }


plt.scatter(x=abs(umap_coordinates['umap-1']), y=umap_coordinates['umap-2'], s=1)
plt.scatter(x=umap_coordinates['umap-1'], y=umap_coordinates['umap-2'], s=1)


adata = cosmx_to_scimap (data_path,
                     unique_CellId=True,
                     remove_NegPrb=True,
                     random_sample=50000,
                     log=False)

adata.write('D:/cosmx/combined.h5ad')
adata.write('D:/cosmx/combined_subset.h5ad')
adata.write('D:/cosmx/Lung9_Rep2.h5ad')

# read
adata = sc.read('D:/cosmx/combined_subset.h5ad')
adata = sc.read('D:/cosmx/Lung9_Rep2.h5ad')

#%% Process


# animate
plotly (adata,phenotype='kmeans',image_id=None,x='CenterX_global_px',y='CenterY_global_px',size=4)
# umap
adata = sm.tl.umap (adata)
# cluster
adata = sm.tl.cluster (adata, k= 7, method = 'kmeans', label='kmeans', use_raw=False)

animate (adata, color='kmeans',subsample=None, s=2, figsize=(25,25),
         n_frames=100, pltStyle='dark_background',title=True, fontsize=35,
         x_coordinate='CenterX_global_px', y_coordinate='CenterY_global_px',
               save_animation = "C:/Users/ajn16/Downloads/cosmx")





adata.obs['fov'] = adata.obs['fov'].astype('category')

# basic filtering
sc.pp.filter_cells(adata, min_genes=30)
sc.pp.filter_genes(adata, min_cells=10)



# Norm
sc.pp.normalize_total(adata, target_sum=10000, exclude_highly_expressed=True)

# take square root of data
#adata.X = np.sqrt(adata.X)
sc.pp.log1p(adata)

# save raw
adata.raw = adata

# Find MVGs
sc.pp.highly_variable_genes(adata, n_top_genes=150)
#sc.pl.highly_variable_genes(adata)
adata = adata[:, adata.var.highly_variable]


# custom normalization
#adata.raw = adata
#adata.X = adata.X / adata.X.sum(axis=1, keepdims=True) # Divide each element by total


# PCA
sc.tl.pca(adata, svd_solver='arpack')
#sc.pl.pca(adata, color='CD74')
#sc.pl.pca_variance_ratio(adata, log=True)

# Network
sc.pp.neighbors(adata, n_neighbors=50, n_pcs=10)
adata.obs['log_CD3D'] = np.log1p(adata.obs['Mean.CD3'])

# find genes
#adata.var.index[adata.var.index.str.contains('CEACAM8', case=False)]

#UMAP
#sc.tl.tsne(adata)
#sc.pl.tsne(adata, use_raw=False, color=['log_CD3D', 'ITGAX', 'CD163', 'CXCL8', 'CD3D', 'CD4','KRT7'], cmap='vlag')


sc.tl.umap(adata,min_dist=0.4, spread=2.5)
sc.pl.umap(adata, use_raw=True, color=['log_CD3D', 'ITGAX', 'CD163', 'CXCL8', 'CD3D', 'FCGR3A','KRT7','leiden'], cmap='vlag')

#cluster
sc.tl.leiden(adata, resolution=0.3)
sc.pl.umap(adata, color=['leiden', 'imageid'])

#marker genes
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


# real space
plotly (adata,phenotype='leiden',image_id=None,x='CenterX_global_px',y='CenterY_global_px',size=3)

# spatial lag analysisi
adata = sm.tl.spatial_expression (adata, x_coordinate='CenterX_global_px',
                                  y_coordinate='CenterY_global_px',
                                  method='knn', knn=50, imageid='imageid', 
                                  use_raw=True,subset=None,
                                  label='spatial_expression_knn')


# cluster
adata = sm.tl.spatial_cluster (adata, k= 5, method = 'kmeans', df_name='spatial_expression_knn') # results will be saved under adata.obs['spatial_kmeans']

#marker genes
sc.tl.rank_genes_groups(adata, 'spatial_kmeans', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


# view
plotly (adata,phenotype='spatial_kmeans',image_id=None,x='CenterX_global_px',y='CenterY_global_px',size=3)












































