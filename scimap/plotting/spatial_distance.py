#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Sat Oct 24 13:07:38 2020
# @author: Ajit Johnson Nirmal
"""
!!! abstract "Short Description"
    `sm.pl.spatial_distance`: This function enables the visualization of the average 
    shortest distances between selected phenotypes or cell types, offering insights 
    into spatial relationships within biological samples. To accurately generate 
    these visual representations, it's essential to first compute spatial distances 
    using `sm.tl.spatial_distance`. This preparatory step ensures the data necessary 
    for creating comprehensive heatmaps, numeric comparisons, and distribution plots is 
    available, facilitating a deeper understanding of spatial patterning and interactions 
    among cell populations.

## Function
"""

# library
import pandas as pd
import matplotlib
import numpy as np
import seaborn as sns; sns.set(color_codes=True)
sns.set_style("white")

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42


def spatial_distance (adata, 
                      spatial_distance='spatial_distance',
                      phenotype='phenotype',
                      imageid='imageid',
                      log=False,
                      method='heatmap',
                      heatmap_summarize=True,
                      heatmap_na_color='grey',
                      heatmap_cmap='vlag_r',
                      heatmap_row_cluster=False,
                      heatmap_col_cluster=False,
                      heatmap_standard_scale=0,
                      distance_from=None,
                      distance_to=None,
                      x_axis = None,
                      y_axis = None,
                      facet_by = None,
                      plot_type = None,
                      return_data = False, 
                      subset_col=None, 
                      subset_value=None,
                      **kwargs):
    """
Parameters:
        adata (anndata.AnnData):  
            The annotated data matrix with spatial distance calculations.

        spatial_distance (str, optional):  
            Key in `adata.uns` where spatial distance data is stored, typically the output of `sm.tl.spatial_distance`.

        phenotype (str):  
            Column in `adata.obs` containing phenotype or cell type annotations.

        imageid (str, optional):  
            Column in `adata.obs` identifying different images or samples.

        log (bool, optional):  
            If True, applies log transformation to the distance data.

        method (str, optional):  
            Visualization method: 'heatmap', 'numeric', or 'distribution'.

        heatmap_summarize (bool, optional):  
            If True, summarizes distances across all images or samples for the heatmap.

        heatmap_na_color (str, optional):  
            Color for NA values in the heatmap.

        heatmap_cmap (str, optional):  
            Colormap for the heatmap.

        heatmap_row_cluster, heatmap_col_cluster (bool, optional):  
            If True, clusters rows or columns in the heatmap.

        heatmap_standard_scale (int, optional):  
            Standardizes rows (0) or columns (1) in the heatmap.

        distance_from, distance_to (str, optional):  
            Phenotypes of interest for distance calculation in 'numeric' or 'distribution' plots.

        x_axis, y_axis (str, optional):  
            Axes labels for 'numeric' or 'distribution' plots.

        facet_by (str, optional):  
            Categorizes plots into subplots based on this column.

        plot_type (str, optional):  
            For 'numeric' plots: options include 'box', 'violin', etc. For 'distribution' plots: 'hist', 'kde', etc.

        subset_col (str, optional):  
            Column name for subsetting data before plotting.

        subset_value (list, optional):  
            Values in `subset_col` to include in the plot.

        **kwargs:  
            Additional keyword arguments for plotting functions.

Returns:
    Plot and dataFrame (matplotlib, pandasDF):
        If `return_data` is True, returns the data frame used for plotting; otherwise, displays the plot.

Example:
    ```python
    
    # Generate a heatmap of spatial distances
    sm.pl.spatial_distance(adata, method='heatmap', phenotype='cell_type', imageid='sample_id')

    # Numeric plot showing distance from one phenotype to all others
    sm.pl.spatial_distance(adata, method='numeric', distance_from='Tumor', phenotype='cell_type', plot_type='boxen')

    # Distribution plot comparing distances between two specific phenotypes
    sm.pl.spatial_distance(adata, method='distribution', distance_from='Tumor', distance_to='Stroma', 
                     plot_type='kde', x_axis='distance', y_axis='group')
    
    ```
    """


    # set color for heatmap
    cmap_updated = matplotlib.cm.get_cmap(heatmap_cmap)
    cmap_updated.set_bad(color=heatmap_na_color)


    # Copy the spatial_distance results from anndata object
    try:
        diatance_map = adata.uns[spatial_distance].copy()
    except KeyError:
        raise ValueError('spatial_distance not found- Please run sm.tl.spatial_distance first')

    # subset the data if user requests
    if subset_col is not None:
        if isinstance(subset_value, str):
            subset_value = [subset_value]
        # find the cell names to be subsetted out
        obs = adata.obs[[subset_col]]
        cells_to_subset = obs[obs[subset_col].isin(subset_value)].index

        # subset the diatance_map
        diatance_map = diatance_map.loc[diatance_map.index.intersection(cells_to_subset)]
        #diatance_map = diatance_map.loc[cells_to_subset]


    # Convert distance to log scale if user requests
    if log is True:
        diatance_map = np.log1p(diatance_map)

    # Method
    if method=='heatmap':
        if heatmap_summarize is True:
            # create the necessary data
            data = pd.DataFrame({'phenotype': adata.obs[phenotype]})
            data = pd.merge(data, diatance_map, how='outer',left_index=True, right_index=True) # merge with the distance map
            k = data.groupby(['phenotype']).mean() # collapse the whole dataset into mean expression
            d = k[k.index]
        else:
            # create new naming scheme for the phenotypes
            non_summary = pd.DataFrame({'imageid': adata.obs[imageid], 'phenotype': adata.obs[phenotype]})
            non_summary['imageid'] = non_summary['imageid'].astype(str) # convert the column to string
            non_summary['phenotype'] = non_summary['phenotype'].astype(str) # convert the column to string
            non_summary['image_phenotype'] = non_summary['imageid'].str.cat(non_summary['phenotype'],sep="_")
            # Merge distance map with phenotype
            data = pd.DataFrame(non_summary[['image_phenotype']])
            data = pd.merge(data, diatance_map, how='outer',left_index=True, right_index=True)
            k = data.groupby(['image_phenotype']).mean()
            d = k.sort_index(axis=1)
        # Generate the heatmap
        mask = d.isnull() # identify the NAN's for masking 
        d = d.fillna(0) # replace nan's with 0 so that clustering will work
        # Heatmap
        sns.clustermap(d, cmap=heatmap_cmap, row_cluster=heatmap_row_cluster,
                       col_cluster=heatmap_col_cluster, mask=mask,
                       standard_scale=heatmap_standard_scale, **kwargs)
    else:

        # condition-1
        if distance_from is None and distance_to is None:
            raise ValueError('Please include distance_from and/or distance_to parameters to use this method')

        # condition-2
        if distance_from is None and distance_to is not None:
            raise ValueError('Please `distance_from` parameters to use this method')

        # condition-3
        if distance_to is not None:
            # convert input to list if needed
            if isinstance(distance_to, str):
                distance_to = [distance_to]

        # Start
        pheno_df = pd.DataFrame({'imageid': adata.obs[imageid], 'phenotype': adata.obs[phenotype]}) #image id and phenotype
        data = pd.merge(pheno_df, diatance_map, how='outer',left_index=True, right_index=True) # merge with the distance map
        data = data[data['phenotype'] == distance_from] # subset the pheno of interest

        if distance_to is not None:
            data = data[distance_to] # drop columns that are not requested in distance_to
        else:
            data = data.drop(['phenotype','imageid'], axis=1) # drop the phenotype column before stacking

        d = data.stack().reset_index() # collapse everything to one column
        d.columns = ['cellid', 'group', 'distance']
        d = pd.merge(d, pheno_df, left_on='cellid', right_index=True) # bring back the imageid and phenotype

        # Convert columns to str
        for col in ['imageid', 'group','phenotype']:
            d[col] = d[col].astype(str)

        # Convert columns to categorical so that it drops unused categories
        for col in ['imageid', 'group','phenotype']:
            d[col] = d[col].astype('category')

        # re arrange the order based on from and to list provided
        if distance_to is not None:
            d['group'] = d['group'].cat.reorder_categories(distance_to)
            d = d.sort_values('group')

        # Plotting
        if method=='numeric':
            if x_axis is None and y_axis is None and facet_by is None and plot_type is None:
                sns.catplot(data=d, x="distance", y="group", col="imageid", kind="boxen", **kwargs)
            else:
                sns.catplot(data=d, x=x_axis, y=y_axis, col=facet_by, kind=plot_type, **kwargs)

        if method=='distribution':
            if x_axis is None and y_axis is None and facet_by is None and plot_type is None:
                sns.displot(data=d, x="distance", hue="imageid",  col="group", kind="kde", **kwargs)
            else:
                sns.displot(data=d, x=x_axis, hue=y_axis, col=facet_by, kind=plot_type,**kwargs)

    # return
    if return_data is True:
        return d
