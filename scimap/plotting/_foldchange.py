#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Wed Apr  7 17:46:29 2021
# @author: Ajit Johnson Nirmal
"""
!!! abstract "Short Description"
    `sm.pl.foldchange`: The Function allows users to visualize foldchange in abundance of celltypes between samples/ROI's. 
    Run `sm.tl.foldchange` first to compute the foldchange.

## Function
"""

# lib
import seaborn as sns; sns.set(color_codes=True)
import matplotlib
import matplotlib.pyplot as plt 
import numpy as np
from pandas.plotting import parallel_coordinates
sns.set_style("white")

# Function
def foldchange (adata, label='foldchange', 
                p_val=0.05, nonsig_color='grey',subset_xaxis=None,subset_yaxis=None,
                cmap = 'vlag', log=True,center=0, 
                method='heatmap', invert_axis=None,
                parallel_coordinates_color=None,matplotlib_bbox_to_anchor=(1.04,1),
                matplotlib_legend_loc='upper left',xticks_rotation=90,
                return_data = False,
                **kwargs):
    """
Parameters:
    adata : Anndata object

    label : strong, optional  
        label used when running `sm.tl.foldchange`.

    p_val : float, optional  
        p_val cut-off above which is considered not-significant. The cells containing
        non-significant changes will be highlighted in the heatmap.

    nonsig_color : string, optional  
        Color used to highlight non-significant fold changes in the heatmap.

    subset_xaxis : list, optional  
        Subset x-axis before plotting. Pass in a list of categories. eg- subset_xaxis = ['CelltypeA', 'CellTypeB']. 

    subset_yaxis : list, optional  
        Subset y-axis before plotting. Pass in a list of categories. eg- subset_yaxis = ['ROI_1', 'ROI_5']. 

    cmap : string, optional  
        Color map. Can be a name or a Colormap instance (e.g. 'magma', 'viridis').

    log : bool, optional  
        Convert foldchange to log2 scale.

    center : float, optional  
        The center value to be used in heatmap.

    method : string, optional  
        Two methods are available for plotting the foldchanges  
        a) Heatmap: Use `heatmap`  
        b) parallel coordinates plot : Use `parallel_coordinates`  

    invert_axis : bool, optional  
        Flip the axis of the plot.

    parallel_coordinates_color : list, optional  
        Custom colors for each category.

    matplotlib_bbox_to_anchor : tuple, optional  
        Bounding box argument used along with matplotlib_legend_loc to control
        the legend location when using the matplotlib method.

    matplotlib_legend_loc : TYPE, optional  
        Location of legend used along with matplotlib_bbox_to_anchor to control
        the legend location when using the matplotlib method.

    xticks_rotation : int, optional  
        Angle the x-axis ticks.

    return_data: bool, optional  
        Return the final data used for plotting.

    **kwargs : Additional keyword arguments passed to:  
        a) sns.clustermap  
        b) pandas.parallel_coordinates  

Returns:
    Plot:  
        Data used for the plot if `return_data = True`
    
Example:
```python
    # Heatmap of foldchnage  
    sm.pl.foldchange (adata, label='foldchange', method='heatmap',
                     p_val=0.05, nonsig_color='grey',
                     cmap = 'vlag', log=True, center=0, linecolor='black',linewidths=0.7,
                     vmin=-5, vmax=5, row_cluster=False)
    
    # Parallel_coordinates plot of the foldchanges
    foldchange (adata, label='foldchange', 
                log=True, method='parallel_coordinates', invert_axis=True,
                parallel_coordinates_color=['black','blue','green','red','#000000'],
                matplotlib_bbox_to_anchor=(1.04,1),
                matplotlib_legend_loc='upper left',
                xticks_rotation=90,
                return_data = False
```
    """
        
    # set color for heatmap
    #cmap_updated = copy.copy(matplotlib.cm.get_cmap(cmap))
    cmap_updated = matplotlib.cm.get_cmap(cmap)
    cmap_updated.set_bad(color=nonsig_color)
    
    
    # get the data
    fc = adata.uns[str(label)+'_fc']
    p = adata.uns[str(label)+'_pval']
        
    #fold
    fold = fc.copy()
    p_mask = p.copy()
    
    # reference image
    ref = fold.index.name
       
    # log
    if log is True:
        fold = np.log2(fold)
    
    # create a mask for non-sig values
    p_mask[p_mask > p_val] = np.nan
    
    # subset x axis data
    if subset_xaxis is not None:
        if isinstance (subset_xaxis, str):
            subset_xaxis = [subset_xaxis]
        fold = fold [subset_xaxis]
        p_mask = p_mask [subset_xaxis]
        #reorder
        
    # subset y axis data
    if subset_yaxis is not None:
        if isinstance (subset_yaxis, str):
            subset_yaxis = [subset_yaxis]
        fold = fold.loc [subset_yaxis]
        p_mask = p_mask.loc [subset_yaxis]
        #reorder
        
    # invert axis if user requests
    if invert_axis is True:
        fold = fold.T
        p_mask = p_mask.T
    
    #mask
    mask = p_mask.isnull() # identify the NAN's for masking 
    
    if method == 'heatmap':
        # heatmap of the foldchange
        #g= sns.clustermap(fold, cmap=cmap, mask=mask, center=center, col_cluster=False, row_cluster=False)
        g= sns.clustermap(fold, cmap=cmap, mask=mask, center=center, **kwargs)
        plt.suptitle('reference: '+ str(ref))
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=xticks_rotation)
        plt.tight_layout()
        
    
    if method == 'parallel_coordinates':
        fold['sample'] = fold.index
        # plotting
        fig, axes = plt.subplots()
        if parallel_coordinates_color is not None:
            parallel_coordinates(fold, 'sample', color=parallel_coordinates_color, **kwargs)
        else:
            #parallel_coordinates(fold, 'sample', colormap=cmap_updated)
            parallel_coordinates(fold, 'sample', colormap=cmap_updated, **kwargs)
        axes.grid(False)
        plt.legend(bbox_to_anchor=matplotlib_bbox_to_anchor, loc=matplotlib_legend_loc)
        plt.axhline(y=0, color='black', linestyle='-')
        plt.xticks(rotation = xticks_rotation)
        plt.suptitle('reference: '+ str(ref))
        fig.tight_layout()
        
    # return data
    if return_data is True:
        return fold
    