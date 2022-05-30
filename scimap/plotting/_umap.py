#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Sat May 28 22:04:03 2022
# @author: Ajit Johnson
# UMAP plot

"""
!!! abstract "Short Description"
    `sm.pl.umap`: The function allows users to generate a scatter plot of the UMAP. 
## Function
"""


# import lib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import itertools
import matplotlib.colors as colors
import seaborn as sns
import matplotlib.patches as mpatches

# function
def umap (adata, color=None, use_layer=None, use_raw=False, log=False, label='umap',
          cmap='vlag', palette=None, alpha=0.8, figsize=(5, 5), s=None, ncols=None, 
          tight_layout=False, return_data=False, save_figure=None, **kwargs):
    """
Parameters:

    adata : AnnData Object  

    color : list, optional
        Keys for annotations of observations in `adata.obs.columns` or genes in `adata.var.index`. e.g. ['CD3D', 'SOX10']
        The default is None.
        
    use_layer : string, optional  
        Pass name of any `Layer` in AnnData. The default is `None` and `adata.X` is used.
        
    use_raw : bool, optional  
        If set to `True`, values in `adata.raw.X` will be used to color the plot. The default is False.
        
    log : bool, optional  
        If set to `True`, the data will natural log transformed using `np.log1p()` for coloring. The default is False.
        
    label : string, optional  
        The `label key` used when running `sm.tl.umap()`. The default is 'umap'.
        
    cmap : string, optional  
        Color map to use for continous variables. Can be a name or a Colormap 
        instance (e.g. "magma”, "viridis"). The default is 'vlag'.
        
    palette : dict, optional  
        Colors to use for plotting categorical annotation groups. 
        It accepts a `dict` mapping categories to colors. 
        e.g. `palette = {'T cells': '#000000', 'B cells': '#FFF675'}`.
        Auto color will be generated for other categories that are not specified. The default is None.
        
    alpha : float, optional  
        blending value, between 0 (transparent) and 1 (opaque). The default is 0.8.
        
    figsize : tuple, optional  
        Width, height in inches. The default is (10, 10).
        
    s : int, optional  
        The marker size in points. The default is None.
        
    ncols : int, optional  
        Number of panels per row. The default is None.
        
    tight_layout : bool, optional  
        Adjust the padding between and around subplots. If True it will ensure that
        the legends are visible. The default is False.
        
    return_data : bool, optional  
        Returns the data used for plotting. The default is False.
        
    save_figure : string, optional  
        Pass path to saving figure with file extension.
        e.g `\path\to\directory\figure.pdf` The default is None.
    
    **kwargs : Other `matplotlib` parameters. 

Returns:

    final_data : Dataframe
        If return_data is set to `True`.

Example:
```python

# Run UMAP
adata = sm.tl.umap(adata)
    
# plot results
sm.pl.umap(adata, color=['CD3D', 'SOX10'])


```

    """
    
    # check if umap tool has been run
    try:
        adata.obsm[label]
    except KeyError:
        raise KeyError("Please run `sm.tl.umap(adata)` first")
            
    # identify the coordinates
    umap_coordinates = pd.DataFrame(adata.obsm[label],index=adata.obs.index, columns=['umap-1','umap-2'])
    
    # other data that the user requests
    if color is not None:
        if isinstance(color, str):
            color = [color]
        # identify if all elemets of color are available        
        if set(color).issubset(list(adata.var.index) + list(adata.obs.columns)) is False:
            raise ValueError("Element passed to `color` is not found in adata, please check!")
        
        # organise the data
        if any(item in color for item in list(adata.obs.columns)):
            adataobs = adata.obs.loc[:, adata.obs.columns.isin(color)]
        else:
            adataobs = None
            
        if any(item in color for item in list(adata.var.index)):
            # find the index of the marker
            marker_index = np.where(np.isin(list(adata.var.index), color))[0]
            if use_layer is not None:
                adatavar = adata.layers[use_layer][:, np.r_[marker_index]]
            elif use_raw is True:
                adatavar = adata.raw.X[:, np.r_[marker_index]]
            else:
                adatavar = adata.X[:, np.r_[marker_index]]
            adatavar = pd.DataFrame(adatavar, index=adata.obs.index, columns = list(adata.var.index[marker_index]))
        else:
            adatavar = None

        # combine all color data
        if adataobs is not None and adatavar is not None:
            color_data = pd.concat ([adataobs, adatavar], axis=1)
        elif adataobs is not None and adatavar is None:
            color_data = adataobs
        elif adataobs is None and adatavar is not None:
            color_data = adatavar
    else:
        color_data = None
    
    # combine color data with umap coordinates
    if color_data is not None:
        final_data = pd.concat([umap_coordinates, color_data], axis=1)
    else:
        final_data = umap_coordinates
    
    # create some reasonable defaults
    # estimate number of columns in subpolt
    nplots = len(final_data.columns) - 2 # total number of plots
    if ncols is None:
        if nplots >= 4:
            subplot = [math.ceil(nplots / 4), 4]
        elif nplots == 0:
            subplot = [1, 1]
        else:
            subplot = [math.ceil(nplots / nplots), nplots]
    else:
        subplot = [math.ceil(nplots /ncols), ncols]
    
    if nplots == 0:
        n_plots_to_remove = 0
    else:
        n_plots_to_remove = np.prod(subplot) - nplots # figure if we have to remove any subplots
    
    # size of points
    if s is None:
        if nplots == 0:
            s = (100000 / adata.shape[0])
        else:
            s = (100000 / adata.shape[0]) / nplots
    

    # if there are categorical data then assign colors to them
    if final_data.select_dtypes(exclude=["number","bool_","object_"]).shape[1] > 0:
        # find all categories in the dataframe
        cat_data = final_data.select_dtypes(exclude=["number","bool_","object_"])
        # find all categories
        all_cat = []
        for i in cat_data.columns:
            all_cat.append(list(cat_data[i].cat.categories))
        
        # generate colormapping for all categories
        less_9 = [colors.rgb2hex(x) for x in sns.color_palette('Set1')]
        nineto20 = [colors.rgb2hex(x) for x in sns.color_palette('tab20')]
        greater20 = [colors.rgb2hex(x) for x in sns.color_palette('gist_ncar', max([len(i) for i in all_cat]))]
        
        all_cat_colormap = dict()
        for i in range(len(all_cat)):
            if len(all_cat[i]) <= 9:
                dict1 = dict(zip(all_cat[i] , less_9[ : len(all_cat[i]) ]   ))
            elif len(all_cat[i]) > 9 and len(all_cat[i]) <= 20:
                dict1 = dict(zip(all_cat[i] , nineto20[ : len(all_cat[i]) ]   ))
            else:
                dict1 = dict(zip(all_cat[i] , greater20[ : len(all_cat[i]) ]   ))
            all_cat_colormap.update(dict1)
        
        # if user has passed in custom colours update the colors
        if palette is not None:
            all_cat_colormap.update(palette)
    else:
        all_cat_colormap = None
        

    # plot
    fig, ax = plt.subplots(subplot[0],subplot[1], figsize=figsize)
    plt.rcdefaults()
    #plt.rcParams['axes.facecolor'] = 'white'
    
    # remove unwanted axes
    #fig.delaxes(ax[-1])
    if n_plots_to_remove > 0:
        for i in range(n_plots_to_remove):
            fig.delaxes(ax[-1][  (len(ax[-1])-1)-i :   (len(ax[-1]))-i   ][0])
    
    # to make sure the ax is always 2x2
    if any(i > 1 for i in subplot):
        if any(i == 1 for i in subplot):
            ax = ax.reshape(subplot[0],subplot[1])
    
    if nplots == 0:
        ax.scatter(x = final_data['umap-1'], y = final_data['umap-2'], s=s, cmap=cmap, alpha=alpha, **kwargs)
        plt.xlabel("UMAP-1"); plt.ylabel("UMAP-2") 
        plt.tick_params(right= False,top= False,left= False, bottom= False)
        ax.get_xaxis().set_ticks([]); ax.get_yaxis().set_ticks([])
        if tight_layout is True:
            plt.tight_layout()
    
    elif all(i == 1 for i in subplot):
        column_to_plot = [e for e in list(final_data.columns) if e not in ('umap-1', 'umap-2')][0]
        if all_cat_colormap is None:
            im = ax.scatter(x = final_data['umap-1'], y = final_data['umap-2'], s=s, 
                       c=final_data[column_to_plot],
                       cmap=cmap, alpha=alpha, **kwargs)
            plt.colorbar(im, ax=ax)
        else:
            ax.scatter(x = final_data['umap-1'], y = final_data['umap-2'], s=s, 
                       c=final_data[column_to_plot].map(all_cat_colormap),
                       cmap=cmap, alpha=alpha, **kwargs)
            # create legend
            patchList = []
            for key in list(final_data[column_to_plot].unique()):
                data_key = mpatches.Patch(color=all_cat_colormap[key], label=key)
                patchList.append(data_key)    
                ax.legend(handles=patchList,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        
        plt.xlabel("UMAP-1"); plt.ylabel("UMAP-2") 
        plt.title(column_to_plot)
        plt.tick_params(right= False,top= False,left= False, bottom= False)
        ax.set(xticklabels = ([])); ax.set(yticklabels = ([]))
        if tight_layout is True:
            plt.tight_layout()
        
    else: 
        column_to_plot = [e for e in list(final_data.columns) if e not in ('umap-1', 'umap-2')]
        k=0
        for i,j in itertools.product(range(subplot[0]), range(subplot[1])):
            
            if final_data[column_to_plot[k]].dtype == 'category':
                ax[i,j].scatter(x = final_data['umap-1'], y = final_data['umap-2'], s=s,
                                c=final_data[column_to_plot[k]].map(all_cat_colormap),
                                cmap=cmap, alpha=alpha, **kwargs)
                # create legend
                patchList = []
                for key in list(final_data[column_to_plot[k]].unique()):
                    data_key = mpatches.Patch(color=all_cat_colormap[key], label=key)
                    patchList.append(data_key)    
                    ax[i,j].legend(handles=patchList,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            else:               
                im = ax[i,j].scatter(x = final_data['umap-1'], y = final_data['umap-2'], s=s,
                                c=final_data[column_to_plot[k]],
                                cmap=cmap, alpha=alpha, **kwargs)
                plt.colorbar(im, ax=ax[i, j])

            ax[i,j].tick_params(right= False,top= False,left= False, bottom= False)
            ax[i,j].set_xticklabels([]); ax[i,j].set_yticklabels([])
            ax[i,j].set_xlabel("UMAP-1"); ax[i,j].set_ylabel("UMAP-2")
            ax[i,j].set_title(column_to_plot[k])
            if tight_layout is True:
                plt.tight_layout()
            k = k+1 # iterator
            
    # if save figure is requested
    if save_figure is not None:
        plt.savefig(save_figure)  
    
    # return data if needed
    if return_data is True:
        return final_data
            
    
    