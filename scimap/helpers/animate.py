#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Sat May 28 17:31:18 2022
# @author: Ajit Johnson Nirmal
# Animation in matplotlib

"""
!!! abstract "Short Description"
    `sm.hl.animate`: This function creates dynamic animations transitioning between 
    UMAP embeddings and physical X and Y coordinates, offering a visual bridge between 
    abstract and physical data representations. Given varying computer configurations and 
    execution environments (such as Jupyter notebooks), real-time animation playback may 
    experience performance issues or may not display properly. Therefore, it is strongly advised 
    to save animations to disk. Note that saving requires `imagemagick` installed on your system. 
    Installation instructions for `imagemagick` can be found 
    at: https://imagemagick.org/script/download.php

## Function
"""


# libs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd
import matplotlib.colors as colors
import seaborn as sns
import matplotlib.patches as mpatches


# function
def animate (adata, color=None,
             palette=None,
             embedding='umap', 
             x_coordinate='X_centroid', 
             y_coordinate='Y_centroid',
             flip_y=True,
             imageid='imageid', 
             subset=None,
             layer=None, 
             use_raw=False, 
             log=False,
             subsample=None,
             random_state=0,
             n_frames=50, 
             interval=50,
             reverse=True,
             final_frame=5, 
             s=None, 
             alpha=1,  
             cmap='vlag',
             tight_layout=True,
             plot_legend=False,
             title=None, 
             fontsize=20,
             watermark=True,
             figsize=(5,5), 
             pltStyle=None,
             verbose=True,
             save_animation=None,**kwargs):
    """
Parameters:
    adata (anndata.AnnData):   
        The annotated data object to be visualized.
    
    color (list, optional):  
        Identifiers for annotations or genes to color the animation. Accepts a single annotation name. 
    
    palette (dict, optional):  
        Custom color mapping for categorical annotations. Unspecified categories are automatically colored. 
    
    embedding (str, optional):  
        Specifies the UMAP embedding label in `adata.obsm`. 
    
    x_coordinate, y_coordinate (str, optional):  
        Columns in `adata.obs` for spatial coordinates. Defaults are 'X_centroid' and 'Y_centroid', respectively.
    
    flip_y (bool, optional):  
        Whether to invert the Y-axis, useful if the Y-coordinates appear flipped. 
    
    imageid (str, optional):  
        Column name in `adata.obs` that identifies unique images, for datasets containing multiple images. 
    
    subset (list, optional):  
        Identifiers for specific images to visualize, used in conjunction with `imageid`.
    
    layer (str, optional):  
        Specifies a layer in `adata.layers` to use for the animation. Default is None, using `adata.X`.
    
    use_raw (bool, optional):  
        Whether to use data from `adata.raw.X` for coloring. 
    
    log (bool, optional):  
        Applies natural log transformation to the data if True. 
    
    subsample (float, optional):  
        Fraction of data to randomly subsample for large datasets, between 0-1. 
    
    random_state (int, optional):  
        Seed for the random number generator, ensuring reproducibility. 
    
    n_frames (int, optional):  
        Number of frames between UMAP and spatial coordinates, affecting animation smoothness. 
    
    interval (int, optional):  
        Time interval between frames in milliseconds. 
    
    reverse (bool, optional):  
        If True, includes a reverse transition from UMAP to spatial coordinates. 
    
    final_frame (int, optional):  
        Number of frames to display the final frame, enhancing visualization. 
    
    s (int, optional):  
        Marker size in points. 
    
    alpha (float, optional):  
        Opacity of markers, between 0 (transparent) and 1 (opaque). 
    
    cmap (str, optional):  
        Colormap for continuous variables. 
    
    tight_layout (bool, optional):  
        Adjusts subplot padding, ensuring visibility of legends. 
    
    plot_legend (bool, optional):  
        Whether to display the legend. 
    
    title (bool or str, optional):  
        Adds a title to the plot. Custom titles can be specified. 
    
    fontsize (int, optional):  
        Font size for the title. 
    
    watermark (bool, optional):  
        Displays a 'made with scimap' watermark. 
    
    figsize (tuple, optional):  
        Figure dimensions in inches. 
    
    pltStyle (str, optional):  
        Matplotlib plot style to use. 
    
    save_animation (str, optional):  
        File path to save the animation. Saving is recommended for optimal viewing. 
    
    **kwargs : Additional matplotlib parameters.
 

Returns:
        Animation:  
            An interactive animation or saved file illustrating the dynamic transition.

Example:
    ```python
    
    # Run UMAP
    adata = sm.tl.umap(adata)
    
    # Basic animation with default UMAP and spatial coordinates
    sm.hl.animate(adata)

    # Customized animation with specific cell-type coloring and reverse transition
    sm.hl.animate(adata, color='cell_type', reverse=False, save_animation='umap_to_spatial.gif')

    # Animation with a subset of images, using custom palette and increased frame interval
    sm.hl.animate(adata, color='condition', palette={'Control': '#1f77b4', 'Treated': '#ff7f0e'}, 
            subset='image_01', interval=100, save_animation='custom_animation.gif')
    
    ```

    """
    
    

    # intrapolation function between co-ordinate sytems
    def tween(e1, e2, n_frames, final_frame):
        
        # number of frame to pop
        #n_frames = int(n_frames + (n_frames*0.3))
        for i in range(5):
            yield e1
        for i in range(n_frames):
            alpha = i / float(n_frames - 1)
            yield (1 - alpha) * e1 + alpha * e2
        for i in range(final_frame):
            yield e2
            
        return
    
        
    # check if umap tool has been run
    try:
        adata.obsm[embedding]
    except KeyError:
        raise KeyError("Please run `sm.tl.umap(adata)` first")
    
    # identify the coordinates
    umap_coordinates = pd.DataFrame(adata.obsm[embedding],index=adata.obs.index, columns=['umap-1','umap-2'])
    real_coordinates = adata.obs[[x_coordinate,y_coordinate]]
    
    # other data that the user requests
    if color is not None:
        if isinstance(color, str):
            color = [color]
            
        # identify if all elemets of color are available        
        if len(color) > 1:
            raise ValueError("Only a single value in `color` is supported")
            
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
            if layer is not None:
                adatavar = adata.layers[layer][:, np.r_[marker_index]]
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
            # convert to string
            color_data[color] = color_data[color].astype('category')
        elif adataobs is None and adatavar is not None:
            color_data = adatavar    
        
    else:
        color_data = None
    
    # combine color data with umap coordinates
    if color_data is not None:
        final_data = pd.concat([umap_coordinates, real_coordinates, color_data], axis=1)
    else:
        final_data = umap_coordinates
    
    # subset the final data if nedded
    if subset is not None:
        if isinstance(subset, str):
            subset = [subset]
        cell_to_keep = adata[adata.obs[imageid].isin(subset)].obs.index
        final_data = final_data.loc[cell_to_keep]
    
    # subsample the data if user requests
    if subsample is not None:
        final_data = final_data.sample(frac=subsample, replace=False, random_state=random_state)
    
    # extract the spaces
    e1 = final_data[['umap-1', 'umap-2']].values.astype(float)
    e2 = final_data[[x_coordinate,y_coordinate]].values.astype(float)


    # rescale to same co-ordinates system
    e1[:, 0] -= (max(e1[:, 0]) + min(e1[:, 0])) / 2
    e1[:, 1] -= (max(e1[:, 1]) + min(e1[:, 1])) / 2
    # scale
    scale = max(max(e1[:, 0]) - min(e1[:, 0]), max(e1[:, 1]) - min(e1[:, 1]))
    e1[:, 0] /= scale
    e1[:, 1] /= scale
    # Translate
    e1[:, 0] += 0.5
    e1[:, 1] += 0.5
    
    # rescale co-ordinates
    e2[:, 0] -= (max(e2[:, 0]) + min(e2[:, 0])) / 2
    e2[:, 1] -= (max(e2[:, 1]) + min(e2[:, 1])) / 2
    # scale
    scale = max(max(e2[:, 0]) - min(e2[:, 0]), max(e2[:, 1]) - min(e2[:, 1]))
    e2[:, 0] /= scale
    e2[:, 1] /= scale
    # Translate
    e2[:, 0] += 0.5
    e2[:, 1] += 0.5
    
    # remove the identified indeces
    def delete_multiple_element(list_object, indices):
        indices = sorted(indices, reverse=True)
        for idx in indices:
            if idx < len(list_object):
                list_object.pop(idx)
    
    # run the interpolation
    interpolation = list(tween(e1, e2, n_frames=n_frames, final_frame=final_frame))
    # drop x number of frames
    top_frames = int(n_frames + 5)
    
    l = np.percentile(range(5,top_frames),30); h = np.percentile(range(5,top_frames),80)
    index_between = list(range(int(l), int(h)))
    numElems = int(len(index_between) * 0.5)
    drop = np.round(np.linspace(0, len(index_between) - 1, numElems)).astype(int)
    drop_index = [index_between[i] for i in drop] 
    
    # delete frames
    delete_multiple_element(interpolation, drop_index)
    
    top20 = np.percentile(range(5,top_frames),20); top30 = np.percentile(range(5,top_frames),30)
    bottom80 = np.percentile(range(5,top_frames),80); bottom90 = np.percentile(range(5,top_frames),90)
    
    ib_top = list(range(int(top20), int(top30)))
    ib_bottom = list(range(int(bottom80), int(bottom90)))
    ib = ib_top + ib_bottom
    numElems2 = int(len(ib) * 0.20)
    drop2 = np.round(np.linspace(0, len(ib) - 1, numElems2)).astype(int)
    di = [ib[i] for i in drop2] 
    # delete frames
    delete_multiple_element(interpolation, di)
    
    top10 = np.percentile(range(5,top_frames),10); top19 = np.percentile(range(5,top_frames),19)
    bottom91 = np.percentile(range(5,top_frames),91); bottom95 = np.percentile(range(5,top_frames),95)
    
    ib_top = list(range(int(top10), int(top19)))
    ib_bottom = list(range(int(bottom91), int(bottom95)))
    ib = ib_top + ib_bottom
    numElems2 = int(len(ib) * 0.10)
    drop2 = np.round(np.linspace(0, len(ib) - 1, numElems2)).astype(int)
    di = [ib[i] for i in drop2] 
    # delete frames
    delete_multiple_element(interpolation, di)
    
    


    if reverse is True:
        interpolation = interpolation + interpolation[::-1]
    
    # generate colors
    if s is None:
        s = 130000 / final_data.shape[0]
    
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
        
    # number of plots
    nplots = len(final_data.columns) - 4 # total number of plots
    if nplots > 0:
        column_to_plot = [e for e in list(final_data.columns) if e not in ('umap-1', 'umap-2',x_coordinate,y_coordinate)][0]
        if all_cat_colormap is not None:
            custom_color = list(final_data[column_to_plot].map(all_cat_colormap).values)


    # plot
    plt.rcdefaults()
    if pltStyle is not None:
        plt.style.use(pltStyle)
    fig, ax = plt.subplots(figsize=figsize)
    
    
    ax.set(xlim=(-0.1, 1.1), ylim=(-0.1, 1.1))
    if flip_y is True:
        ax.invert_yaxis()
    
    
    
    if nplots == 0:
        scat = ax.scatter(x = interpolation[0][:, 0], y = interpolation[0][:, 1], s=s, cmap=cmap, alpha=alpha, **kwargs)
        plt.tick_params(right= False,top= False,left= False, bottom= False)
        ax.get_xaxis().set_ticks([]); ax.get_yaxis().set_ticks([])
        if watermark is True:
            ax.text(1.08, 1.08, "made with scimap.xyz",horizontalalignment="right",
            verticalalignment="bottom", alpha=0.5,fontsize=fontsize * 0.4)
        if title is True: 
            plt.title(column_to_plot, fontsize=fontsize)
        elif isinstance(title, str):
            plt.title(title, fontsize=fontsize)  
        if tight_layout is True:
            plt.tight_layout()
    
    if nplots > 0:
        if all_cat_colormap is None:
            scat = ax.scatter(x = interpolation[0][:, 0], y = interpolation[0][:, 1], s=s, 
                           c=final_data[column_to_plot],
                           cmap=cmap, alpha=alpha, **kwargs)
            if plot_legend is True:
                plt.colorbar(scat, ax=ax)
        else:
            scat = ax.scatter(x = interpolation[0][:, 0], y = interpolation[0][:, 1], s=s, 
                           c=custom_color,
                           cmap=cmap, alpha=alpha, **kwargs)
            # create legend
            if plot_legend is True:
                patchList = []
                for key in list(final_data[column_to_plot].unique()):
                    data_key = mpatches.Patch(color=all_cat_colormap[key], label=key)
                    patchList.append(data_key)    
                    ax.legend(handles=patchList,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    
        if title is True: 
            plt.title(column_to_plot, fontsize=fontsize)
        elif isinstance(title, str):
            plt.title(title, fontsize=fontsize) 
        if watermark is True:
            ax.text(1.08, 1.08, "made with scimap.xyz",horizontalalignment="right",
            verticalalignment="bottom", alpha=0.5,fontsize=fontsize * 0.4)
        plt.tick_params(right= False,top= False,left= False, bottom= False)
        ax.set(xticklabels = ([])); ax.set(yticklabels = ([]))
        if tight_layout is True:
            plt.tight_layout()
        

 
    def animate(i):
        scat.set_offsets(interpolation[i])
        
    anim = FuncAnimation(fig, animate, interval=interval, frames=len(interpolation)-1)
    
    
     
    if save_animation is not None:
        if verbose:
            print ('Saving file- This can take several minutes to hours for large files')
        anim.save( save_animation + '_scimap.gif', writer='imagemagick', fps=24)

    # save animation
    #anim.save('/Users/aj/Downloads/filename.mp4')
    
    return plt.show(anim, block=False)
    
    
