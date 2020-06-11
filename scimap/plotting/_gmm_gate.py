#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 19:06:40 2020
@author: Ajit Johnson Nirmal
Manual Gate finder
"""

import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture 
from scipy import stats
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt



def gmm_gate (adata, marker_of_interest):
    """
    

    Parameters
    ----------
    adata : Ann Data Object
    marker_of_interest : string
        Marker of interest.

    Returns
    -------
    Distribution plot of the marker of interest along with GMM overlaid.
    
    Example
    -------
    sm.pl.gmm_dist_plot (adata, marker_of_interest='CD45')

    """
    
    # If no raw data is available make a copy
    if adata.raw is None:
        adata.raw = adata
    
    # Copy of the raw data if it exisits
    if adata.raw is not None:
        adata.X = adata.raw.X
    
    # Make a copy of the data with the marker of interest
    data = pd.DataFrame(np.log1p(adata.X), columns = adata.var.index, index= adata.obs.index)[[marker_of_interest]]
    
    # Clip of the top and bottom outliers before applying the model
    data = data.clip(lower =np.percentile(data,1), upper=np.percentile(data,99))
    
    # Apply Gaussian mixture model
    m = data[marker_of_interest].values
    data_gm = m.reshape(-1, 1)
    gmm = GaussianMixture(n_components=2)
    gmm.fit(data_gm)
    gate = np.mean(gmm.means_)
    
    mean = gmm.means_
    cov = gmm.covariances_
    std = [ np.sqrt(  np.trace(cov[i])/2) for i in range(0,2) ]
    #std = np.sqrt(np.trace(cov)/2) 
    
    
    # Merge with image ID
    if len(adata.obs['ImageId'].unique()) > 1:
        # Identify the most and least positivity of the given marker of interest
        # Generate a dataframe with various gates
        dd = data.values
        dd = np.where(dd < gate, np.nan, dd)
        np.warnings.filterwarnings('ignore')
        dd = np.where(dd > gate, 1, dd)
        dd = pd.DataFrame(dd, index = data.index, columns = ['gate'])
        dd = dd.merge(pd.DataFrame(adata.obs['ImageId']), how='outer', left_index=True, right_index=True)
        dd = dd[dd['gate'] == 1]
        image_positivity = dd.groupby('ImageId').count().sort_values('gate')
        # Create a string for overlaying on image
        low = list(image_positivity.head(3).index)
        low = ','.join(low)
        high = list(image_positivity.tail(3).index)
        high = ','.join(high)
        # Dor image
        textstr = '\n'.join((
            ("Most positive in: ", 
             high,
             "Least positive in ",
             low
            )))

    # Plot
    sns.set_style("white")
    fig, ax = plt.subplots()
    x = np.linspace(mean[0] - 3*std[0], mean[0] + 3*std[0], 1000)
    y = np.linspace(mean[1] - 3*std[1], mean[1] + 3*std[1], 1000)
    sns.distplot(data[marker_of_interest], color = 'grey')
    plt.axvline(x= mean[0], c='#000000', linestyle='dotted')
    plt.axvline(x= mean[1], c='#000000', linestyle='dotted')
    plt.plot(x, stats.norm.pdf(x, mean[0], std[0]), linestyle='dashed')
    plt.plot(y, stats.norm.pdf(y, mean[1], std[1]), linestyle='dashed')
    plt.xticks(np.arange(min(data[marker_of_interest])-1, max(data[marker_of_interest])+1, 0.5))
    plt.title(marker_of_interest, fontsize=20)    
    if len(adata.obs['ImageId'].unique()) > 1:
        ax.text(0.05, 0.95, textstr, fontsize=14, transform=ax.transAxes,
            verticalalignment='top')


