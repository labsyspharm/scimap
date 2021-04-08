#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 22:05:14 2021
@author: Ajit Johnson Nirmal
Pie plot for categories of interest.
A pie plot is a proportional representation of the numerical data in a column. 
"""

# Lib
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np





# Function
def pie (adata, pie_column='imageid', group_by=None, ncols=None, 
         label='auto', title='auto', colors=None, autopct='%1.1f%%',
         wedgeprops = {'linewidth': 0},**kwargs):
    
    # create copy of the required data
    data = adata.obs
    
    # calculate the proportion
    if group_by is None:
        prop = data[pie_column].value_counts().reset_index(inplace=False)
        prop.columns = [pie_column, 'value']
        prop['group_by'] = pie_column
        labels = np.unique(prop[pie_column])

    else:
        # if group_by is provided
        prop = pd.DataFrame(data.groupby([group_by,pie_column]).size()).reset_index(inplace=False)
        prop.columns = ['group_by',pie_column,'value']
        labels = np.unique(prop[pie_column])

        #
        if ncols is not None:
            g = prop.groupby('group_by')
            rows = int(np.ceil(len(g)/ncols))
        else:
            g = prop.groupby('group_by')
            rows = 1
            ncols = len(g)
    
    # remove label if requested 
    if label == 'auto':
        label = labels
    elif label is None:
        label = None
    else:
        label = label
    
    
    # plot
    if group_by is None:
        fig, ax = plt.subplots()
        ax.pie(prop.value, labels=label,colors=colors, wedgeprops = {'linewidth': 0},**kwargs)
        if title is None:
            pass
        else:
            ax.set_title(pie_column)
    else:
        # plot the figure
        fig, axes = plt.subplots(ncols=ncols, nrows=rows)
        for (c, grp), ax in zip(g, axes.flat):
            ax.pie(grp.value, labels=label, colors=colors, wedgeprops = {'linewidth': 0},**kwargs)
            if title is None:
                pass
            else:
                ax.set_title(c)     
    plt.show()
    
    
    
    
    
