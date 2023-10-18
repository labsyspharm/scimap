# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 14:17:09 2019
@author: Ajit Johnson Nirmal
Program: Lasso Selector for selecting regions of interest
"""

# Load packages
from __future__ import print_function
#from six.moves import input
import numpy as np
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'auto')
import matplotlib.pyplot as plt

class SelectFromCollection(object):
    """Select indices from a matplotlib collection using `LassoSelector`.

    Selected indices are saved in the `ind` attribute. This tool highlights
    selected points by fading them out (i.e., reducing their alpha values).
    If your collection has alpha < 1, this tool will permanently alter them.

    Note that this tool selects collection objects based on their *origins*
    (i.e., `offsets`).

    Parameters
    ----------
    ax : :class:`~matplotlib.axes.Axes`
        Axes to interact with.

    collection : :class:`matplotlib.collections.Collection` subclass
        Collection you want to select from.

    alpha_other : 0 <= float <= 1
        To highlight a selection, this tool sets all selected points to an
        alpha value of 1 and non-selected points to `alpha_other`.
    
    Example:        
        marker = 'CD45'
        x = adata.obs['X_position']
        y = adata.obs['Y_position']
        
        m_idx = adata.var.index.tolist().index(marker) # Get the index of marker of interest
        tmp_dataframe = pd.DataFrame(adata.X)
        hue = np.array(tmp_dataframe[m_idx])
        
        # Plotting    
        fig, ax = plt.subplots()
        pts = ax.scatter(x, y, s=1,c=hue,cmap='viridis')
        ax.invert_yaxis()
        
        # Function call to do the Lasso selection
        selector = SelectFromCollection(ax, pts)
        # Return indeces of the selected points
        tumor_idx= selector.ind
        len(tumor_idx)
        
        # Update adata
        adata.obs.loc[adata.obs.index[tumor_idx], 'roi-1'] = "roi-1"
        adata.obs['roi-1'].value_counts() # Checking
        
        """

    def __init__(self, ax, collection,alpha_other=0.3):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other

        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)

        # Ensure that we have separate colors for each object
        self.fc = collection.get_facecolors()
        if len(self.fc) == 0:
            raise ValueError('Collection must have a facecolor')
        elif len(self.fc) == 1:
            self.fc = np.tile(self.fc, self.Npts).reshape(self.Npts, -1)

        self.lasso = LassoSelector(ax, onselect=self.onselect)
        self.ind = []

    def onselect(self, verts):
        path = Path(verts)
        self.ind = np.nonzero([path.contains_point(xy) for xy in self.xys])[0]
        self.fc[:, -1] = self.alpha_other
        self.fc[self.ind, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

    def disconnect(self):
        self.lasso.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()