#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Fri Nov 27 09:34:22 2020
# @author: Ajit Johnson Nirmal

"""
!!! abstract "Short Description"
    `sm.pl.stacked_barplot`: This function creates stacked bar plots to visualize the 
    distribution and proportions of categories within a specified categorical column 
    across different groups or samples in an AnnData object. It supports both `matplotlib` 
    for generating static plots and `Plotly` for interactive, browser-based visualizations. 
    The flexibility to choose between plotting libraries caters to diverse analysis needs, 
    from detailed publication-ready figures to dynamic exploration of complex datasets, 
    enhancing the interpretability of spatial and phenotypic compositions.

## Function
"""

# Required library
import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns

sns.set(color_codes=True)
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import os


pio.renderers.default = 'browser'
sns.set(style="white")

plt.rcParams['pdf.fonttype'] = 42


# Function
def stacked_barplot(
    adata,
    x_axis='imageid',
    y_axis='phenotype',
    subset_xaxis=None,
    subset_yaxis=None,
    order_xaxis=None,
    order_yaxis=None,
    method='percent',
    plot_tool='matplotlib',
    matplotlib_cmap=None,
    matplotlib_bbox_to_anchor=(1, 1.02),
    matplotlib_legend_loc=2,
    fileName='stacked_barplot.pdf',
    saveDir=None,
    return_data=False,
    **kwargs,
):
    """
    Parameters:
            adata (anndata.AnnData):
                The annotated data matrix.

            x_axis (str):
                Column in `adata.obs` to be used as x-axis categories.

            y_axis (str):
                Column in `adata.obs` representing categories to stack.

            subset_xaxis (list, optional):
                Subsets categories in x_axis before plotting.

            subset_yaxis (list, optional):
                Subsets categories in y_axis before plotting.

            order_xaxis (list, optional):
                Specifies custom ordering for x-axis categories.

            order_yaxis (list, optional):
                Specifies custom ordering for y-axis categories.

            method (str, optional):
                Plotting method; 'percent' for percentage proportions, 'absolute' for actual counts.

            plot_tool (str, optional):
                Choice of plotting library; 'matplotlib' for static plots, 'plotly' for interactive plots.

            matplotlib_cmap (str, optional):
                Matplotlib colormap for coloring the bars.

            matplotlib_bbox_to_anchor (tuple, optional):
                Adjusts the legend's bounding box location in matplotlib plots.

            matplotlib_legend_loc (int, optional):
                Sets the legend location in matplotlib plots.

            return_data (bool, optional):
                If True, returns a DataFrame used for plotting instead of displaying the plot.

            fileName (str, optional):
                Name of the file to save the plot. Relevant only if `saveDir` is not None.

            saveDir (str, optional):
                Directory to save the generated plot. If None, the plot is not saved.

            **kwargs:
                Additional arguments passed to the plotting function (matplotlib or plotly).

    Returns:
        Plot (matplotlib):
            If `return_data` is True, returns a DataFrame containing the data used for plotting.
            Otherwise, displays the stacked bar plot.

    Example:
        ```python

        # Default stacked bar plot showing percentage composition
        sm.pl.stacked_barplot(adata, x_axis='sample_id', y_axis='cell_type', method='percent')

        # Stacked bar plot using absolute counts with matplotlib customization
        sm.pl.stacked_barplot(adata, x_axis='region', y_axis='phenotype', method='absolute', plot_tool='matplotlib',
                        matplotlib_cmap='tab20', figsize=(12, 6), edgecolor='white')

        # Interactive stacked bar plot using Plotly with subset and custom order
        sm.pl.stacked_barplot(adata, x_axis='condition', y_axis='cell_state', subset_xaxis=['Control', 'Treated'],
                        order_yaxis=['State1', 'State2', 'State3'], method='percent', plot_tool='plotly',
                        color_discrete_map={'State1': '#1f77b4', 'State2': '#ff7f0e', 'State3': '#2ca02c'})

        # Retrieve data used for plotting
        data_df = sm.pl.stacked_barplot(adata, x_axis='batch', y_axis='cell_type', return_data=True)

        ```
    """

    # create the dataframe with details
    data = pd.DataFrame(adata.obs)[[x_axis, y_axis]].astype(str)

    # subset the data if needed
    # if subset_data is not None:data = data[data[list(subset_data.keys())[0]].isin(list(subset_data.values())[0])]

    if subset_xaxis is not None:
        if isinstance(subset_xaxis, str):
            subset_xaxis = [subset_xaxis]
        data = data[data[x_axis].isin(subset_xaxis)]
    if subset_yaxis is not None:
        if isinstance(subset_yaxis, str):
            subset_yaxis = [subset_yaxis]
        data = data[data[y_axis].isin(subset_yaxis)]

    # Method: Absolute or Percentile
    if method == 'percent':
        total = data.groupby([x_axis, y_axis]).size().unstack().fillna(0).sum(axis=1)
        rg = pd.DataFrame(
            data.groupby([x_axis, y_axis])
            .size()
            .unstack()
            .fillna(0)
            .div(total, axis=0)
            .stack()
        )
    elif method == 'absolute':
        rg = pd.DataFrame(
            data.groupby([x_axis, y_axis]).size().unstack().fillna(0).stack()
        )
    else:
        raise ValueError('method should be either percent or absolute')

    # change column name
    rg.columns = ['count']

    # Add the index as columns in the data frame
    rg.reset_index(inplace=True)

    # re-order the x oy y axis if requested by user
    if order_xaxis is not None:
        rg[x_axis] = rg[x_axis].astype('category')
        rg[x_axis] = rg[x_axis].cat.reorder_categories(order_xaxis)
        rg = rg.sort_values(x_axis)
    if order_yaxis is not None:
        rg[y_axis] = rg[y_axis].astype('category')
        rg[y_axis] = rg[y_axis].cat.reorder_categories(order_yaxis)
        rg = rg.sort_values(y_axis)
    if order_xaxis and order_yaxis is not None:
        rg = rg.sort_values([x_axis, y_axis])

    pivot_df = rg.pivot(index=x_axis, columns=y_axis, values='count')

    # Plotting tool
    if plot_tool == 'matplotlib':

        if matplotlib_cmap is None:
            if len(rg[y_axis].unique()) <= 9:
                matplotlib_cmap = "Set1"
            elif len(rg[y_axis].unique()) > 9 and len(rg[y_axis].unique()) <= 20:
                matplotlib_cmap = plt.cm.tab20  # tab20
            else:
                matplotlib_cmap = plt.cm.gist_ncar

        # Plotting
        # add width if not passed via parameters
        try:
            width
        except NameError:
            width = 0.9
        # actual plotting
        # p = pivot_df.plot.bar(stacked=True, cmap=matplotlib_cmap, width=width,  **kwargs)
        # handles, labels = p.get_legend_handles_labels() # for reversing the order of the legend
        # p.legend(reversed(handles), reversed(labels), bbox_to_anchor=matplotlib_bbox_to_anchor, loc=matplotlib_legend_loc)

        # Actual plotting
        ax = pivot_df.plot.bar(
            stacked=True, cmap=matplotlib_cmap, width=width, **kwargs
        )
        fig = ax.get_figure()  # Get the Figure object to save
        handles, labels = (
            ax.get_legend_handles_labels()
        )  # for reversing the order of the legend
        ax.legend(
            reversed(handles),
            reversed(labels),
            bbox_to_anchor=matplotlib_bbox_to_anchor,
            loc=matplotlib_legend_loc,
        )

        # Saving the figure if saveDir and fileName are provided
        if saveDir:
            if not os.path.exists(saveDir):
                os.makedirs(saveDir)
            full_path = os.path.join(saveDir, fileName)
            fig.savefig(full_path, dpi=300)  # Use fig.savefig instead of p.savefig
            plt.close(fig)  # Close the figure properly
            print(f"Saved plot to {full_path}")
        else:
            plt.show()

    elif plot_tool == 'plotly':

        fig = px.bar(rg, x=x_axis, y="count", color=y_axis, **kwargs)
        fig.update_layout(
            {'plot_bgcolor': 'rgba(0, 0, 0, 0)', 'paper_bgcolor': 'rgba(0, 0, 0, 0)'},
            xaxis=dict(tickmode='linear'),  # type = 'category'
        )
        fig.show()

    else:

        raise ValueError('plot_tool should be either matplotlib or plotly')

    # Return data
    if return_data is True:
        return pivot_df
