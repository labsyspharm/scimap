# -*- coding: utf-8 -*-
# Created on Fri Mar  6 12:13:22 2020
# @author: Ajit Johnson Nirmal

"""
!!! abstract "Short Description"
    `sm.pp.rescale`: The function allows users to rescale the data. This step is often performed to standardize the 
    the expression of all markers to a common scale. The rescaling can be either performed automatically or manually. 
    User defined gates can be passed to rescale the data manually, else the algorithm fits a GMM (gaussian mixed model) to 
    identify the cutoff point. The resultant data is between 0-1 where values below 0.5 are considered non-expression while 
    above 0.5 is considered positive. 

## Function
"""

# Import library
import os
import pandas as pd
import numpy as np
import itertools
from sklearn.preprocessing import MinMaxScaler
from sklearn.mixture import GaussianMixture
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import seaborn as sns; sns.set(color_codes=True)
from matplotlib.lines import Line2D


def rescale (adata, gate=None, return_gates=False, imageid='imageid', failed_markers=None, method='all',save_fig=False):

    """
Parameters:

    adata : AnnData object

    gate : dataframe, optional  
        DataFrame with first column as markers and second column as the gate values in log1p scale.
        Note: If a marker is not included, the function will try to automatically identify a gate
        based on gaussian mixture modeling. If a marker is included in the `gate` dataframe but
        no values are passed, the marker is simply scaled between 0-1 but does not alter the undelying
        distribution.

    return_gates : boolian, optional  
        Internal parameter for checking.

    failed_markers : list, optional  
        list of markers that are not expressed at all in any cell. pass in as ['CD20', 'CD3D'].

    method : string, optional  
        Two avialble option are- 'all' or 'by_image'. In the event that multiple images were loaded in with distinct 'imageid',
        users have the option to scale all data togeather or each image independently. Please be aware of batch effects when
        passing 'all' with multiple images.

    imageid : string, optional  
        Column name of the column containing the image id.

    save_fig : boolian, optional  
        If True, the gates identified by the GMM method will be saved in a subdirectory
        within your working directory.

Returns:
    AnnData
        Modified AnnData Object.

Example:
```python
    manual_gate = pd.DataFrame({'marker': ['CD3D', 'KI67'], 'gate': [7, 8]})
    adata = sm.pp.rescale (adata, gate=manual_gate, failed_markers=['CD20', 'CD21'])
```
    """


    def rescale_independent (adata, gate, return_gates, failed_markers,save_fig):

        print('Scaling Image '+ str(adata.obs[imageid].unique()))

        # Copy of the raw data if it exisits
        if adata.raw is not None:
            adata.X = adata.raw.X

        data = pd.DataFrame(adata.X, columns = adata.var.index, index= adata.obs.index)
        # Merging the manual gates and non-working markers togeather if any
        if gate is not None:
            m_markers = list(gate.iloc[:,0])
            manual_gate_markers = gate
        if failed_markers != None:
            manual_gate_markers = pd.DataFrame(data[failed_markers].quantile(0.9999999))
            manual_gate_markers['markers'] = failed_markers
            # move column to front
            cols = manual_gate_markers.columns.tolist()
            cols.insert(0, cols.pop(cols.index('markers')))
            manual_gate_markers = manual_gate_markers.reindex(columns= cols)
            manual_gate_markers.columns = ['marker', 'gate']
            m_markers = failed_markers
        if gate is not None and failed_markers != None:
            m_markers = list(gate.iloc[:,0]) + list(manual_gate_markers.iloc[:,0])
            gate.columns = ['marker', 'gate']
            manual_gate_markers = pd.concat([gate, manual_gate_markers])
        if gate is None and failed_markers == None:
            m_markers = []

        # Find markers to send to gmm modelling
        if gate is not None or failed_markers is not None:
            gmm_markers = list(np.setdiff1d(data.columns, m_markers))
        else:
            gmm_markers = list(data.columns)

        # If manual gate is not provided scale the data
        if len(gmm_markers) != 0:
            gmm_data = data[gmm_markers]
            # Clip off the 99th percentile
            def clipping (x):
                clip = x.clip(lower =np.percentile(x,0.01), upper=np.percentile(x,99.99)).tolist()
                return clip
            # Run the function
            gmm_data = gmm_data.apply(clipping)

            # Scaling the data old
            sum_data = gmm_data.sum(axis=1) # Calculate total count for each cell
            n_count = gmm_data.div(sum_data, axis=0) # Divide genes by total count for every cell
            med = np.median(list(itertools.chain(*gmm_data.values.tolist()))) # Calculate median count of the entire dataset
            n_count_med = n_count*med # Multiply by scaling fator (median count of entire dataset)
            n_log = np.log1p(n_count_med) # Log transform data
            scaler = MinMaxScaler(feature_range=(0, 1))
            s = scaler.fit_transform(n_log)
            normalised_data = pd.DataFrame(s, columns = gmm_data.columns, index= gmm_data.index)


            # Gaussian fit to identify the gate for each marker and scale based on the gate
            # Empty data frame to hold the results
            all_gmm_data = pd.DataFrame()
            def gmm_gating (data, marker, return_gates, save_fig, gmm_data):
                # Print
                print('Finding the optimal gate for ' + str(marker))
                # Identify the marker to fit the model
                m = data[marker].values
                # Perform GMM
                data_gm = m.reshape(-1, 1)
                #gmm = GaussianMixture(n_components=2, means_init=[[0],[1]],covariance_type='tied')
                gmm = GaussianMixture(n_components=2)
                gmm.fit(data_gm)
                gate = np.mean(gmm.means_)

                # Find the closest value to the gate
                absolute_val_array = np.abs(m - gate)
                smallest_difference_index = absolute_val_array.argmin()
                closest_element = m[smallest_difference_index]

                # If Save Figure is true
                if save_fig == True:
                    m_ndata = gmm_data[marker].values
                    # generate a linear regression for prediction
                    x = data[marker].values.reshape(-1, 1)
                    y = gmm_data[marker].values.reshape(-1, 1)
                    reg = LinearRegression().fit(x, y)
                    #reg.score(x, y)
                    #SVR
                    #reg = SVR(kernel='poly').fit(x, y)
                    #regr.score(x, y)

                    # the three lines
                    g1 = m[np.abs(m - gmm.means_[0]).argmin()]
                    mg = m[np.abs(m - np.mean(gmm.means_)).argmin()]
                    g2 = m[np.abs(m - gmm.means_[1]).argmin()]
                    # predicted gates on log scale
                    g1 = np.log1p(reg.predict(np.array([[g1]])))
                    mg = np.log1p(reg.predict(np.array([[mg]])))
                    g2 = np.log1p(reg.predict(np.array([[g2]])))

                    # saving figure folder
                    if not os.path.exists('auto_gating'):
                        os.makedirs('auto_gating')

                    # generate figure
                    plt.ioff()
                    lines = [Line2D([0], [0], color='b', linestyle='--'),
                             Line2D([1], [1], color='red'),
                             Line2D([0], [0], color='b',linestyle='--')]
                    labels = ['Neg Gaussian', 'Gate', 'Pos Gaussian']
                    sns.set_style("white")
                    fig, ax = plt.subplots( nrows=1, ncols=1 )
                    sns.distplot(np.log1p(m_ndata),color="grey")
                    plt.axvline(g1,ls='--')
                    plt.axvline(g2,ls='--')
                    plt.axvline(mg,color="red")
                    plt.legend(lines, labels)
                    plt.title(marker, fontsize=30)
                    fig.savefig('auto_gating/' + str(marker) + '.png')
                    plt.clf()
                    plt.close('all')
                    plt.ion()

                # rescale the data based on the identified gate
                marker_study = pd.DataFrame(m, index= data.index)
                marker_study = marker_study.sort_values(0)

                # Find the index of the gate
                gate_index = marker_study.index[marker_study[0] == closest_element][0]

                # Split into high and low groups
                high = marker_study.loc[gate_index:,:]
                low = marker_study.loc[:gate_index,:]

                # Prepare for scaling the high and low dataframes
                scaler_high = MinMaxScaler(feature_range=(0.5, 1))
                scaler_low = MinMaxScaler(feature_range=(0, 0.5))

                # Scale it
                h = pd.DataFrame(scaler_high.fit_transform(high), index = high.index)
                l = pd.DataFrame(scaler_low.fit_transform(low), index = low.index)

                # Merge the high and low and resort it
                scaled_data = pd.concat([l,h])
                scaled_data = scaled_data.loc[~scaled_data.index.duplicated(keep='first')]
                scaled_data = scaled_data.reindex(data.index)

                #return scaled_data
                if return_gates == True:
                    return gate
                else:
                    return scaled_data

            # Apply the function
            r_gmm_gating = lambda x: gmm_gating(data=normalised_data, marker=x,return_gates=return_gates,save_fig=save_fig,gmm_data=gmm_data) # Create lamda function
            all_gmm_data = list(map(r_gmm_gating, gmm_markers)) # Apply function
            all_gmm_data = pd.concat(all_gmm_data, axis=1, sort=False)
            all_gmm_data.columns = gmm_markers
        else:
            all_gmm_data = pd.DataFrame()

        # Empty data frame to hold the results
        all_manual_data = pd.DataFrame()
        if len(m_markers) != 0:
            m_data = np.log1p(data[m_markers])
            # Clip the data
            def clipping (x):
                clip = x.clip(lower =np.percentile(x,1), upper=np.percentile(x,99)).tolist()
                return clip
            # Run the function
            #m_data = m_data.apply(clipping)

            def manual_gating (data,marker,gate):

                # Work on processing manual gates
                m = gate[gate.iloc[:,0] == marker].iloc[:,1].values[0] # gate of the marker passed in

                if np.isnan(m):
                    # Find the mean value of the marker so that it is scaled right at the middle
                    # in other it retains the original scale
                    m = np.mean(data[marker].values)
                    print('Warning: No manual gate was found for ' + str(marker) + '. Scaling it between 0 and 1')
                else:
                    print('Scaling ' + str(marker))


                # Find the closest value to the gate
                absolute_val_array = np.abs(data[marker].values - float(m))
                
                # throw error if the array has nan values
                if np.isnan(absolute_val_array).any():
                    raise ValueError ("An exception occurred: " + str(marker) + ' has nan values')
                
                #absolute_val_array = data[marker].values - float(m)
                # if the gate is above the largest value (essentially the marker has failed)
                #if absolute_val_array.min() < 0:
                #    smallest_difference_index = absolute_val_array.argmax()
                #else:
                #    smallest_difference_index = absolute_val_array.argmin()
                smallest_difference_index = absolute_val_array.argmin()
                closest_element = data[marker].values[smallest_difference_index]

                # rescale the data based on the identified gate
                marker_study = data[marker]
                marker_study = marker_study.sort_values(0)

                # Find the index of the gate
                gate_index = marker_study.index[marker_study == closest_element][0]

                # Split into high and low groups
                high = marker_study[gate_index:]
                low = marker_study[:gate_index]

                # Prepare for scaling the high and low dataframes
                scaler_high = MinMaxScaler(feature_range=(0.5, 1))
                scaler_low = MinMaxScaler(feature_range=(0, 0.5))

                # Scale it
                h = pd.DataFrame(scaler_high.fit_transform(high.values.reshape(-1, 1)), index = high.index)
                l = pd.DataFrame(scaler_low.fit_transform(low.values.reshape(-1, 1)), index = low.index)

                # Merge the high and low and resort it
                scaled_data = pd.concat([l,h])
                scaled_data = scaled_data.loc[~scaled_data.index.duplicated(keep='first')]
                scaled_data = scaled_data.reindex(data.index)

                # Return
                return scaled_data

            # Apply the function
            r_manual_gating = lambda x: manual_gating(data=m_data, marker=x, gate=manual_gate_markers) # Create lamda function
            all_manual_data = list(map(r_manual_gating, m_markers)) # Apply function
            all_manual_data = pd.concat(all_manual_data, axis=1, sort=False)
            all_manual_data.columns = m_markers

        else:
            all_manual_data = pd.DataFrame()


        # If both manual and automatic gating was used, combine them into a single result
        if not all_manual_data.empty:
            all_scaled_data = all_manual_data
        if not all_gmm_data.empty:
            all_scaled_data = all_gmm_data
        if not all_manual_data.empty and not all_gmm_data.empty:
            all_scaled_data = all_gmm_data.merge(all_manual_data, how='outer', left_index=True, right_index=True)

        # re index the columns
        all_scaled_data = all_scaled_data.reindex(columns= data.columns)
        return all_scaled_data

    # Apply method of choice
    if method == 'all':
        all_scaled_data = rescale_independent (adata, gate=gate, return_gates=return_gates, failed_markers=failed_markers, save_fig=save_fig)
    if method == 'by_image':
        adata_list = [adata[adata.obs[imageid] == i] for i in adata.obs[imageid].unique()]
        r_rescale_independent = lambda x: rescale_independent(adata=x, gate=gate, return_gates=return_gates, failed_markers=failed_markers,save_fig=save_fig) # Create lamda function
        scaled_data = list(map(r_rescale_independent, adata_list)) # Apply function
        all_scaled_data = pd.concat(scaled_data)

    # Create a copy of the raw data
    if adata.raw is None:
        adata.raw = adata

    # Replace with normalized data
    adata.X = all_scaled_data

    # Return data
    return adata
