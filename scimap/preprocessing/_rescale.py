# -*- coding: utf-8 -*-
# Created on Fri Mar  6 12:13:22 2020
# @author: Ajit Johnson Nirmal

"""
!!! abstract "Short Description"
    `sm.pp.rescale`: The function allows users to rescale the data. This step is often performed to standardize the 
    the expression of all markers to a common scale. The rescaling can be either performed automatically or manually. 
    User defined gates can be passed to rescale the data manually, else the algorithm fits a GMM (gaussian mixed model) to 
    identify the cutoff point. The resultant data is between 0-1 where values below 0.5 are considered non-expressing while 
    above 0.5 is considered positive. 

## Function
"""

# Import library
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.mixture import GaussianMixture


# Function
def rescale (adata, gate=None, log=True,
             imageid='imageid', failed_markers=None,
              method='all',random_state=0):
    """
Parameters:

    adata : AnnData Object  

    gate : dataframe, optional  
        DataFrame with first column as markers and subsequent column with gate values for each image in the dataset.
        The column names should correspond to the unique `imageid`. If only one column of gate is provied 
        to a dataset with multiple images, the same gate will be applied to all images.
        Note: If gates are not provided or left out for a particular marker, the function will try to 
        automatically identify a gate based on applying gaussian mixture modeling algorithm (GMM). The default is None.
        
    log : bool, optional  
        By default the data stored in `adata.raw.X` is extracted for scaling. If the user wishes to log transform (log1p)
        it before applying the gates, this parameter can be set to True. Please note if the function is used to 
        identify gates based on GMM, it is recommended for the data to be log transformed. The default is True.
        
    imageid : string, optional  
        The column containing the Image IDs. When passing manual gates the columns of the dataframe need to match 
        to the elements within the passed `imageid` column. The default is 'imageid'.
        
    failed_markers : dict, optional  
        Markers that were deemed to have failed based on prior visual inspection. This parameter accepts a python 
        dictionary with `key` as `imageid` and `value` as markers that failed in that particular `imageid`. 
        Example: `failed_markers = {'image_1': ['failed_marker_1'], 'image_2' : ['failed_marker_1', 'failed_marker_2']}`. 
        To make it easier to allow specifying markers that failed in `all` images within the dataset, the parameter also 
        recognizes the special keyword `all`. For example, `failed_markers = {'all': ['failed_marker_X'], 'image_2' : ['failed_marker_1', 'failed_marker_2']}`. 
        The default is None.
        
    method : string, optional  
        Two avialble option are- 'all' or 'by_image'. In the event that multiple images were loaded in with distinct 'imageid',
        users have the option to apply GMM by pooling all data togeather or to apply it to each image independently. 
        Please be aware of batch effects when passing 'all' to multiple images. In contrast, if there are not enough variation 
        within individual images, the GMM cannot reliably distinguish between the negative and positive populations as well.  
        
    random_state : int, optional  
        Seed for GMM. The default is 0.

Returns:

    Modified AnnData Object
        The values in `adata.X` are replaced with the scaled data.
        The final gates used for saving the data is also stored in `adata.uns['gates']`

Example:
```python
# create a df with manual gates
manual_gate = pd.DataFrame({'marker': ['CD3D', 'KI67'], 'gate': [7, 8]}) 
adata = sm.pp.rescale (adata, gate=manual_gate, failed_markers={'all':['CD20', 'CD21']})
    
# you could also import the gates as a pandas dataframe without index
manual_gate = pd.read_csv('manual_gates.csv')
adata = sm.pp.rescale (adata, gate=manual_gate, failed_markers={'all':['CD20', 'CD21']})
    
# The function can also be run without providing manual gates. This will trigger the GMM mode
adata = sm.pp.rescale (adata, gate=None, failed_markers={'all':['CD20', 'CD21']})
    
```

    """
    
    # make a copy to raw data if raw is none
    if adata.raw is None:
        adata.raw = adata
    
    # Mapping between markers and gates in the given dataset
    dataset_markers = list(adata.var.index)
    dataset_images = list(adata.obs[imageid].unique())
    m= pd.DataFrame(index=dataset_markers, columns=dataset_images).reset_index()
    m= pd.melt(m, id_vars=[m.columns[0]])
    m.columns = ['markers', 'imageid', 'gate']
    # Manipulate m with and without provided manual fates
    if gate is None:
        gate_mapping = m.copy()
    elif bool(set(list(gate.columns)) & set(dataset_images)) is False:
        global_manual_m = pd.melt(gate, id_vars=[gate.columns[0]])
        global_manual_m.columns = ['markers', 'imageid', 'm_gate']
        gate_mapping = m.copy()
        gate_mapping.gate = gate_mapping.gate.fillna(gate_mapping.markers.map(dict(zip(global_manual_m.markers, global_manual_m.m_gate))))
    else:
        manual_m = pd.melt(gate, id_vars=[gate.columns[0]])
        manual_m.columns = ['markers', 'imageid', 'm_gate']
        gate_mapping = pd.merge(m, manual_m,  how='left', left_on=['markers','imageid'], right_on = ['markers','imageid'])
        gate_mapping['gate'] = gate_mapping['gate'].fillna(gate_mapping['m_gate'])
        gate_mapping = gate_mapping.drop(columns='m_gate')
    
    # Addressing failed markers
    def process_failed (adata_subset, foramted_failed_markers):
        print('Processing Failed Marker in ' + str(adata_subset.obs[imageid].unique()[0]))
        # prepare data
        data_subset = pd.DataFrame(adata_subset.raw.X, columns=adata_subset.var.index, index=adata_subset.obs.index)
        if log is True:
            data_subset = np.log1p(data_subset)
        
        # subset markers in the subset
        fm_sub = foramted_failed_markers[adata_subset.obs[imageid].unique()].dropna()

            
        def process_failed_internal (fail_mark, data_subset):
            return data_subset[fail_mark].max()
        r_process_failed_internal = lambda x: process_failed_internal (fail_mark=x,data_subset=data_subset)
        f_g = list(map(r_process_failed_internal, [ x[0] for x in fm_sub.values]))
        subset_gate = pd.DataFrame( {'markers': [ x[0] for x in fm_sub.values],  
                       'imageid': adata_subset.obs[imageid].unique()[0],
                       'gate': f_g,})     
        # return
        return subset_gate
    
    # Identify the failed markers
    if failed_markers is not None:
        # check if failed marker is a dict
        if isinstance(failed_markers, dict) is False:
            raise ValueError ('`failed_markers` should be a python dictionary, please refer documentation')
        # create a copy 
        fm = failed_markers.copy()
        # seperate all from the rest
        if 'all' in failed_markers:
            all_failed = failed_markers['all']
            if isinstance(all_failed, str):
                all_failed = [all_failed]
            failed_markers.pop('all', None)
            
            df = pd.DataFrame(columns = adata.obs[imageid].unique())
            for i in range(len(all_failed)):
                df.loc[i] = np.repeat(all_failed[i], len(df.columns))
            #for i in  range(len(df.columns)):
            #    df.loc[i] = all_failed[i]
        # rest of the failed markers
        #fail = pd.DataFrame.from_dict(failed_markers)        
        fail = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in failed_markers.items() ]))
        # merge
        if 'all' in fm:
            foramted_failed_markers = pd.concat([fail, df], axis=0)
        else: 
            foramted_failed_markers = fail
        
        # send the adata objects that need to be processed
        # Check if any image needs to pass through the GMM protocol
        adata_list = [adata[adata.obs[imageid] == i] for i in foramted_failed_markers.columns]
        # apply the process_failed function
        r_process_failed = lambda x: process_failed (adata_subset=x,foramted_failed_markers=foramted_failed_markers)
        failed_gates = list(map(r_process_failed, adata_list))    
        # combine the results and merge with gate_mapping
        result = []
        for i in range(len(failed_gates)):
            result.append(failed_gates[i])
        result = pd.concat(result, join='outer')
        # use this to merge with gate_mapping
        x1 = gate_mapping.set_index(['markers', 'imageid'])['gate']
        x2 = result.set_index(['markers', 'imageid'])['gate']
        x1.update(x2)
        gate_mapping = x1.reset_index()
        
    # trim the data before applying GMM
    def clipping (x):
        clip = x.clip(lower =np.percentile(x,0.01), upper=np.percentile(x,99.99)).tolist()
        return clip
            
    # Find GMM based gates
    def gmm_gating (marker, data):
        print('Finding the optimal gate by GMM for ' + str(marker))
        data_gm = data[marker].values.reshape(-1, 1)
        gmm = GaussianMixture(n_components=2, random_state=random_state).fit(data_gm)
        gate = np.mean(gmm.means_)
        return gate
    
    # Running gmm_gating on the dataset
    def gmm_gating_internal (adata_subset, gate_mapping, method):
        print('GMM for ' + str(adata_subset.obs[imageid].unique()))
        data_subset = pd.DataFrame(adata_subset.raw.X, columns=adata_subset.var.index, index=adata_subset.obs.index)      
        # find markers
        if method == 'all':
            image_specific = gate_mapping.copy()
            marker_to_gate = list(gate_mapping[gate_mapping.gate.isnull()].markers.unique())
        else:        
            image_specific = gate_mapping[gate_mapping['imageid'].isin(adata_subset.obs[imageid].unique())]
            marker_to_gate = image_specific[image_specific.gate.isnull()].markers.values   
        # Apply clipping
        data_subset_clipped = data_subset.apply(clipping)
        # log transform data
        if log is True:
            data_subset_clipped = np.log1p(data_subset_clipped)
        # identify the gates for the markers
        r_gmm_gating = lambda x: gmm_gating(marker=x, data=data_subset_clipped) 
        gates = list(map(r_gmm_gating, marker_to_gate))     
        # create a df with results
        result = image_specific[image_specific.gate.isnull()]
        mapping = dict(zip(marker_to_gate, gates))
        for i in result.index:
            result.loc[i, 'gate'] = mapping[result.loc[i, 'markers']]
        #result['gate'] = result['gate'].fillna(result['markers'].map(dict(zip(marker_to_gate, gates))))        
        # return
        return result
    
    
    # Create a list of image IDs that need to go through the GMM
    gmm_images = gate_mapping[gate_mapping.gate.isnull()].imageid.unique()  
    
    # Check if any image needs to pass through the GMM protocol
    if len(gmm_images) > 0 :
        # Create a list of adata that need to go through the GMM
        if method == 'all':
            adata_list = [adata]
        else:
            adata_list = [adata[adata.obs[imageid] == i] for i in gmm_images]
        # run function
        r_gmm_gating_internal = lambda x: gmm_gating_internal (adata_subset=x, 
                                                               gate_mapping=gate_mapping,
                                                               method=method) 
        all_gates = list(map(r_gmm_gating_internal, adata_list))
        
        # combine the results and merge with gate_mapping
        result = []
        for i in range(len(all_gates)):
            result.append(all_gates[i])
        result = pd.concat(result, join='outer')
        # use this to merge with gate_mapping
        gate_mapping.gate = gate_mapping.gate.fillna(gate_mapping.markers.map(dict(zip(result.markers, result.gate))))
            
    
    # Rescaling function
    def data_scaler (adata_subset, gate_mapping):
        print('Scaling Image ' + str(adata_subset.obs[imageid].unique()[0]))
        # Organise data
        data_subset = pd.DataFrame(adata_subset.raw.X, columns=adata_subset.var.index, index=adata_subset.obs.index)
        if log is True:
            data_subset = np.log1p(data_subset)
        # subset markers in the subset
        gate_mapping_sub = gate_mapping[gate_mapping['imageid'] == adata_subset.obs[imageid].unique()[0]]
        
        # organise gates
        def data_scaler_internal (marker, gate_mapping_sub):
            print('Scaling ' + str(marker))
            # find the gate
            moi = gate_mapping_sub[gate_mapping_sub.markers == marker]['gate'].values[0]
            
            # Find the closest value to the gate
            absolute_val_array = np.abs(data_subset[marker].values - float(moi))
            # throw error if the array has nan values
            if np.isnan(absolute_val_array).any():
                raise ValueError ("An exception occurred: " + str(marker) + ' has nan values')
            # smallest diff
            smallest_difference_index = absolute_val_array.argmin()
            closest_element = data_subset[marker].values[smallest_difference_index]
            
            # rescale the data based on the identified gate
            marker_study = data_subset[marker]
            marker_study = marker_study.sort_values(axis=0)
            # Find the index of the gate
            # account for 0
            if all(marker_study) == 0:
                gate_index = pd.DataFrame(marker_study).tail(2).index[0]
            else:
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
            scaled_data = scaled_data.reindex(data_subset.index)
            # return
            return scaled_data
        
        # run internal function
        r_data_scaler_internal = lambda x: data_scaler_internal (marker=x, gate_mapping_sub=gate_mapping_sub) 
        scaled_subset = list(map(r_data_scaler_internal, gate_mapping_sub.markers.values))
            
        # combine the results and merge with gate_mapping
        scaled_subset_result = []
        for i in range(len(scaled_subset)):
            scaled_subset_result.append(scaled_subset[i])
        scaled_subset_result = pd.concat(scaled_subset_result, join='outer', axis=1)
        scaled_subset_result.columns = gate_mapping_sub.markers.values
        
        # return
        return scaled_subset_result
    
    # pass each dataset seperately
    adata_list = [adata[adata.obs[imageid] == i] for i in adata.obs[imageid].unique()]
    
    # Run the scaler function
    r_data_scaler = lambda x: data_scaler (adata_subset=x, gate_mapping=gate_mapping) 
    scaled_subset = list(map(r_data_scaler, adata_list))  
            
    # combine the results and merge with gate_mapping
    final_result = []
    for i in range(len(scaled_subset)):
        final_result.append(scaled_subset[i])
    final_result = pd.concat(final_result, join='outer')
    
    # reindex the final_results
    final_result = final_result.reindex(adata.obs.index)
    
    # save final gates
    adata.uns['gates'] = gate_mapping.pivot_table(index=['markers'], columns=['imageid']).droplevel(0, axis=1)#.reset_index()
    
    # add to the anndata
    adata.X = final_result
    
    # return adata
    return adata
