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
import argparse
from sklearn.preprocessing import MinMaxScaler
from sklearn.mixture import GaussianMixture
import anndata as ad
import json


# Function
def rescale(
    adata,
    gate=None,
    log=True,
    imageid='imageid',
    failed_markers=None,
    method='all',
    verbose=True,
    random_state=0,
    gmm_components=3,
):
    """
    Parameters:
        adata (AnnData Object, required):
            An annotated data object that contains single-cell expression data.

        gate (DataFrame, optional):
            A pandas DataFrame where the first column lists markers, and subsequent columns contain gate values
            for each image in the dataset. Column names must correspond to unique `imageid` identifiers, and the marker column must be named "markers".
            If a single column of gate values is provided for a dataset with multiple images, the same gate will be uniformly applied to all images.
            In this case, ensure that the columns are named exactly "markers" and "gates".
            If no gates are provided for specific markers, the function attempts to automatically determine gates using a Gaussian Mixture Model (GMM).

            Note: If you have used `napariGater()`, the gates are stored within `adata.uns['gates']`.
            You can directly pass `adata.uns['gates']` to use these pre-defined gates.

        log (bool, optional):
            If `True`, the data in `adata.raw.X` will be log-transformed (using log1p) before gate application. This transformation is recommended when automatic gate identification through GMM is performed, as it helps in normalizing data distributions.

        imageid (str, optional):
            The name of the column in `adata` that contains Image IDs. This is necessary for matching manual gates specified in the `gate` DataFrame to their respective images.

        failed_markers (dict, optional):
            A dictionary mapping `imageid` to markers that failed quality control. This allows for the exclusion of specific markers from the analysis based on prior visual inspection or other criteria. The dictionary can use `all` as a key to specify markers that failed across all images.

        method (str, optional):
            Specifies the gating strategy: `all` to pool data from all images for GMM application, or `by_image` to apply GMM separately for each image. `all` may introduce batch effects, while `by_image` requires sufficient variation within each image to distinguish negative from positive populations effectively.

        random_state (int, optional):
            The seed used by the random number generator for GMM. Ensures reproducibility of results.

        verbose (bool, optional):
            If `True`, detailed progress updates and diagnostic messages will be printed during the function's execution.

        gmm_components (int, optional):
            Number of components to use in the Gaussian Mixture Model for automatic gating.
            Must be at least 2. Gate will be placed between the highest two components.
            Default is 3.

    Returns:
        Modified AnnData Object (AnnData):
            Returns the input `adata` object with updated expression data (`adata.X`) after rescaling. The gates applied, either provided manually or determined automatically, are stored within `adata.uns['gates']`.

    Example:
        ```python

        # Example with manual gates
        manual_gate = pd.DataFrame({'marker': ['CD3D', 'KI67'], 'gate': [7, 8]})
        adata = sm.pp.rescale(adata, gate=manual_gate, failed_markers={'all': ['CD20', 'CD21']})

        # Importing gates from a CSV
        manual_gate = pd.read_csv('manual_gates.csv')
        adata = sm.pp.rescale(adata, gate=manual_gate, failed_markers={'all': ['CD20', 'CD21']})

        # Running without manual gates to use GMM for automatic gate determination
        adata = sm.pp.rescale(adata, gate=None, failed_markers={'all': ['CD20', 'CD21']})

        ```

    """

    # log=True; imageid='imageid'; failed_markers=None; method='all'; random_state=0

    # make a copy to raw data if raw is none
    if adata.raw is None:
        adata.raw = adata

    # Mapping between markers and gates in the given dataset
    dataset_markers = adata.var.index.tolist()
    dataset_images = adata.obs[imageid].unique().tolist()
    m = pd.DataFrame(index=dataset_markers, columns=dataset_images).reset_index()
    m = pd.melt(m, id_vars=[m.columns[0]])
    m.columns = ['markers', 'imageid', 'gate']

    # Manipulate m with and without provided manual gates
    if gate is None:
        gate_mapping = m.copy()
    else:
        # Check overlap between gate columns and dataset images
        matching_images = set(gate.columns) & set(dataset_images)

        # link to make sure index name is markers as we use reset_index later
        if gate.index.name != 'markers' and 'markers' not in gate.columns:
            gate.index.name = 'markers'

        if len(matching_images) == 0 and len(gate.columns) > 0:
            # Case 1: No matching images and single value column - apply globally
            gate = gate.reset_index()  # Convert index to column
            gate_mapping = m.copy()
            gate_mapping.gate = gate_mapping.gate.fillna(
                gate_mapping.markers.map(
                    dict(
                        zip(gate['markers'], gate['gates'])
                    )  # these columns are hardcoded in CSV
                )
            )
        else:
            # Case 2: handles both if all imageid matches with gate columns or if they partially match
            gate = gate.reset_index()
            manual_m = pd.melt(gate, id_vars=gate[['markers']])
            manual_m.columns = ['markers', 'imageid', 'm_gate']
            gate_mapping = pd.merge(
                m,
                manual_m,
                how='left',
                left_on=['markers', 'imageid'],
                right_on=['markers', 'imageid'],
            )
            gate_mapping['gate'] = gate_mapping['gate'].fillna(gate_mapping['m_gate'])
            gate_mapping = gate_mapping.drop(columns='m_gate')

    # Addressing failed markers
    def process_failed(adata_subset, foramted_failed_markers):
        if verbose:
            print(
                'Processing Failed Marker in '
                + str(adata_subset.obs[imageid].unique()[0])
            )
        # prepare data
        data_subset = pd.DataFrame(
            adata_subset.raw.X,
            columns=adata_subset.var.index,
            index=adata_subset.obs.index,
        )
        if log is True:
            data_subset = np.log1p(data_subset)

        # subset markers in the subset
        fm_sub = foramted_failed_markers[adata_subset.obs[imageid].unique()].dropna()

        def process_failed_internal(fail_mark, data_subset):
            return data_subset[fail_mark].max()

        r_process_failed_internal = lambda x: process_failed_internal(
            fail_mark=x, data_subset=data_subset
        )
        f_g = list(map(r_process_failed_internal, [x[0] for x in fm_sub.values]))
        subset_gate = pd.DataFrame(
            {
                'markers': [x[0] for x in fm_sub.values],
                'imageid': adata_subset.obs[imageid].unique()[0],
                'gate': f_g,
            }
        )
        # return
        return subset_gate

    # Identify the failed markers
    if failed_markers is not None:
        # check if failed marker is a dict
        if isinstance(failed_markers, dict) is False:
            raise ValueError(
                '`failed_markers` should be a python dictionary, please refer documentation'
            )
        # create a copy
        fm = failed_markers.copy()
        # seperate all from the rest
        if 'all' in failed_markers:
            all_failed = failed_markers['all']
            if isinstance(all_failed, str):
                all_failed = [all_failed]
            failed_markers.pop('all', None)

            df = pd.DataFrame(columns=adata.obs[imageid].unique())
            for i in range(len(all_failed)):
                df.loc[i] = np.repeat(all_failed[i], len(df.columns))
            # for i in  range(len(df.columns)):
            #    df.loc[i] = all_failed[i]
        # rest of the failed markers
        # fail = pd.DataFrame.from_dict(failed_markers)
        fail = pd.DataFrame(
            dict([(k, pd.Series(v)) for k, v in failed_markers.items()])
        )
        # merge
        if 'all' in fm:
            foramted_failed_markers = pd.concat([fail, df], axis=0)
        else:
            foramted_failed_markers = fail

        # send the adata objects that need to be processed
        # Check if any image needs to pass through the GMM protocol
        adata_list = [
            adata[adata.obs[imageid] == i] for i in foramted_failed_markers.columns
        ]
        # apply the process_failed function
        r_process_failed = lambda x: process_failed(
            adata_subset=x, foramted_failed_markers=foramted_failed_markers
        )
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
    def clipping(x):
        clip = x.clip(
            lower=np.percentile(x, 0.01), upper=np.percentile(x, 99.99)
        ).tolist()
        return clip

    # Find GMM based gates
    def gmm_gating(marker, data, gmm_components):
        """Internal function to identify gates using GMM

        Parameters:
            marker: marker name
            data: expression data
            gmm_components: number of components for GMM (minimum 2)
        """
        # Ensure minimum of 2 components
        gmm_components = max(2, gmm_components)

        # Prepare data for GMM
        data_gm = data[marker].values.reshape(-1, 1)
        data_gm = data_gm[~np.isnan(data_gm), None]

        # Fit GMM with gmm_components
        gmm = GaussianMixture(
            n_components=gmm_components, random_state=random_state
        ).fit(data_gm)

        # Sort components by their means
        means = gmm.means_.flatten()
        sorted_idx = np.argsort(means)
        sorted_means = means[sorted_idx]

        # Calculate gate as midpoint between the second-to-last and last components
        gate = np.mean([sorted_means[-2], sorted_means[-1]])

        return gate

    # Running gmm_gating on the dataset
    def gmm_gating_internal(adata_subset, gate_mapping, method):
        if verbose:
            print(
                'Running GMM for image: ' + str(adata_subset.obs[imageid].unique()[0])
            )
        data_subset = pd.DataFrame(
            adata_subset.raw.X,
            columns=adata_subset.var.index,
            index=adata_subset.obs.index,
        )
        # find markers
        if method == 'all':
            image_specific = gate_mapping.copy()
            marker_to_gate = list(
                gate_mapping[gate_mapping.gate.isnull()].markers.unique()
            )
        else:
            image_specific = gate_mapping[
                gate_mapping['imageid'].isin(adata_subset.obs[imageid].unique())
            ]
            marker_to_gate = image_specific[image_specific.gate.isnull()].markers.values

        if verbose and len(marker_to_gate) > 0:
            print('Applying GMM to markers: ' + ', '.join(marker_to_gate))

        # Apply clipping
        data_subset_clipped = data_subset.apply(clipping)
        # log transform data
        if log is True:
            data_subset_clipped = np.log1p(data_subset_clipped)
        # identify the gates for the markers
        r_gmm_gating = lambda x: gmm_gating(
            marker=x, data=data_subset_clipped, gmm_components=gmm_components
        )
        gates = list(map(r_gmm_gating, marker_to_gate))
        # create a df with results
        result = image_specific[image_specific.gate.isnull()]
        mapping = dict(zip(marker_to_gate, gates))
        for i in result.index:
            result.loc[i, 'gate'] = mapping[result.loc[i, 'markers']]
        # result['gate'] = result['gate'].fillna(result['markers'].map(dict(zip(marker_to_gate, gates))))
        # return
        return result

    # Create a list of image IDs that need to go through the GMM
    gmm_images = gate_mapping[gate_mapping.gate.isnull()].imageid.unique()

    # Check if any image needs to pass through the GMM protocol
    if len(gmm_images) > 0:
        # Create a list of adata that need to go through the GMM
        if method == 'all':
            adata_list = [adata]
        else:
            adata_list = [adata[adata.obs[imageid] == i] for i in gmm_images]
        # run function
        r_gmm_gating_internal = lambda x: gmm_gating_internal(
            adata_subset=x, gate_mapping=gate_mapping, method=method
        )
        all_gates = list(map(r_gmm_gating_internal, adata_list))

        # combine the results and merge with gate_mapping
        result = []
        for i in range(len(all_gates)):
            result.append(all_gates[i])
        result = pd.concat(result, join='outer')
        # use this to merge with gate_mapping
        gate_mapping.gate = gate_mapping.gate.fillna(
            gate_mapping.markers.map(dict(zip(result.markers, result.gate)))
        )

    # Rescaling function
    def data_scaler(adata_subset, gate_mapping):
        if verbose:
            print('\nScaling Image: ' + str(adata_subset.obs[imageid].unique()[0]))
        # Organise data
        data_subset = pd.DataFrame(
            adata_subset.raw.X,
            columns=adata_subset.var.index,
            index=adata_subset.obs.index,
        )
        if log is True:
            data_subset = np.log1p(data_subset)
        # subset markers in the subset
        gate_mapping_sub = gate_mapping[
            gate_mapping['imageid'] == adata_subset.obs[imageid].unique()[0]
        ]

        # organise gates
        def data_scaler_internal(marker, gate_mapping_sub):
            if verbose:
                gate_value = gate_mapping_sub[gate_mapping_sub.markers == marker][
                    'gate'
                ].values[0]
                print(f'Scaling {marker} (gate: {gate_value:.3f})')
            # find the gate
            moi = gate_mapping_sub[gate_mapping_sub.markers == marker]['gate'].values[0]

            # Find the closest value to the gate
            absolute_val_array = np.abs(data_subset[marker].values - float(moi))
            # throw error if the array has nan values
            if np.isnan(absolute_val_array).any():
                raise ValueError(
                    "An exception occurred: " + str(marker) + ' has nan values'
                )
            # smallest diff
            smallest_difference_index = absolute_val_array.argmin()
            closest_element = data_subset[marker].values[smallest_difference_index]

            # rescale the data based on the identified gate
            marker_study = data_subset[marker]
            marker_study = marker_study.sort_values(axis=0)
            # Find the index of the gate
            # account for 0
            if all(marker_study == 0):
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
            h = pd.DataFrame(
                scaler_high.fit_transform(high.values.reshape(-1, 1)), index=high.index
            )
            l = pd.DataFrame(
                scaler_low.fit_transform(low.values.reshape(-1, 1)), index=low.index
            )
            # Merge the high and low and resort it
            scaled_data = pd.concat([l, h])
            scaled_data = scaled_data.loc[~scaled_data.index.duplicated(keep='first')]
            scaled_data = scaled_data.reindex(data_subset.index)
            # scaled_data[scaled_data > 0.5].count(axis=1).sum()
            # return
            return scaled_data

        # run internal function
        r_data_scaler_internal = lambda x: data_scaler_internal(
            marker=x, gate_mapping_sub=gate_mapping_sub
        )
        scaled_subset = list(
            map(r_data_scaler_internal, gate_mapping_sub.markers.values)
        )

        # combine the results and merge with gate_mapping
        scaled_subset_result = []
        for i in range(len(scaled_subset)):
            scaled_subset_result.append(scaled_subset[i])
        scaled_subset_result = pd.concat(scaled_subset_result, join='outer', axis=1)
        scaled_subset_result.columns = gate_mapping_sub.markers.values
        # scaled_subset_result[scaled_subset_result['CD3E'] > 0.5]['CD3E'].count(axis=1).sum()

        # return
        return scaled_subset_result

    # pass each dataset seperately
    adata_list = [adata[adata.obs[imageid] == i] for i in adata.obs[imageid].unique()]

    # Run the scaler function
    r_data_scaler = lambda x: data_scaler(adata_subset=x, gate_mapping=gate_mapping)
    scaled_subset = list(map(r_data_scaler, adata_list))

    # combine the results and merge with gate_mapping
    final_result = []
    for i in range(len(scaled_subset)):
        final_result.append(scaled_subset[i])
    final_result = pd.concat(final_result, join='outer')

    # reindex the final_results
    final_result = final_result.reindex(adata.obs.index)

    # save final gates
    adata.uns['gates'] = gate_mapping.pivot_table(
        index=['markers'], columns=['imageid']
    ).droplevel(
        0, axis=1
    )  # .reset_index()

    # add to the anndata
    adata.X = final_result

    # return adata
    return adata


# Make the Function CLI compatible
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='The function allows users to rescale the data.'
    )

    parser.add_argument('--adata', type=str, help='Path to AnnData object file (.h5ad)')

    parser.add_argument(
        '--gate',
        type=str,
        help='Path to gates CSV file. First column should be markers, subsequent columns should be gate values',
    )

    parser.add_argument(
        '--log',
        type=bool,
        default=True,
        help='If True, log transform (log1p) before applying gates',
    )

    parser.add_argument(
        '--imageid',
        type=str,
        default='imageid',
        help='Column name containing Image IDs',
    )

    parser.add_argument(
        '--failed_markers',
        type=str,
        help='Path to JSON file containing failed markers dictionary',
    )

    parser.add_argument(
        '--method',
        type=str,
        default='all',
        choices=['all', 'by_image'],
        help='Method for GMM: all or by_image',
    )

    parser.add_argument(
        '--verbose', type=bool, default=True, help='Print detailed progress messages'
    )

    parser.add_argument(
        '--random_state', type=int, default=0, help='Random seed for GMM'
    )

    parser.add_argument(
        '--gmm_components',
        type=int,
        default=3,
        help='Number of components for GMM (minimum 2)',
    )

    args = parser.parse_args()

    # Load the AnnData object
    adata = ad.read_h5ad(args.adata)

    # Load gates if provided
    gate = None
    if args.gate:
        gate = pd.read_csv(args.gate)

    # Load failed markers if provided
    failed_markers = None
    if args.failed_markers:
        with open(args.failed_markers, 'r') as f:
            failed_markers = json.load(f)

    # Call the function
    adata = rescale(
        adata=adata,
        gate=gate,
        log=args.log,
        imageid=args.imageid,
        failed_markers=failed_markers,
        method=args.method,
        verbose=args.verbose,
        random_state=args.random_state,
        gmm_components=args.gmm_components,
    )

    # Save the modified AnnData object back to disk
    adata.write_h5ad(args.adata)
