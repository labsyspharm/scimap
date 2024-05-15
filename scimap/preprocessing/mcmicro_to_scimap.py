# -*- coding: utf-8 -*-
# Created on Mon Mar  2 09:12:35 2020
# @author: Ajit Johnson Nirmal

"""
!!! abstract "Short Description"
    `sm.pp.mcmicro_to_scimap`: The function allows users to directly import the output from [mcmicro](https://mcmicro.org/) 
    into `scimap`.

## Function
"""

# Import library
import numpy as np
import anndata as ad
import pandas as pd
import argparse
import sys
import pathlib


def mcmicro_to_scimap(
    feature_table_path,
    remove_dna=True,
    remove_string_from_name=None,
    log=True,
    drop_markers=None,
    random_sample=None,
    unique_CellId=True,
    CellId='CellID',
    split='X_centroid',
    custom_imageid=None,
    min_cells=None,
    verbose=True,
    output_dir=None,
):
    """
    Parameters:
        feature_table_path (list of str):
            A list containing file paths to single-cell spatial feature tables. Each path corresponds to a unique image's feature table.

        remove_dna (bool):
            If set to `True`, channels identified by containing the substring 'dna' are excluded from the final dataset. This parameter is useful for omitting DNA-related features.

        remove_string_from_name (str):
            Specifies a substring to be removed from all channel names, aiding in the normalization of marker names across datasets.

        log (bool):
            If `True`, applies a log transformation (specifically log1p, which adds 1 before taking the logarithm) to the data. This is often used to normalize data distributions.

        drop_markers (list of str):
            A list of marker names to exclude from the analysis. For example, to remove specific markers, you would list them here: ["CD3D", "CD20"].

        random_sample (int):
            Specifies the number of cells to include in a subsampled dataset. This parameter is used for downsampling to a specified cell count.

        CellId (str):
            The name of the column in the input data that contains cell identifiers. This is used to track individual cells across analyses.

        unique_CellId (bool):
            Determines whether to automatically generate a unique identifier for each cell by combining the `CellId` with an `imageid`. Set to `False` to use the original `CellId` values directly. This is crucial for maintaining unique cell identifiers, especially when integrating multiple datasets.

        split (str):
            The name of the column that demarcates the boundary between quantitative marker data and additional metadata in the input CSV. Specifying this column allows for the separation of data into counts and metadata components.

        custom_imageid (str):
            Allows for the specification of a custom Image ID for each dataset. By default, the Image ID is derived from the name of the input CSV file.

        min_cells (int):
            Sets a threshold for the minimum number of cells required for an image to be included in the analysis. Images with fewer cells than this number are excluded. This is useful for filtering out datasets with sparse cellular data.

        verbose (bool):
            If set to `True`, the function will print detailed messages about its progress and the steps being executed.

        output_dir (str):
            The file path to the directory where output files will be saved. This parameter specifies the destination for any generated files.


    Returns:
        AnnData Object (anndata):
            An annotated data object containing the processed single-cell spatial features, ready for downstream analysis.


    Example:
            ```

            feature_table_path = ['/Users/aj/scimapExampleData/quantification/exemplar-001--unmicst_cell.csv',
                                  '/Users/aj/scimapExampleData/quantification/exemplar-002--unmicst_cell.csv']

            adata = sm.pp.mcmicro_to_scimap(feature_table_path, drop_markers=['ELANE'], random_sample=5000)

            ```

    """

    # feature_table_path list or string
    if isinstance(feature_table_path, str):
        feature_table_path = [feature_table_path]
    feature_table_path = [pathlib.Path(p) for p in feature_table_path]

    # Import data based on the location provided
    def load_process_data(image):
        # Print the data that is being processed
        if verbose:
            print(f"Loading {image.name}")
        d = pd.read_csv(image)
        # If the data does not have a unique image ID column, add one.
        if 'imageid' not in d.columns:
            if custom_imageid is not None:
                imid = custom_imageid
            else:
                # imid = random.randint(1000000,9999999)
                imid = image.stem
            d['imageid'] = imid
        # Unique name for the data
        if unique_CellId is True:
            d.index = d['imageid'].astype(str) + '_' + d[CellId].astype(str)
        else:
            d.index = d[CellId]

        # move image id and cellID column to end
        cellid_col = [col for col in d.columns if col != CellId] + [CellId]
        d = d[cellid_col]
        imageid_col = [col for col in d.columns if col != 'imageid'] + ['imageid']
        d = d[imageid_col]
        # If there is INF replace with zero
        d = d.replace([np.inf, -np.inf], 0)
        # Return data
        return d

    # Apply function to all images and create a master dataframe
    r_load_process_data = lambda x: load_process_data(image=x)  # Create lamda function
    all_data = list(
        map(r_load_process_data, list(feature_table_path))
    )  # Apply function

    # Merge all the data into a single large dataframe
    for i in range(len(all_data)):
        all_data[i].columns = all_data[0].columns
    entire_data = pd.concat(all_data, axis=0, sort=False)

    # Randomly sample the data
    if random_sample is not None:
        entire_data = entire_data.sample(n=random_sample, replace=False)

    # Remove the images that contain less than a defined threshold of cells (min_cells)
    if min_cells is not None:
        to_drop = (
            entire_data['imageid']
            .value_counts()[entire_data['imageid'].value_counts() < min_cells]
            .index
        )
        entire_data = entire_data[~entire_data['imageid'].isin(to_drop)]
        print(
            'Removed Images that contained less than '
            + str(min_cells)
            + ' cells: '
            + str(to_drop.values)
        )

    # Split the data into expression data and meta data
    # Step-1 (Find the index of the column with name Area)
    split_idx = entire_data.columns.get_loc(split)
    meta = entire_data.iloc[:, split_idx:]
    # Step-2 (select only the expression values)
    entire_data = entire_data.iloc[:, :split_idx]

    # Rename the columns of the data
    if remove_string_from_name is not None:
        entire_data.columns = entire_data.columns.str.replace(
            remove_string_from_name, ''
        )

    # Save a copy of the column names in the uns space of ANNDATA
    markers = list(entire_data.columns)

    # Remove DNA channels
    if remove_dna is True:
        entire_data = entire_data.loc[
            :, ~entire_data.columns.str.contains('dna', case=False)
        ]

    # Drop unnecessary markers
    if drop_markers is not None:
        if isinstance(drop_markers, str):
            drop_markers = [drop_markers]
        entire_data = entire_data.drop(columns=drop_markers)

    # Create an anndata object
    adata = ad.AnnData(entire_data)
    adata.obs = meta
    adata.uns['all_markers'] = markers

    # Add log data
    if log is True:
        adata.raw = adata
        adata.X = np.log1p(adata.X)
        adata.layers['log'] = np.log1p(adata.X)

    # Save data if requested
    if output_dir is not None:
        output_dir = pathlib.Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        imid = feature_table_path[0].stem
        adata.write(output_dir / f'{imid}.h5ad')
        # adata.write(str(output_dir) + '/' + imid + '.h5ad')
    else:
        # Return data
        return adata


def main(argv=sys.argv):
    parser = argparse.ArgumentParser(
        description='The function allows users to directly import the output from mcmicro'
    )
    parser.add_argument(
        '--feature_table_path',
        nargs='*',
        required=True,
        help='List of path to the single-cell feature tables. Each Image should have a unique path supplied.',
    )
    parser.add_argument(
        '--remove_dna',
        action='store_true',
        required=False,
        default=True,
        help='Remove the DNA channels from the final output. Looks for channels with the string dna in it.',
    )
    parser.add_argument(
        '--remove_string_from_name',
        type=str,
        required=False,
        default=None,
        help='Used to celan up channel names. If a string is given, that particular string will be removed from all marker names. If multiple images are passed, just use the string that appears in the first image.',
    )
    parser.add_argument(
        '--log',
        required=False,
        default=True,
        help='Log the data (log1p transformation will be applied).',
    )
    parser.add_argument(
        '--drop_markers',
        nargs='*',
        required=False,
        default=None,
        help='List of markers to drop from the analysis. e.g. ["CD3D", "CD20"]',
    )
    parser.add_argument(
        '--random_sample',
        type=int,
        required=False,
        default=None,
        help='Randomly sub-sample the data with the desired number of cells.',
    )
    parser.add_argument(
        '--unique_CellId',
        required=False,
        default=True,
        help='By default, the function creates a unique name for each cell/row by combining the `CellId` and `imageid`. If you wish not to perform this operation please pass `False`. The function will use whatever is under `CellId`. In which case, please be careful to pass unique `CellId` especially when loading multiple datasets togeather.',
    )
    parser.add_argument(
        '--CellId',
        type=str,
        required=False,
        default='CellID',
        help='Name of the column that contains the cell ID.',
    )
    parser.add_argument(
        '--split',
        type=str,
        required=False,
        default='X_centroid',
        help='To split the CSV into counts table and meta data, pass in the name of the column that immediately follows the marker quantification.',
    )
    parser.add_argument(
        '--custom_imageid',
        type=str,
        required=False,
        default=None,
        help='Pass a user defined Image ID. By default the name of the CSV file is used.',
    )
    parser.add_argument(
        '--min_cells',
        type=int,
        required=False,
        default=None,
        help='If these many cells are not in the image, the image will be dropped. Particulary useful when importing multiple images.',
    )
    parser.add_argument(
        '--verbose',
        required=False,
        default=True,
        help='The function will print detailed messages about its progress.',
    )
    parser.add_argument(
        '--output_dir',
        type=str,
        required=False,
        default=None,
        help='Path to output directory.',
    )
    args = parser.parse_args(argv[1:])
    print(vars(args))
    mcmicro_to_scimap(**vars(args))


if __name__ == '__main__':
    main()
