# -*- coding: utf-8 -*-
# Created on Wed Apr 19 16:53:02 2023
# @author: Ajit Johnson Nirmal


"""
!!! abstract "Short Description"
    ComBat is a well-established method for correcting batch effects in high-dimensional data, such as single-cell RNA-seq.
    This implementation uses the `combat` function to correct batch effects across multiple slides.  

## Function
"""

# import libs
from combat.pycombat import pycombat
import pandas as pd
import anndata as ad
import numpy as np
import argparse


def combat(
    adata,
    batch='imageid',
    layers=None,
    log=False,
    replaceOriginal=False,
    label='combat',
):
    """
    Parameters:

        adata (AnnData object):
            Annotated data matrix.

        batch (str, optional):
            The batch key or column in `adata.obs` that indicates the batches for each cell.

        layers (str or None, optional):
            The layer in `adata.layers` that contains the expression data to correct.
            If None, `adata.X` is used. use `raw` to use the data stored in `adata.raw.X`

        log (bool, optional):
            Whether to log transform the data before applying ComBat. Generally use it with `raw`.

        replaceOriginal (bool, optional):
            Whether to replace the original expression data in `adata` with the corrected data.

        label (str, optional):
            The prefix for the key in `adata` that will contain the corrected data.
            If `replaceOriginal` is True, this parameter has no effect.

        Returns:

        adata (anndata):
            The corrected expression data is stored in a new layer `adata.layers['combat']`.

        Examples:

            ```python

            # applying batch correction using raw data
            adata = combat (adata,
                        batch='imageid',
                        layers='raw',
                        log=True,
                        replaceOriginal=False,
                        label='combat')

            # results will be available in adata.layers['combat']

            ```
    """

    # isolate the data
    if layers is None:
        data = pd.DataFrame(adata.X, index=adata.obs.index, columns=adata.var.index)
    elif layers == 'raw':
        data = pd.DataFrame(adata.raw.X, index=adata.obs.index, columns=adata.var.index)
    else:
        data = pd.DataFrame(
            adata.layers[layers], index=adata.obs.index, columns=adata.var.index
        )

    # log the data if requested
    if log is True:
        data = np.log1p(data)

    # isolate batch
    batchData = adata.obs[batch]

    # convert to category
    batchData = batchData.astype('category')

    # make sure there are atleast two batches
    if len(batchData.unique()) < 2:
        raise Exception(
            "Sorry a minimum of 2 batches is required. Please check the '"
            + str(batch)
            + "' column"
        )

    # perform batch correction
    batchCorrected = pycombat(data.T, batchData).T

    # add as a specific layer
    adata.layers[label] = batchCorrected

    # replace original
    if replaceOriginal is True:
        if layers is None:
            adata.X = batchCorrected
        elif layers == 'raw':
            del adata.raw
            adata.raw = ad.AnnData(batchCorrected, obs=adata.obs)
        else:
            adata.layers[layers] = batchCorrected

    # return adata
    return adata


# Make the Function CLI compatable
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run ComBat batch correction.')
    parser.add_argument('--adata', type=str, help='Annotated data matrix.')
    parser.add_argument(
        '--batch',
        type=str,
        default='imageid',
        help='The batch key or column in `adata.obs` that indicates the batches for each cell.',
    )
    parser.add_argument(
        '--layers',
        type=str,
        default=None,
        help='The layer in `adata.layers` that contains the expression data to correct. If None, `adata.X` is used. use `raw` to use the data stored in `adata.raw.X`',
    )
    parser.add_argument(
        '--log',
        type=bool,
        default=False,
        help='Whether to log transform the data before applying ComBat. Generally use it with `raw`.',
    )
    parser.add_argument(
        '--replaceOriginal',
        type=bool,
        default=False,
        help='Whether to replace the original expression data in `adata` with the corrected data.',
    )
    parser.add_argument(
        '--label',
        type=str,
        default='combat',
        help='The prefix for the key in `adata` that will contain the corrected data. If `replaceOriginal` is True, this parameter has no effect.',
    )
    args = parser.parse_args()
    combat(
        adata=args.adata,
        batch=args.batch,
        layers=args.layers,
        log=args.log,
        replaceOriginal=args.replaceOriginal,
        label=args.label,
    )
