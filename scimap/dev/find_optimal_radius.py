import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.neighbors import BallTree
from anndata import AnnData, read_h5ad
from typing import Union, Optional, List

def plot_neighbor_distribution(
    adata: Union[str, AnnData],
    x_coordinate: str = 'X_centroid',
    y_coordinate: str = 'Y_centroid',
    z_coordinate: Optional[str] = None,
    imageid: str = 'imageid',
    subset: Optional[str] = None,
    radius_max: float = 100.0,
    radius_step: float = 10.0,
    return_counts: bool = False
) -> Optional[pd.DataFrame]:
    """
    Plot distribution of neighbor counts within increasing radii using BallTree,
    computed per image and optionally for a subset image.
    
    Parameters
    ----------
    adata : Union[str, AnnData]
        AnnData object or path to a .h5ad file.
    x_coordinate, y_coordinate, z_coordinate : str
        Coordinate column names in adata.obs.
    imageid : str
        Column in adata.obs that defines unique image IDs.
    subset : Optional[str]
        If set, only process this imageid.
    radius_max : float
        Maximum radius to consider for neighbor counting.
    radius_step : float
        Incremental step for radius.
    return_counts : bool
        If True, returns the neighbor count DataFrame.
        
    Returns
    -------
    Optional[pd.DataFrame]
        DataFrame of neighbor counts if return_counts=True.
    """

    # Load from file if path
    if isinstance(adata, str):
        adata = read_h5ad(adata)

    # Determine list of image-level AnnData objects
    if subset is not None:
        adata_list = [adata[adata.obs[imageid] == subset].copy()]
    else:
        adata_list = [adata[adata.obs[imageid] == img].copy() for img in adata.obs[imageid].unique()]

    def _neighbor_distribution_per_image(
        ad: AnnData,
        radii: np.ndarray
    ) -> pd.DataFrame:
        obs = ad.obs
    
        if z_coordinate is not None:
            coord_cols = [x_coordinate, y_coordinate, z_coordinate]
        else:
            coord_cols = [x_coordinate, y_coordinate]
    
        # Drop missing coordinate rows
        valid_obs = obs.dropna(subset=coord_cols)
        coords = valid_obs[coord_cols].to_numpy()
    
        if len(coords) == 0:
            raise ValueError("No valid coordinates found after dropping NA values.")
    
        tree = BallTree(coords, leaf_size=40)
        counts = []
    
        for r in radii:
            ind = tree.query_radius(coords, r=r)
            count = np.array([len(i) - 1 for i in ind])  # exclude self
            counts.append(count)
    
        count_df = pd.DataFrame(
            np.vstack(counts).T,
            index=valid_obs.index,
            columns=[f"{r:.1f}" for r in radii]
        )
        count_df[imageid] = valid_obs[imageid].values
        return count_df


    # Compute neighbor counts across all images
    radii = np.arange(radius_step, radius_max + radius_step, radius_step)
    all_counts: List[pd.DataFrame] = [
        _neighbor_distribution_per_image(ad, radii)
        for ad in adata_list
    ]

    final_df = pd.concat(all_counts)

    # Prepare for plotting
    plot_df = final_df.reset_index().melt(id_vars=['index', imageid], var_name='Radius', value_name='Neighbor Count')
    plot_df.rename(columns={'index': 'Cell'}, inplace=True)

    # Plot
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=plot_df, x='Radius', y='Neighbor Count', hue=imageid)
    plt.xticks(rotation=45)
    plt.title('Distribution of Neighbor Counts by Radius (per image)')
    plt.tight_layout()
    plt.legend(title=imageid, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.show()

    if return_counts:
        return final_df
