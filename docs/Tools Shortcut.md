# Tools Shortcut

```python
import scimap as sm
```

### Preprocessing: `pp`

`Scimap` provides a suite of tools to preprocess the data for subsequent analysis.

| Function                                                                                      | Short Description                            |
|-----------------------------------------------------------------------------------------------|----------------------------------------------|
| [`sm.pp.mcmicro_to_scimap`](All%20Functions/A.%20Pre%20Processing/sm.pp.mcmicro_to_scimap.md) | `mcmicro` output to scimap compatible object |
| [`sm.pp.rescale`](All%20Functions/A.%20Pre%20Processing/sm.pp.rescale.md)                     | Manual/Auto gate based scaling of data       |


### Tools: `tl`

|                                                                                        |                                                                   |
|----------------------------------------------------------------------------------------|-------------------------------------------------------------------|
| [`sm.tl.cluster`](All%20Functions/B.%20Tools/sm.tl.cluster.md)                         | Cluster or sub-cluster single-cells using a variety of algorithms |
| [`sm.tl.foldchange`](All%20Functions/B.%20Tools/sm.tl.foldchange.md)                   | Compute foldchange in phenotypes between samples/ROIs             |
| [`sm.tl.phenotype_cells`](All%20Functions/B.%20Tools/sm.tl.phenotype_cells.md)         | Probability distribution based cell phenotyping                   |
| [`sm.tl.spatial_aggregate`](All%20Functions/B.%20Tools/sm.tl.spatial_aggregate.md)     | Aggregates of cell-types with local neighborhood                  |
| [`sm.tl.spatial_count`](All%20Functions/B.%20Tools/sm.tl.spatial_count.md)             | Distribution of cell-types with local neighborhood                |
| [`sm.tl.spatial_distance`](All%20Functions/B.%20Tools/sm.tl.spatial_distance.md)       | Computes nearest distance between all phenotypes for every cell   |
| [`sm.tl.spatial_expression`](All%20Functions/B.%20Tools/sm.tl.spatial_expression.md)   | Distribution of spatial expression with local neighborhood        |
| [`sm.tl.spatial_interaction`](All%20Functions/B.%20Tools/sm.tl.spatial_interaction.md) | cell–cell interactions analysis                                   |
| [`sm.tl.spatial_lda`](All%20Functions/B.%20Tools/sm.tl.spatial_lda.md)                 | Latent Dirichlet Allocation (LDA) modelling for spatial motifs    |
| [`sm.tl.spatial_pscore`](All%20Functions/B.%20Tools/sm.tl.spatial_pscore.md)           | Scoring proximity between user defined cell types                 |

### Plotting: `pl`

|                                                                                        |                                                                     |
|----------------------------------------------------------------------------------------|---------------------------------------------------------------------|
| [`sm.pl.gate_finder`](All%20Functions/C.%20Plotting/sm.pl.gate_finder.md)              | Overlays gating positivity on the image for manual gating           |
| [`sm.pl.image_viewer`](All%20Functions/C.%20Plotting/sm.pl.image_viewer)               | Opens the image with overlays on a `napari` browser                 |
| [`sm.pl.stacked_barplot`](All%20Functions/C.%20Plotting/sm.pl.stacked_barplot)         | Generate a stacked barplot from any two columns of categorical data |
| [`sm.pl.pie`](All%20Functions/C.%20Plotting/sm.pl.pie)                                 | Generate a pieplot of cell-type proportion or any categorical data  |
| [`sm.pl.foldchange`](All%20Functions/C.%20Plotting/sm.pl.foldchange)                   | vizualize foldchange in phenotypes between samples/ROIs             |
| [`sm.pl.voronoi`](All%20Functions/C.%20Plotting/sm.pl.voronoi)                         | Generate a voronoi diagram and color it with categorical data       |
| [`sm.pl.spatial_interaction`](All%20Functions/C.%20Plotting/sm.pl.spatial_interaction) | Heatmap of cell–cell interaction analysis                           |
| [`sm.pl.spatial_distance`](All%20Functions/C.%20Plotting/sm.pl.spatial_distance)       | Visualize distance between phenotypes                               |
| [`sm.pl.spatial_pscore`](All%20Functions/C.%20Plotting/sm.pl.spatial_pscore)           | Bar plot of the derived Spatial Proximity Scores                    |
|                                                                |                                                                                  |

### Helper Functions: `hl`

|                                                                                         |                                                                   |
|-----------------------------------------------------------------------------------------|-------------------------------------------------------------------|
| [`sm.hl.classify`](All%20Functions/D.%20Helper%20Functions/sm.hl.classify.md)           | Quickly classify cells based on pos/negativity of list of markers |
| [`sm.hl.scimap_to_csv`](All%20Functions/D.%20Helper%20Functions/sm.hl.scimap_to_csv.md) | Export scimap object to CSV                                       |
| [`sm.hl.add_roi_omero`](All%20Functions/D.%20Helper%20Functions/sm.hl.add_roi_omero.md) | Add ROI's extracted from Omero to Scimap object                   |
