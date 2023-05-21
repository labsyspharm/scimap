---
title: Shortcut
description: Spatial Single-Cell Analysis Toolkit
hide:
  - navigation
---

# Tools Shortcut

```python
import scimap as sm
```

### Preprocessing: `pp`

`Scimap` provides a suite of tools to preprocess the data for subsequent analysis.

| Function                                                                                      | Short Description                            |
|-----------------------------------------------------------------------------------------------|----------------------------------------------|
| [`sm.pp.mcmicro_to_scimap`](Functions/pp/mcmicro_to_scimap.md) | `mcmicro` output to scimap compatible object |
| [`sm.pp.rescale`](Functions/pp/rescale.md)                     | Manual/Auto gate based scaling of data       |
| [`sm.pp.combat`](Functions/pp/combat.md)                       | Batch correction tool                        |


### Tools: `tl`

|                                                                                                    |                                                                   |
|----------------------------------------------------------------------------------------------------|-------------------------------------------------------------------|
| [`sm.tl.phenotype_cells`](Functions/tl/phenotype_cells.md)                     | Probability distribution based cell phenotyping                   |
| [`sm.tl.cluster`](Functions/tl/cluster.md)                                     | Cluster or sub-cluster single-cells using a variety of algorithms |
| [`sm.tl.umap`](Functions/tl/umap.md)                                           | Dimensionality Reduction using UMAP                               |
| [`sm.tl.foldchange`](Functions/tl/foldchange.md)                               | Compute foldchange in phenotypes between samples/ROIs             |
| [`sm.tl.spatial_distance`](Functions/tl/spatial_distance.md)                   | Computes nearest distance between all phenotypes for every cell   |
| [`sm.tl.spatial_interaction`](Functions/tl/spatial_interaction.md)             | cell–cell interactions analysis                                   |
| [`sm.tl.spatial_count`](Functions/tl/spatial_count.md)                         | Distribution of cell-types with local neighborhood                |
| [`sm.tl.spatial_lda`](Functions/tl/spatial_lda.md)                             | Latent Dirichlet Allocation (LDA) modelling for spatial motifs    |
| [`sm.tl.spatial_expression`](Functions/tl/spatial_expression.md)               | Distribution of spatial expression with local neighborhood        |
| [`sm.tl.spatial_cluster`](Functions/tl/spatial_cluster.md)               | Distribution of spatial expression with local neighborhood        |
| [`sm.tl.spatial_pscore`](Functions/tl/spatial_pscore.md)                       | Scoring proximity between user defined cell types                 |
| [`sm.tl.spatial_aggregate`](Functions/tl/spatial_aggregate.md)                 | Aggregates of cell-types with local neighborhood                  |
| [`sm.tl.spatial_similarity_search`](Functions/tl/spatial_similarity_search.md) | Search for similar looking regions within and across images       |

### Plotting: `pl`

|                                                                                           |                                                                     |
|-------------------------------------------------------------------------------------------|---------------------------------------------------------------------|
| [`sm.pl.image_viewer`](Functions/pl/image_viewer.md)               | Opens the image with overlays on a `napari` browser                 |
| [`sm.pl.addROI_image`](Functions/pl/addROI_image.md)               | Add ROI's with  `napari` browser                                    |
| [`sm.pl.gate_finder`](Functions/pl/gate_finder.md)                 | Overlays gating positivity on the image for manual gating           |
| [`sm.pl.distPlot`](Functions/pl/distPlot.md)               | Distribution plot of a given marker                                    |
| [`sm.pl.densityPlot2D`](Functions/pl/densityPlot2D.md)               | 2D density plotting of marker expression                                    |
| [`sm.pl.cluster_plots`](Functions/pl/cluster_plots.md)               | Meta function that outputs umap, heatmap and ranked makers for each group                                     |
| [`sm.pl.umap`](Functions/pl/umap.md)                               | Overlays markers on UMAP                                            |
| [`sm.pl.foldchange`](Functions/pl/foldchange.md)                   | vizualize foldchange in phenotypes between samples/ROIs             |
| [`sm.pl.spatial_scatterPlot`](Functions/pl/spatial_scatterPlot.md)       | Scatter plots of spatially resolved data                               |
| [`sm.pl.spatial_distance`](Functions/pl/spatial_distance.md)       | Visualize distance between phenotypes                               |
| [`sm.pl.spatial_interaction`](Functions/pl/spatial_interaction.md) | Heatmap of cell–cell interaction analysis                           |
| [`sm.pl.spatial_pscore`](Functions/pl/spatial_pscore.md)           | Bar plot of the derived Spatial Proximity Scores                    |
| [`sm.pl.stacked_barplot`](Functions/pl/stacked_barplot.md)         | Generate a stacked barplot from any two columns of categorical data |
| [`sm.pl.pie`](Functions/pl/pie.md)                                 | Generate a pieplot of cell-type proportion or any categorical data  |
| [`sm.pl.voronoi`](Functions/pl/voronoi.md)                         | Generate a voronoi diagram and color it with categorical data       |
|                                                                                           |                                                                     |

### Helper Functions: `hl`

|                                                                                         |                                                                    |
|-----------------------------------------------------------------------------------------|--------------------------------------------------------------------|
| [`sm.hl.classify`](Functions/hl/classify.md)           | Quickly classify cells based on pos/negativity of list of markers  |
| [`sm.hl.rename`](Functions/hl/rename.md)               | Quickly rename values within columns based on `dict` mapping       |
| [`sm.hl.addROI_omero`](Functions/hl/addROI_omero.md)   | Add ROI's extracted from Omero to Scimap object                    |
| [`sm.hl.dropFeatures`](Functions/hl/dropFeatures.md)   | Handy Function to subset the `adata` object                        |
| [`sm.hl.animate`](Functions/hl/animate.md)             | Create a animated scatter plot of `embedding -> physical location` |
| [`sm.hl.scimap_to_csv`](Functions/hl/scimap_to_csv.md) | Export scimap object to CSV                                        |

