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
| [`sm.pp.mcmicro_to_scimap`](All%20Functions/A.%20Pre%20Processing/sm.pp.mcmicro_to_scimap.md) | `mcmicro` output to scimap compatible object |
| [`sm.pp.rescale`](All%20Functions/A.%20Pre%20Processing/sm.pp.rescale.md)                     | Manual/Auto gate based scaling of data       |
| [`sm.pp.combat`](All%20Functions/A.%20Pre%20Processing/sm.pp.combat.md)                       | Batch correction tool                        |


### Tools: `tl`

|                                                                                                    |                                                                   |
|----------------------------------------------------------------------------------------------------|-------------------------------------------------------------------|
| [`sm.tl.umap`](All%20Functions/B.%20Tools/sm.tl.umap.md)                                           | Dimensionality Reduction using UMAP                               |
| [`sm.tl.cluster`](All%20Functions/B.%20Tools/sm.tl.cluster.md)                                     | Cluster or sub-cluster single-cells using a variety of algorithms |
| [`sm.tl.foldchange`](All%20Functions/B.%20Tools/sm.tl.foldchange.md)                               | Compute foldchange in phenotypes between samples/ROIs             |
| [`sm.tl.phenotype_cells`](All%20Functions/B.%20Tools/sm.tl.phenotype_cells.md)                     | Probability distribution based cell phenotyping                   |
| [`sm.tl.spatial_aggregate`](All%20Functions/B.%20Tools/sm.tl.spatial_aggregate.md)                 | Aggregates of cell-types with local neighborhood                  |
| [`sm.tl.spatial_count`](All%20Functions/B.%20Tools/sm.tl.spatial_count.md)                         | Distribution of cell-types with local neighborhood                |
| [`sm.tl.spatial_distance`](All%20Functions/B.%20Tools/sm.tl.spatial_distance.md)                   | Computes nearest distance between all phenotypes for every cell   |
| [`sm.tl.spatial_expression`](All%20Functions/B.%20Tools/sm.tl.spatial_expression.md)               | Distribution of spatial expression with local neighborhood        |
| [`sm.tl.spatial_interaction`](All%20Functions/B.%20Tools/sm.tl.spatial_interaction.md)             | cell–cell interactions analysis                                   |
| [`sm.tl.spatial_lda`](All%20Functions/B.%20Tools/sm.tl.spatial_lda.md)                             | Latent Dirichlet Allocation (LDA) modelling for spatial motifs    |
| [`sm.tl.spatial_pscore`](All%20Functions/B.%20Tools/sm.tl.spatial_pscore.md)                       | Scoring proximity between user defined cell types                 |
| [`sm.tl.spatial_similarity_search`](All%20Functions/B.%20Tools/sm.tl.spatial_similarity_search.md) | Search for similar looking regions within and across images       |

### Plotting: `pl`

|                                                                                           |                                                                     |
|-------------------------------------------------------------------------------------------|---------------------------------------------------------------------|
| [`sm.pl.umap`](All%20Functions/C.%20Plotting/sm.pl.umap.md)                               | Overlays markers on UMAP                                            |
| [`sm.pl.gate_finder`](All%20Functions/C.%20Plotting/sm.pl.gate_finder.md)                 | Overlays gating positivity on the image for manual gating           |
| [`sm.pl.image_viewer`](All%20Functions/C.%20Plotting/sm.pl.image_viewer.md)               | Opens the image with overlays on a `napari` browser                 |
| [`sm.pl.stacked_barplot`](All%20Functions/C.%20Plotting/sm.pl.stacked_barplot.md)         | Generate a stacked barplot from any two columns of categorical data |
| [`sm.pl.pie`](All%20Functions/C.%20Plotting/sm.pl.pie.md)                                 | Generate a pieplot of cell-type proportion or any categorical data  |
| [`sm.pl.foldchange`](All%20Functions/C.%20Plotting/sm.pl.foldchange.md)                   | vizualize foldchange in phenotypes between samples/ROIs             |
| [`sm.pl.voronoi`](All%20Functions/C.%20Plotting/sm.pl.voronoi.md)                         | Generate a voronoi diagram and color it with categorical data       |
| [`sm.pl.spatial_interaction`](All%20Functions/C.%20Plotting/sm.pl.spatial_interaction.md) | Heatmap of cell–cell interaction analysis                           |
| [`sm.pl.spatial_distance`](All%20Functions/C.%20Plotting/sm.pl.spatial_distance.md)       | Visualize distance between phenotypes                               |
| [`sm.pl.spatial_pscore`](All%20Functions/C.%20Plotting/sm.pl.spatial_pscore.md)           | Bar plot of the derived Spatial Proximity Scores                    |
| [`sm.pl.addROI_image`](All%20Functions/C.%20Plotting/sm.pl.addROI_image.md)               | Add ROI's with  `napari` browser                                    |
|                                                                                           |                                                                     |

### Helper Functions: `hl`

|                                                                                         |                                                                    |
|-----------------------------------------------------------------------------------------|--------------------------------------------------------------------|
| [`sm.hl.classify`](All%20Functions/D.%20Helper%20Functions/sm.hl.classify.md)           | Quickly classify cells based on pos/negativity of list of markers  |
| [`sm.hl.rename`](All%20Functions/D.%20Helper%20Functions/sm.hl.rename.md)               | Quickly rename values within columns based on `dict` mapping       |
| [`sm.hl.animate`](All%20Functions/D.%20Helper%20Functions/sm.hl.animate.md)             | Create a animated scatter plot of `embedding -> physical location` |
| [`sm.hl.scimap_to_csv`](All%20Functions/D.%20Helper%20Functions/sm.hl.scimap_to_csv.md) | Export scimap object to CSV                                        |
| [`sm.hl.addROI_omero`](All%20Functions/D.%20Helper%20Functions/sm.hl.addROI_omero.md)   | Add ROI's extracted from Omero to Scimap object                    |
| [`sm.hl.dropFeatures`](All%20Functions/D.%20Helper%20Functions/sm.hl.dropFeatures.md)   | Handy Function to subset the `adata` object                        |

