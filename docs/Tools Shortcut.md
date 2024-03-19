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

# Preprocessing (`pp`)

`Scimap` offers an array of preprocessing functions, meticulously designed to optimize high-dimensional datasets for advanced analysis.

| Function | Detailed Description |
|---|---|
| [`sm.pp.mcmicro_to_scimap`](Functions/pp/mcmicro_to_scimap.md) | Facilitates the transformation of `mcmicro` outputs into `scimap`-compatible formats, enabling the integration of microscopy data for comprehensive analysis. |
| [`sm.pp.log1p`](Functions/pp/log1p.md) | Applies a log1p transformation to the raw data within an AnnData object. |
| [`sm.pp.rescale`](Functions/pp/rescale.md) | Employs both manual and automated gating-based strategies for data rescaling, enhancing measurement scale sensitivity and dynamic range for optimal data representation. |
| [`sm.pp.combat`](Functions/pp/combat.md) | Implements an advanced batch correction algorithm to effectively neutralize batch-related variabilities, ensuring data consistency across different experimental batches. |

# Tools (`tl`)

`Scimap` presents a broad spectrum of analytical tools, each engineered to extract nuanced insights from single-cell datasets through sophisticated algorithms.

| Function | Description |
|---|---|
| [`sm.tl.phenotype_cells`](Functions/tl/phenotype_cells.md) | Leverages probability distribution models for precise cell phenotyping, enabling detailed characterization of cellular identities based on marker expressions. |
| [`sm.tl.cluster`](Functions/tl/cluster.md) | Provides a flexible clustering framework to delineate cellular subpopulations using a range of algorithms, facilitating the discovery of previously unrecognized cell types. |
| [`sm.tl.umap`](Functions/tl/umap.md) | Applies the UMAP algorithm for dimensionality reduction, affording a more interpretable visualization of complex single-cell data landscapes. |
| [`sm.tl.foldchange`](Functions/tl/foldchange.md) | Calculates fold changes in phenotype expressions between samples or Regions of Interest (ROIs), enabling quantitative comparisons of cellular characteristics. |
| [`sm.tl.spatial_distance`](Functions/tl/spatial_distance.md) | Computes the nearest distances between all phenotypes for each cell, offering insights into spatial arrangements and distributions within tissue contexts. |
| [`sm.tl.spatial_interaction`](Functions/tl/spatial_interaction.md) | Analyzes cell-cell interactions to uncover patterns of spatial organization and intercellular communication within microenvironments. |
| [`sm.tl.spatial_count`](Functions/tl/spatial_count.md) | Evaluates the distribution of cell types within local neighborhoods, providing a quantitative measure of cellular diversity and density. |
| [`sm.tl.spatial_lda`](Functions/tl/spatial_lda.md) | Utilizes Latent Dirichlet Allocation (LDA) modeling to identify spatial motifs, facilitating the understanding of complex spatial patterns and their biological implications. |
| [`sm.tl.spatial_expression`](Functions/tl/spatial_expression.md) | Investigates the distribution of spatial expression patterns within local neighborhoods, offering insights into the spatial heterogeneity of gene expression. |
| [`sm.tl.spatial_cluster`](Functions/tl/spatial_cluster.md) | Identifies clusters based on spatial expression patterns, enabling the elucidation of spatially defined cellular networks and communities. |
| [`sm.tl.spatial_pscore`](Functions/tl/spatial_pscore.md) | Scores the proximity between predefined cell types, quantifying the spatial relationships and potential functional interactions between distinct cellular populations. |
| [`sm.tl.spatial_aggregate`](Functions/tl/spatial_aggregate.md) | Summarizes aggregates of cell types within local neighborhoods, providing a macroscopic view of cellular organization and tissue architecture. |
| [`sm.tl.spatial_similarity_search`](Functions/tl/spatial_similarity_search.md) | Searches for regions within and across images that exhibit similar spatial patterns, aiding in the identification of recurring spatial motifs and their biological relevance. |

# Plotting (`pl`)

`Scimap` incorporates a comprehensive suite of plotting functions designed for the intuitive visualization and interpretation of spatial and phenotypic data.

| Function | Description |
|---|---|
| [`sm.pl.image_viewer`](Functions/pl/image_viewer.md) | Integrates with `napari` to offer an interactive platform for enhanced image viewing and annotation with data overlays. |
| [`sm.pl.addROI_image`](Functions/pl/addROI_image.md) | Facilitates the addition of Regions of Interest (ROIs) through `napari`, enriching spatial analyses with precise locational data. |
| [`sm.pl.gate_finder`](Functions/pl/gate_finder.md) | Aids in the manual gating process by overlaying marker positivity on images, simplifying the identification and analysis of cellular subsets. |
| [`sm.pl.heatmap`](Functions/pl/heatmap.md) | Creates heatmaps to visually explore marker expression or feature distributions across different groups. |
| [`sm.pl.markerCorrelation`](Functions/pl/markerCorrelation.md) | Computes and visualizes the correlation among selected markers. |
| [`sm.pl.groupCorrelation`](Functions/pl/groupCorrelation.md) | Calculates and displays the correlation between the abundances of groups across user defined conditions. |
| [`sm.pl.distPlot`](Functions/pl/distPlot.md) | Generates distribution plots for specific markers, allowing for the visual comparison of marker expression across different conditions or cell types. |
| [`sm.pl.densityPlot2D`](Functions/pl/densityPlot2D.md) | Creates two-dimensional density plots of marker expressions, facilitating the visualization of expression patterns and densities in a spatial context. |
| [`sm.pl.cluster_plots`](Functions/pl/cluster_plots.md) | Provides a meta-function that outputs a combination of UMAP, heatmap, and ranked markers for each group, offering a comprehensive view of clustering results. |
| [`sm.pl.umap`](Functions/pl/umap.md) | Overlays markers on UMAP projections, enhancing the interpretation of dimensional reduction analyses with annotated data points. |
| [`sm.pl.foldchange`](Functions/pl/foldchange.md) | Visualizes fold changes in phenotypes between samples or ROIs, enabling the graphical comparison of cellular expression profiles. |
| [`sm.pl.spatial_scatterPlot`](Functions/pl/spatial_scatterPlot.md) | Produces scatter plots of spatially resolved data, illustrating the distribution and organization of cells within tissue sections. |
| [`sm.pl.spatial_distance`](Functions/pl/spatial_distance.md) | Visualizes the spatial distances between phenotypes, providing insights into the physical separation and clustering of cell types. |
| [`sm.pl.spatial_interaction`](Functions/pl/spatial_interaction.md) | Displays heatmaps of cell-cell interaction analyses, highlighting the complex interplay between different cellular populations. |
| [`sm.pl.spatialInteractionNetwork`](Functions/pl/spatialInteractionNetwork.md) | Displays  cell-cell interaction analyses as a Network plot. |
| [`sm.pl.spatial_pscore`](Functions/pl/spatial_pscore.md) | Generates bar plots of Spatial Proximity Scores, quantifying and visualizing the proximity between selected cell types within a spatial context. |
| [`sm.pl.stacked_barplot`](Functions/pl/stacked_barplot.md) | Creates stacked barplots from any two columns of categorical data, offering a clear visualization of proportions and relationships between categories. |
| [`sm.pl.pie`](Functions/pl/pie.md) | Produces pie charts to represent the proportions of cell types or any categorical data, facilitating the quick assessment of composition within datasets. |
| [`sm.pl.voronoi`](Functions/pl/voronoi.md) | Generates Voronoi diagrams colored by categorical data, providing a unique visual representation of spatial distributions and territories. |

# Helper Functions (`hl`)

`Scimap` also features a collection of helper functions designed to facilitate efficient data manipulation and enhance the analytical workflow.

| Function | Description |
|---|---|
| [`sm.hl.classify`](Functions/hl/classify.md) | Streamlines the classification of cells based on the positivity or negativity of specified markers, simplifying the assignment of cellular identities. |
| [`sm.hl.rename`](Functions/hl/rename.md) | Enables quick and efficient renaming within data columns through the application of dictionary-based mappings, enhancing data clarity and consistency. |
| [`sm.hl.addROI_omero`](Functions/hl/addROI_omero.md) | Allows for the seamless integration of Regions of Interest (ROIs) extracted from Omero into `scimap` objects, bridging imaging and analytical platforms. |
| [`sm.hl.dropFeatures`](Functions/hl/dropFeatures.md) | Provides a convenient function for subsetting the `adata` object by removing specified features, aiding in the focus on relevant data. |
| [`sm.hl.animate`](Functions/hl/animate.md) | Creates animated scatter plots transitioning from embedding to physical location, offering dynamic visual insights into spatial data relationships. |
| [`sm.hl.merge_adata_obs`](Functions/hl/merge_adata_obs.md) | Facilitates the merging of multiple AnnData objects, ensuring cohesive analysis across disparate datasets. |
| [`sm.hl.scimap_to_csv`](Functions/hl/scimap_to_csv.md) | Enables the export of `scimap` objects to CSV format, providing flexibility in data sharing and further analysis with external tools. |


