# Tools


``` python
import scimap as sm
```

### Preprocessing: `pp`

`Scimap` provides a suite of tools to preprocess the data for subsequent analysis.

|                                  Function                                   |              Short Description               |
|:---------------------------------------------------------------------------:|:--------------------------------------------:|
| [`sm.pp.mcmicro_to_scimap (data_path, ...)`](pp/sm.pp.mcmicro_to_scimap.md) | `mcmicro` output to scimap compatible object |
|             [`sm.pp.rescale (adata, ...)`](../pp/sm.pp.rescale)             |    Manual/Auto gate based scaling of data    |
|                                                                             |                                              |

### Tools: `tl`

|                                                                     |                                                 |
|:-------------------------------------------------------------------:|:-----------------------------------------------:|
| [`sm.tl.phenotype_cells (adata, ...)`](tl/sm.tl.phenotype_cells.md) | Probability distribution based cell phenotyping |
| [`sm.tl.spatial_count (adata, ...)`](tl/sm.tl.spatial_count.md)     | Distribution of cell-types with local neighbourhood |
| [`sm.tl.spatial_expression (adata, ...)`](tl/sm.tl.spatial_expression.md)     | Distribution of spatial expression with local neighbourhood |
| [`sm.tl.spatial_aggregate (adata, ...)`](tl/sm.tl.spatial_aggregate.md)     | Agrregates of cell-types with local neighbourhood |
| [`sm.tl.cluster (adata, ...)`](tl/sm.tl.cluster.md)     | Cluster or sub-cluster single-cells using a variety of algorithms |

### Plotting: `pl`

|                                                                    |                                                           |
|:------------------------------------------------------------------:|:---------------------------------------------------------:|
|  [`sm.pl.gate_finder (image_path, ...)`](pl/sm.pl.gate_finder.md)  | Overlays gating positivity on the image for manual gating |
| [`sm.pl.image_viewer (image_path, ...)`](../pl/sm.pl.image_viewer) |    Opens the image with overlays on a `napari` browser    |
|                                                                    |                                                           |
