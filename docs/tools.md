# Tools


``` python
import scimap as sm
```

# Preprocessing: `pp`

`Scimap` provides a suite of tools to preprocess the data for subsequent analysis.

|                                Function                                 |              Short Description               |
|:-----------------------------------------------------------------------:|:--------------------------------------------:|
| [`sm.pp.mcmicro_to_scimap (adata, ...)`](pp/sm.pp.mcmicro_to_scimap.md) | `mcmicro` output to scimap compatible object |
|            [`sm.pp.rescale (adata, ...)`](pp/sm.pp.rescale)             |    Manual/Auto gate based scaling of data    |
|                                                                         |                                              |

# Tools: `tl`

|                              Function                               |                Short Description                |
|:-------------------------------------------------------------------:|:-----------------------------------------------:|
| [`sm.tl.phenotype_cells (adata, ...)`](tl/sm.tl.phenotype_cells.md) | Probability distribution based cell phenotyping |


# Plotting: `pl`

|                          Function                           |                     Short Description                     |
|:-----------------------------------------------------------:|:---------------------------------------------------------:|
| [`sm.pl.gate_finder (adata, ...)`](pp/sm.pl.gate_finder.md) | Overlays gating positivity on the image for manual gating |
| [`sm.pl.image_viewer (adata, ...)`](pp/sm.pl.image_viewer)  |    Opens the image with overlays on a `napari` browser    |
|                                                             |                                                           |
