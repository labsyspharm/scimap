# 0.11.0 (2021-01-30)

 - Included `spatial_pscore` function
 
# 0.10.0 (2020-11-27)

 - Included `stacked_barplot` function to generate a stacked barplot from any two cloumns
   
# 0.9.0 (2020-11-21)

 - Updated `pl.image_viewer` and `pl.gate_finder` functions - 
   Implemeted Zarr functionality for napari viz.
   
# 0.8.7 (2020-11-18)

 - Updated `pl.spatial_interaction` function to include two additional parameters - 
   `subset_phenotype` and `subset_neighbour_phenotype`. 

# 0.8.3 (2020-11-16)

 - Added `hl.add_roi` function. Used to incorporate ROI's extracted from Omero into the scimap object.
 
 
# 0.8.0 (2020-11-09)

 - Updated `pp.mcmicro_to_scimap` function. Added a new parameter `unique_CellId`
 - Added a helper function `scimap_to_csv` to save the andata object as a CSV.
 - Added documentation and tests for `scimap_to_csv`
 
# 0.7.10 (2020-10-30)

 - Updated `pp.rescale` function. If a gate is included in the `manual_gate.csv` 
   file but no gate value is provided, the algorithm simply scales the data between
   0-1 without changing the undelying structure.
   
# 0.7.6 (2020-10-27)

 - Updated `hl.spatial_distance` to include option to convert to 
   log scale and also pass multiple `distance_to` parameter.
 
# 0.7.5 (2020-10-27)

 - Updated `hl.classify` to improve speed
 
# 0.7.3 (2020-10-27)

 - Addition of binary view in `pl.spatial_interaction`

 
# 0.7.2 (2020-10-26)

 - Addition of `hl.classify` function.
 - Documentation for `hl.classify` function.
 - Readme file modification
