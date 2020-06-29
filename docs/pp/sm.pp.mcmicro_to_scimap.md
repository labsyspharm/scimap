# scimap.pp.mcmicro_to_scimap

!!! note "Function"
    `scimap.pp.mcmicro_to_scimap` (
      image_path,
      remove_dna=True,
      remove_string_from_name=None,
      log=True,
      drop_markers=None,
      random_sample=None,
      CellId='CellID',
      split='X_centroid',
      custom_imageid=None,
      min_cells=None)

##### Short description

Function helps to convert `mcmicro` output to AnnData object, ready for analysis.

##### Parameters

`image_path` : list  
    List of path to the image or images. Each Image should have a unique path supplied.  
`remove_dna` : bool, optional (The default is True)  
    Remove the DNA channels from the final output. Looks for channels with the string 'dna' in it.    
`remove_string_from_name` : string, optional (The default is None)  
    Used to celan up channel names. If a string is given, that particular string will be removed from all marker names.  
    If multiple images are passed, just use the string that appears in the first image.  
`log` : bool, optional (The default is True)  
    Log the data (log1p transformation will be applied).  
`drop_markers` : list, optional (The default is None)  
    List of markers to drop from the analysis. e.g. ["CD3D", "CD20"].  
`random_sample` : int, optional (he default is None)  
    Randomly sub-sample the data with the desired number of cells.  
`CellId` : string, optional (The default is CellID)  
    Name of the column that contains the cell ID.  
`split` : string, optional (The default is 'X_centroid')  
    To split the CSV into counts table and meta data, pass in the name of the column  
    that immediately follows the marker quantification.  
`custom_imageid`: string, optional (Default is the name of the CSV file)  
    Pass a user defined Image ID.  
`min_cells`: int, optional (Default is None)  
    If these many cells are not in the image, the image will be dropped.  
    Particularly useful when importing multiple images.  

##### Example

```
image_path = ['/Users/aj/whole_sections/PTCL1_450.csv',
             '/Users/aj/whole_sections/PTCL2_552.csv']

adata = sm.pp.mcmicro_to_scimap (image_path, drop_markers= ['CD21', 'ACTIN'], random_sample=5000)
```
