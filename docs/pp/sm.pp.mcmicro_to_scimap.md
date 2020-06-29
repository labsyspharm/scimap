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

#### Short description

Function helps to convert `mcmicro` output to AnnData object, ready for analysis.

### Parameters

`image_path` : list  
    List of path to the image or images. Each Image should have a unique path supplied.  
`remove_dna` : bool, optional  
    Remove the DNA channels from the final output. Looks for channels with the string 'dna' in it. The default is True.  
`remove_string_from_name` : string, optional  
    Used to celan up channel names. If a string is given, that particular string will be removed from all marker names.  
    If multiple images are passed, just use the string that appears in the first image. The default is None.  
`log` : bool, optional  
    Log the data (log1p transformation will be applied). The default is True.  
`drop_markers` : list, optional  
    List of markers to drop from the analysis. e.g. ["CD3D", "CD20"]. The default is None.  
`random_sample` : int, optional  
    Randomly sub-sample the data with the desired number of cells. The default is None.  
`CellId` : string, optional  
    Name of the column that contains the cell ID. The default is CellID.  
`split` : string, optional  
    To split the CSV into counts table and meta data, pass in the name of the column  
    that immediately follows the marker quantification. The default is 'X_centroid'.  
`custom_imageid`: string, optional  
    Pass a user defined Image ID. By default the name of the CSV file is used.  
`min_cells`: int, optional  
    If these many cells are not in the image, the image will be dropped.  
    Particularly useful when importing multiple images.  

### Example

```
image_path = ['/Users/aj/whole_sections/PTCL1_450.csv',
             '/Users/aj/whole_sections/PTCL2_552.csv']

adata = sm.pp.mcmicro_to_scimap (image_path, drop_markers= ['CD21', 'ACTIN'], random_sample=5000)
```
