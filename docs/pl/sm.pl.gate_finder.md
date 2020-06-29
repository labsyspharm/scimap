**scimap.pl.gate_finder**

!!! note "Function Call"
    `scimap.pl.gate_finder` (
      **image_path,
      adata,
      marker_of_interest,
      from_gate = 6, to_gate = 8, increment = 0.1,
      markers = None,
      channel_names = 'default',
      x='X_centroid',
      y='Y_centroid',
      point_size=10,
      image_id=None,
      seg_mask=None**)

**Short description**

The function helps to identify manual gates for each marker by overlaying the predicted positive cells on the image. For each marker `from_gate`, `to_gate` and an `increment` can be passed to identify predicted positive cells for multiple gates. All of these are overlayed on the image to identify the best gate visually.

**Parameters**

`image_path` : string  
Location to the image file.  

`adata` : AnnData Object  

`marker_of_interest` : string  
Marker for which gate is to be defined e.g. 'CD45'.  

`from_gate` : int, optional *(The default is 6)*  
Start value gate of interest.  

`to_gate` : int, optional *(The default is 8)*  
End value of the gate of interest.  

`increment` : float, optional *(The default is 0.1)*  
Increments between the start and end values.  

`markers` : string, optional *(The default is None)*  
Additional markers to be included in the image for evaluation.  

`channel_names` : list, optional *(The default is `adata.uns['all_markers']`)*  
List of channels in the image in the exact order as image.  

`x` : string, optional *(The default is 'X_centroid')*  
X axis coordinate column name in AnnData object.  

`y` : string, optional *(The default is 'Y_centroid')*  
Y axis coordinate column name in AnnData object.  

`point_size` : int, optional *(The default is 10)*  
point size in the napari plot.  

`image_id` : string, optional *(The default is None)*  
The ID under 'ImageId' for the image of interest.  

`seg_mask` : string, optional *(The default is None)*  
Location to the segmentation mask file.  


**Returns**

`napari` image viewer loads the image with the predicted cells that are positive for the given marker and gate.

**Example**

```
image_path = '/Users/aj/Desktop/ptcl_tma/image.tif'
sm.pl.gate_finder (image_path, adata, marker_of_interest='CD45',
             from_gate = 6, to_gate = 8, increment = 0.1,
             markers=['DNA10', 'CD20'], image_id= '77', seg_mask=None)
```
