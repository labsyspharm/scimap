**scimap.pl.image_viewer**

!!! note "Function Call"
    `scimap.pl.image_viewer` (
      **image_path, 
      adata, 
      overlay=None, 
      overlay_category=None, 
      markers=None, 
      channel_names='default', 
      x_coordinate='X_centroid', 
      y_coordinate='Y_centroid', 
      point_size=10, 
      point_color=None, 
      imageid='imageid', 
      subset=None, 
      seg_mask=None**)

**Short description**

This function spins up a `napari` instance to view an image. For large images, it is advised to load only a few markers at a time by using the `markers` parameter. Any categorical data stored in `adata.obs` can also be overlayed on the image using the `overlay` parameter.

**Parameters**

`image_path` : string  
Location to the image file.  

`seg_mask`: string *(The default is None)*  
Location to the segmentation mask file.  

`adata` : AnnData object  

`overlay` : string, optional *(The default is None)*  
Name of the column with any categorical data such as phenotypes or clusters.  

`overlay_category` : list, optional *(The default is None)*  
If only specfic categories within the overlay column is needed, pass their names as a list. If None, all categories will be used.  

`markers` : list, optional *(The default is None)*  
Markers to be included. If none, all markers will be displayed.  

`imageid` : string, optional *(The default is `imageid`)*   
Column name of the column containing the image id. 

`subset` : string, optional  *(The default is None)*  
imageid of a single image to be subsetted for analyis. Only useful when multiple images are being analyzed together. 

`channel_names` : list, optional *(The default is `adata.uns['all_markers']`)*  
List of channels in the image in the exact order as image.  

`x_coordinate` : string, optional *(The default is 'X_centroid')*  
X axis coordinate column name in AnnData object.  

`y_coordinate` : string, optional *(The default is 'Y_centroid')*  
Y axis coordinate column name in AnnData object.  

`point_size` : int, optional *(The default is 10)*  
point size in the napari plot.  

**Returns**

`napari` image viewer loads the image.

**Example**

```
image_path = '/Users/aj/Desktop/ptcl_tma/image.tif'
sm.pl.image_viewer (image_path, adata, overlay='phenotype',overlay_category=None,
            markers=['CD31','CD19','CD45','CD163','FOXP3'],
            point_size=7, point_color='white')
```
