# scimap.pp.mcmicro_to_scimap

!!! abstract "Function"
    scimap.pp.mcmicro_to_scimap (image_path,remove_dna=True,remove_string_from_name=None,
                        log=True,drop_markers=None,random_sample=None,
                        CellId='CellID',split='X_centroid',custom_imageid=None,
                        min_cells=None)

#### Short description

### Parameters
image_path : string
      Location to the image file..
  adata : Ann Data Object
  marker_of_interest : string
      Marker for which gate is to be defined e.g. 'CD45'.
  from_gate : int, optional
      Start value gate of interest. The default is 6.
  to_gate : int, optional
      End value of the gate of interest. The default is 8.
  increment : float, optional
      Increments between the start and end values. The default is 0.1.
  markers : string, optional
      Additional markers to be included in the plot for evaluation. The default is None.
  channel_names : list, optional
      List of channels in the image in the exact order as image. The default is adata.uns['all_markers'].
  x : string, optional
      X axis coordinate column name in AnnData object. The default is 'X_position'.
  y : string, optional
      Y axis coordinate column name in AnnData object. The default is 'Y_position'.
  point_size : int, optional
      point size in the napari plot. The default is 10.
  image_id : string, optional
      The ID under 'ImageId' for the image of interest. The default is None.
  seg_mask : string, optional
      Location to the segmentation mask file. The default is None.

### Example
image_path = '/Users/aj/Desktop/ptcl_tma/image.tif'
sm.pl.gate_finder (image_path, adata, marker_of_interest='CD45',
               from_gate = 6, to_gate = 8, increment = 0.1,
               markers=['DNA10'], channel_names = 'default',
               x='X_position',y='Y_position',point_size=10,
               image_id= '77', seg_mask=None)
