**scimap.hl.add_roi_omero**

!!! note "Function Call"
    `scimap.hl.add_roi_omero` (
      **adata, 
      roi, 
      x_coordinate='X_centroid', 
      y_coordinate='Y_centroid', 
      label='ROI'**)

**Short description**

Helper function to add ROI's extracted from Omero using this [script.](https://gist.github.com/Yu-AnChen/58754f960ccd540e307ed991bc6901b0) Please note that the function currently does not handle overlapping ROI's. Hence please make sure the ROI's are mutually exclusive.



**Parameters**

`adata` : AnnData object  

`roi` : data frame, required *(The default is None)*  
Pandas dataframe of ROI's that have been extracted from Omero. 

`x_coordinate` : float, required *(The default is 'X_centroid')*  
Column name containing the x-coordinates values.  

`y_coordinate` : float, required *(The default is 'Y_centroid')*  
Column name containing the y-coordinates values.

`label` : string, optional *(The default is 'ROI')*  
Key for the returned data, stored in `adata.obs`. 
 


**Returns**


Modified AnnData object. Check `adata.obs` for an additional column.


**Example**

```
# Load the ROI CSV file
roi = pd.read_csv('ROI/ROI_Z147_1_750.csv')

# Run the function
adata = sm.hl.add_roi_omero (adata, roi, label='aj_ROI')

```
