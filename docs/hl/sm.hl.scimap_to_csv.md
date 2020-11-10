**scimap.hl.scimap_to_csv**

!!! note "Function Call"
    `scimap.hl.scimap_to_csv` (
      **adata, 
       data_type='raw'**)

**Short description**

Helper function to convert andata object to a csv file with the expression matrix and available metadata.

**Parameters**

`adata` : AnnData object  

`pos`: list *(The default is None)*  
Pass a list of markers that should be expressed in the resultant cells.  

`data_type` : string, optional *(The default is 'raw')*  
Three options are available:  
1) 'raw' - The raw data will be returned.  
2) 'log' - The raw data converted to log scale using `np.log1p` will be returned.  
3) 'scaled' - If you have scaled the data using the `sm.pp.rescale`, that will be
returned. Please note, if you have not scaled the data, whatever is within
`adata.X` will be returned.  


**Returns**


A single dataframe containing the expression and metadata.


**Example**

```
data = sm.hl.scimap_to_csv (adata, data_type='raw')

```
