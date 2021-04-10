**scimap.tl.foldchange**

!!! note "Function Call"
    `scimap.tl.foldchange` (
      **adata, 
      from_group, 
      to_group=None, 
      imageid='imageid', 
      phenotype='phenotype', 
      normalize=True, 
      subset_phenotype=None, 
      label='foldchange'**)

**Short description**

The `sm.tl.foldchange` function allows users to compute the foldchange (fc) in cell-type (phenotype) abundance 
between samples or ROI's. <br>
<br>
The reference sample or ROI needs to be passed via the `from_group` parameter. 
The column name of `from_group` should be passed via `imageid`. The function computes the fc 
to all other categories within the same `imageid` column. By default (can be turned off), the cell-abundance will
be normalized for the total number of cells within the sample/ROI to account for difference in area.
A `fisher-exact-test` is performed to compute the p-values. <br>
<br>
The results are stored in `.uns` section of the anndata object. 


**Parameters**

`adata` : AnnData Object  

`from_group` : list, required  
Pass in the name of the sample or ROI that will serve as a reference for calculating fold change.
If multiple sample names or ROI's are passed in as a list e.g. ['ROI1', 'ROI2''], please note that 
they will be combined for calculating the fold change. 

`to_group` : list, optional *(The default is `None`)*  
By default the reference sample/ROI passed via `from_group` will be compared to all other groups 
within the same column. However if users wish to restrict the comparision to a subset of 
samples/ROI's they can be passed though this paramenter as a list. e.g. ['ROI3', 'ROI4']. 

`imageid` : string, optional *(The default is `imageid`)*  
The column that contains the samples/ROI information.

`phenotype` : string, optional *(The default is `phenotype`)*  
The column that contains the cell-type/phenotype information.

`normalize` : bool, optional *(The default is `True`)*  
Inorder to account for the sample/ROI area, the cellular abundance is first normalized 
to the total number of cells within the respective sample/ROI. Please note if you pass values in 
`subset_phenotype`, the abundance normalization is restricted to the total cells of the 
cell types passed in via `subset_phenotype`.

`subset_phenotype` : list, optional *(The default is `None`)*  
If users are interested in only a subset of cell-types, the names of those can be passed in through 
this parameter. The data is subsetted to include only these cell types before computing foldchange.

`label` : string, optional *(The default is `foldchange_fc`)*  
Key for the returned data, stored in `adata.uns`. The default is 'foldchange'. The foldchange and p-values 
are returned seperately with the postfix `_fc` and `_pval` 


**Returns**
Updated `AnnData` object. Check `adata.uns['foldchange_fc']` and `adata.uns['foldchange_pval']` for results.

**Example**

```
# calculate the foldchange of cell-types between image_1 and all other images in the dataset.

adata = sm.tl.foldchange (adata, from_group='image_1', to_group=None, 
                              imageid='imageid', phenotype='phenotype',
                              normalize=True, 
                              subset_phenotype=['Tcells','Bcells','Macs'], 
                              label='foldchange')

```
