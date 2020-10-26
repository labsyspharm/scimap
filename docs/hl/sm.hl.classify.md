**scimap.hl.classify**

!!! note "Function Call"
    `scimap.hl.classify` (
      **adata, 
      pos=None, neg=None, 
      classify_label='passed_classify', 
      phenotype='phenotype', subclassify_phenotype=None, 
      collapse_failed=True, label="classify"**)

**Short description**

The function allows users to annotate cells based on positivity/negativity of defined markers. Users can classify
the entire data or a subset of data that has been previously phenotyped or clustered.

**Parameters**

`adata` : AnnData object  

`pos`: list *(The default is None)*  
Pass a list of markers that should be expressed in the resultant cells.

`neg` : list, optional *(The default is None)*   
Pass a list of markers that should not be expressed in the resultant cells.

`classify_label` : string, optional *(The default is 'passed_classify'.)*  
Provide a name for the calssified cells. 

`phenotype` : string, optional *(The default is 'phenotype'.)*  
Column name of the column containing the phenotype information. 
This is important if `subclassify_phenotype` or `collapse_failed` arguments are used.
        
`subclassify_phenotype` : list, optional *(The default is None.)*  
If only a subset of phenotypes require to classified, pass the name of those phenotypes as a list
through this argument. 

`collapse_failed` : bool, optional *(The default is True)*  
If set to true, the cells that were not classified based on the given criteria will be
binned into a single category named 'failed_classify'. When False, the phenotype
inforamation for other cells will be borrowed from the `phenotype` argument. 

`label` : string, optional *(The default is "classify")*  
Key for the returned data, stored in `adata.obs`. 


**Returns**

`adata` : AnnData  
Updated AnnData Object.

**Example**

```
# Classify all cells with both pos and neg markers (Identify cytotoxic T-cells)
adata = sm.hl.classify(adata, pos=['CD3D','CD8A'], neg=['ASMA'])

# Classify specific sub-types of cells
adata = sm.hl.classify(adata, pos=['CD3D','FOXP3'], neg=['ASMA'], subclassify_phenotype=['T cells','Regulatory T cells'])

# Classify specific sub-types of cells and borrow labels from another column
adata = sm.hl.classify(adata, pos=['CD3D'], neg=['ASMA'], subclassify_phenotype=['T cells'], collapse_failed=False, phenotype='phenotype')

```
