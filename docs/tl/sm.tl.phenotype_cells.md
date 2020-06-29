**scimap.tl.phenotype_cells**

!!! note "Function Call"
    `scimap.tl.phenotype_cells` (
      **adata,
      phenotype,
      gate = 0.5,
      label="phenotype",
      unique_id='imageid',
      pheno_threshold_percent=None,
      pheno_threshold_abs=None**)

**Short description**

The phenotyping function takes in the `scaled data` and a prior knowledge based `phenotype workflow` file to assign phenotype to each cell in the dataset.

*Phenotype workflow file description:*
An example of the `phenotype_workflow.csv` can be found [here](../).

**Parameters**

`adata` : AnnData object  

`gate` : dataframe, optional *(The default is None)*   
DataFrame with first column as markers and second column as the gate values in log1p scale.  

`failed_markers` : list, optional *(The default is None)*  
list of markers that are not expressed at all in any cell. pass in as ['CD20', 'CD3D'].  

`method` : string, optional *(The default is 'all')*  
Two available option are- 'all' or 'by_image'. In the event that multiple images were loaded in with distinct 'imageid', users have the option to scale all data together or each image independently. Please be aware of batch effects when passing 'all' with multiple images.  

`save_fig` : boolian, optional *(The default is False)*  
If True, the gates identified by the GMM method will be saved in a subdirectory within your working directory.  


**Returns**

`AnnData` object with the rescaled data `adata.X`

**Example**

```
manual_gate = pd.DataFrame({'marker': ['CD3D', 'KI67'], 'gate': [7, 8]})
adata = sm.pp.rescale (adata, gate=manual_gate, failed_markers=['CD20', 'CD21'])
```
