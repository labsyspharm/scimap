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
An example of the `phenotype_workflow.csv` can be found [here](https://github.com/ajitjohnson/scimap/blob/master/scimap/tests/_data/phenotype_workflow.csv).  

The `phenotype_workflow` accepts six categories of gating strategy for performing phenotyping.

- all_pos
- all_neg
- any_pos
- any_neg
- pos
- neg

`all_pos`- All of the defined markers should be positive.  
`all_neg`- All of the defined markers should be negative.  
`any_pos`- Any one of the defined marker is sufficient to be positive. (e.g) For defining macrophages, one could use a strategy in which a cell is defined as a macrophage if any of `CD68, CD163 or CD206` is positive.  
`any_neg`- Any of the defined marker is negative.  
`pos`- A given marker is positive. If this argument is passed to multiple markers. (e.g) If regulatory T cell is defined as `CD4+`, `FOXP3+` by passing `pos` to each the markers and the algorithm finds that for a few cells one of the two is not, the algorithm will assign the cell as likely-regulatory T cell and will allow the user to make the decision later.  
`neg`- A given marker is negative.  
*It is always advised to use positive markers over negative markers*  

**Parameters**

`adata` : AnnData Object  

`phenotype` : DataFrame  
A gating strategy for phenotyping the cells.  

`gate` : int, optional *(The default is 0.5)*  
By default rescale function, scales the data such that values above 0.5 are considered positive cells.  

`label` : string, optional *(The default is "phenotype")*  
Name the column under which the final phenotype calling will be saved.  

`unique_id` : string, optional *(The default is 'imageid')*  
Name of the column that contains the unique imageid. This is only utilized when the user uses `pheno_threshold_percent` or `pheno_threshold_abs` parameters.  

`pheno_threshold_percent` : float, optional *(The default is None)*  
Accepts values between (0-100). If any particular phenotype is below the user defined threshold, it is re-categorized as `unknown`. Generally used to deal with low background false positives.  

`pheno_threshold_abs` : int, optional *(The default is None)*  
Serves the same purpose as that of `pheno_threshold_percent`. However, an absolute number can be passed. For example, if user passes in 10- any phenotype that contains less than 10 cells will be re-categorized as `unknown`.  

**Returns**

`AnnData` object with cell phenotypes added as a new column in `adata.obs`

**Example**

```
phenotype = pd.read_csv('path/to/csv/file/')
adata = sm.tl.phenotype_cells (adata, phenotype=phenotype, label="phenotype")
```
