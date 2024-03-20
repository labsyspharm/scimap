# ⬇️ Export data from SCIMAP into csv

After completing some analysis and when you're ready to export the data, you can utilize the specified function for this purpose. However, it's important to note that if you plan to return and continue your analysis later, you should rely on the `.h5ad` file you've saved. While exporting to a CSV file can be handy, only a subset of the data—specifically, what's contained in `adata.obs` and `adata.X`—is exported. Be aware that all other compartments of the data are not preserved during the CSV export process.


```python
# import packages
import scimap as sm
import anndata as ad
```

    Running SCIMAP  1.3.14



```python
# Load the data that we saved in the last tutorial (with ROIs added)
adata = ad.read_h5ad('/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/resources/exemplarData/scimapExampleData/scimapExampleData.h5ad')
```


```python
# by default the raw data is exported
output_dir= '/Users/aj/Dropbox (Partners HealthCare)/nirmal lab/resources/exemplarData/scimapExampleData'
sm.hl.scimap_to_csv(adata, output_dir=output_dir, file_name='scimapProcessed.csv')
```

Refer to the documentation to learn how to export additional layers beyond just the raw data.


```python

```
