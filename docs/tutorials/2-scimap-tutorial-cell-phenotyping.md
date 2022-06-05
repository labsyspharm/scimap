# Cell-Phenotyping using Scimap


```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 18:10:06 2020
@author: Ajit Johnson Nirmal
Scimap Cell Phenotyping Tutorial
"""
```




    '\nCreated on Fri Jun 28 18:10:06 2020\n@author: Ajit Johnson Nirmal\nScimap Cell Phenotyping Tutorial\n'



## Tutorial material

You can download the material for this tutorial from the following [link:](https://www.dropbox.com/s/rra13zir52o9hio/getting_started%20and%20phenotyping.zip?dl=0)  
The presentation files are available [here:](https://github.com/ajitjohnson/Jupyter-Notebooks/blob/master/tutorials/scimap_tutorial/getting_started%20and%20phenotyping/scimap_tutorial.pdf)

## Tutorial video


```python
from IPython.display import HTML
HTML('<iframe width="450" height="250" src="https://www.youtube.com/embed/knh5elRksUk" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>')
```




<iframe width="450" height="250" src="https://www.youtube.com/embed/knh5elRksUk" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>




```python
# Load necessary libraries
import sys
import os
import anndata as ad
import pandas as pd
import scanpy as sc
import seaborn as sns; sns.set(color_codes=True)

# Import Scimap
import scimap as sm
```


```python
# Set the working directory
os.chdir ("/Users/aj/Desktop/scimap_tutorial/")
```


```python
# Load data
adata = ad.read('tutorial_data.h5ad')
```

## Clustering and data exploration

You could use clustering and marker expression analysis within clusters to assign cell types similar to what is carried out with single-cell sequencing data.


```python
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=10) # Computing the neighborhood graph
```


```python
sc.tl.umap(adata) # Build a UMAP to visualize the neighbourhood graph
```


```python
sc.pl.umap(adata, color=['CD3D', 'CD20', 'CD163'], cmap= 'vlag', use_raw=False, s=30) # Plot the UMAP
```


![png](scimap-tutorial-cell-phenotyping_files/scimap-tutorial-cell-phenotyping_13_0.png)


We can already begin to spot issues with carrying out this mode of phenotyping approach. As you can see there is an area of co-expression of CD3D and CD20, which is likely because of segmentation errors. Additionally the boundaries are not distinct between cell-types and it is highly likely that errors will be introduced due to this reason. 


```python
sc.tl.leiden(adata, resolution = 1) # Clustering the neighborhood graph
```


```python
sc.pl.umap(adata, color=['leiden', 'CD3D', 'CD20'],cmap= 'vlag', use_raw=False) # View the clustering
```


![png](scimap-tutorial-cell-phenotyping_files/scimap-tutorial-cell-phenotyping_16_0.png)


### Finding marker genes


```python
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, fontsize=16)
```


![png](scimap-tutorial-cell-phenotyping_files/scimap-tutorial-cell-phenotyping_18_0.png)


From the above plots, it is likely that clusters 1, 2 and 7 could be combined to form a T cell cluster. However, as mentioned earlier the boundaries are not clear and it only get increasingly complex as one would want to perform deeper phenotyping such as CD4 helper T cells, CD8 T cells, regulatory T cells and so on. 

Additionally, marker analsyis suggests that CD30 is most expressed in cluster 8. If you look at the actual image, you will realize that CD30 is not expressed by any cell in this image and the analysis is picking up high background fluorescence. 

## Probability distribution based phenotyping

This approach is more labor intensive, however is significantly more sensitive and much more scalable than clustering based approaches. Takes less than 5 mins to run over a million cells once the gates are identified.

#### In order to run the method, you need 2 things 
- a gating workflow strategy `.csv file`
- manual gates `.csv file`. If manual gates are not provided, the algorithm will attempt to rescale the data by fitting two gaussians on the data. However, it is adviced to perform manual gating as I have found it to be more sensitive.

**The algorithm involves three steps:**
1. Identify the gates using `sm.pl.gate_finder`
2. Rescale the data based on the identified gates using `sm.pp.rescale`
3. Run the phenotyping algorithm on the rescaled data using `sm.tl.phenotype`

#### Define manual gates to rescale data before running the phenotyping algorithm
Instantiating the Qt GUI can take a few seconds and if you create the Viewer before it is finished, the kernel will die and the viewer will not launch. For this reason the %gui qt magic command should always be run in a separate cell from creating the viewer


```python
%gui qt
```

### Step 1: Identify the gates using `sm.pl.gate_finder`


```python
image_path = '/Users/aj/Desktop/scimap_tutorial/reactive_core.tif'
marker_of_interest = 'CD45'
```


```python
sm.pl.gate_finder (image_path, adata, marker_of_interest, 
                   from_gate = 5, to_gate = 9, increment = 0.1, 
                   markers=['ASMA','DNA11', 'CD20', 'CD3D'], point_size=6)
```

### Step 2: Rescale the data based on the identified gates using `sm.pp.rescale`

Note: Below we are passing a `manual_gates.csv` into the `gate` parameter.
This contatins gates that were visually determined using the `sm.pl.gate_finder`
function. For the markers included in the `manual_gates.csv` file, 
the function will scale the data such that cells with expression greater than the gate 
will be considered as positive for that marker and cells
with expression below the gate is considered negative. <br>

For markers that are not included in the `manual_gates.csv` file, the function
will automatically try to determine a gate by running a gaussian mixture model
algorithm on the data. <br>

### Note (for >=v.0.22.0)
Please note that passing manual gates for multiple images has been introduced in `scimap >=v.0.22.0`

```python
# Load the manual gates and rescale the data based on the gates
manual_gate = pd.read_csv('manual_gates.csv')
adata = sm.pp.rescale (adata, gate=manual_gate)
```

    Scaling Image [mcmicro_output]
    Categories (1, object): [mcmicro_output]
    Finding the optimal gate for CD10
    Finding the optimal gate for CD2
    Finding the optimal gate for CD30
    Finding the optimal gate for CD43
    Finding the optimal gate for CD5
    Finding the optimal gate for CD57
    Finding the optimal gate for CD7
    Finding the optimal gate for KI67
    Finding the optimal gate for MHCI
    Finding the optimal gate for PDL1
    Finding the optimal gate for PS6
    Finding the optimal gate for PSTAT3
    Scaling ASMA
    Scaling CD163
    Scaling CD206
    Scaling CD68
    Scaling CD20
    Scaling CD21
    Scaling CD3D
    Scaling CD45
    Scaling CD56
    Scaling CD8A
    Scaling FOXP3
    Scaling CD11B
    Scaling CD11C
    Scaling CD15
    Scaling CD4
    Scaling PD1
    Scaling HLADR
    Scaling CD25



```python
# View the scaled data (note that the log data is replaced with scaled data)
# If you ever want the log data back you will need to run-  np.log1p(adata.raw.X)
adata.X 
```




    array([[0.17841106, 0.45723783, 0.49234127, ..., 0.15973007, 0.1665647 ,
            0.20024123],
           [0.155838  , 0.21377199, 0.34924023, ..., 0.1522421 , 0.0885678 ,
            0.15338667],
           [0.29090098, 0.37456273, 0.50752728, ..., 0.21580293, 0.17788475,
            0.18778977],
           ...,
           [0.33621626, 0.70161347, 0.70359459, ..., 0.26182348, 0.60810172,
            0.22068843],
           [0.15324935, 0.51044092, 0.60877724, ..., 0.11845023, 0.26558105,
            0.08391592],
           [0.18923565, 0.431478  , 0.58019522, ..., 0.09423545, 0.24047052,
            0.12733514]])



### Step 3: Run the phenotyping algorithm on the rescaled data using `sm.tl.phenotype`


```python
# Load the gating workflow
phenotype = pd.read_csv('phenotype_workflow.csv')
adata = sm.tl.phenotype_cells (adata, phenotype=phenotype, label="phenotype") 
```

    Phenotyping Other Immune cells
    Phenotyping ASMA+ cells
    -- Subsetting Other Immune cells
    Phenotyping T cells
    Phenotyping B cells
    Phenotyping Myeloid Lineage
    Phenotyping NK cells
    Phenotyping Granulocytes
    -- Subsetting Myeloid Lineage
    Phenotyping T cells
    Phenotyping B cells
    Phenotyping NK cells
    Phenotyping Granulocytes
    Phenotyping CD68+ Macrophages
    Phenotyping M2 Macrophages
    Phenotyping Myeloid Dendritic cells
    Phenotyping Follicular Dendritic cells
    -- Subsetting T cells
    Phenotyping CD4 T cells
    Phenotyping CD8 T cells
    -- Subsetting CD4 T cells
    Phenotyping Regulatory T cells
    Phenotyping Follicular Helper T cells
    -- Subsetting CD8 T cells
    Phenotyping PD1+ T cells
    -- Subsetting Myeloid Dendritic cells
    Phenotyping CD25+ Dendritic cells
    Consolidating the phenotypes across all groups



```python
# Summary of the phenotyping
adata.obs['phenotype'].value_counts()
```




    B cells                       2037
    CD4 T cells                    502
    ASMA+ cells                    420
    Regulatory T cells             418
    CD8 T cells                    322
    Follicular Helper T cells      282
    T cells                        146
    Unknown                        140
    Other Immune cells             137
    Myeloid Dendritic cells        124
    Follicular Dendritic cells      87
    Myeloid Lineage                 80
    PD1+ T cells                    63
    M2 Macrophages                  55
    Granulocytes                     8
    NK cells                         3
    CD25+ Dendritic cells            1
    Name: phenotype, dtype: int64



It is likely that `CD25+ Dendritic cells, NK cells & Granulocytes` are artifacts. You could set `pheno_threshold_abs= 10` to move these cells into `unknown` category.

**Once the phenotyping is performed, it is adviced to overlay the phenotypes on the image and check if they are correct. If not, alter the `phenotyping workflow` file or the `manual gate` to account for the errors.**


```python
# View phenotypes
sm.pl.image_viewer (image_path, adata, overlay = 'phenotype', point_color='white', point_size=6)
```


```python
# View Leiden clustering
sm.pl.image_viewer (image_path, adata, overlay = 'leiden', point_color='white', point_size=6)
```

#### Heatmap and UMAP of the probability based phenotyping


```python
sc.tl.dendrogram(adata, groupby='phenotype')
```


```python
sc.pl.matrixplot(adata, var_names= adata.var.index, groupby='phenotype', dendrogram=True, use_raw=False, cmap="vlag", standard_scale='var')
```


![png](scimap-tutorial-cell-phenotyping_files/scimap-tutorial-cell-phenotyping_40_0.png)





    GridSpec(2, 3, height_ratios=[0, 10.5], width_ratios=[9.6, 0.8, 0.2])




```python
sc.pl.umap(adata, color=['leiden', 'phenotype']) # View the clustering
```


![png](scimap-tutorial-cell-phenotyping_files/scimap-tutorial-cell-phenotyping_41_0.png)



```python
sns.set(rc={'figure.figsize':(11.7,8.27)})
sc.pl.umap(adata, color=['phenotype'],legend_loc='on data', title='', frameon=False, s = 100) # View the clustering
```


![png](scimap-tutorial-cell-phenotyping_files/scimap-tutorial-cell-phenotyping_42_0.png)



```python
sc.pl.umap(adata, color=['CD3D', 'PD1', 'CD20'],cmap= 'vlag', use_raw=False, frameon=False, s = 100) # View the clustering
```


![png](scimap-tutorial-cell-phenotyping_files/scimap-tutorial-cell-phenotyping_43_0.png)


As it can be seen from above 3 UMAP's it would have been very difficult to find the Follicular helper T cells by a pure clustering approach. Also, the B cells as can be seen above does not from a nice seperate cluster. These among other illustrate the importance of using the probability based algorithm for deep phenotyping.


```python
# Confirm Follicular helper T cells in the image
sm.pl.image_viewer (image_path, adata, 
                    overlay = 'phenotype', overlay_category=['Follicular Helper T cells'], 
                    markers = ['CD3D','CD20','PD1','CD8A','CD4','DNA11'],
                    point_color='white', point_size=6)
```


```python
# Save the results
adata.write('tutorial_data.h5ad')
```

**This concludes this tutorial**
