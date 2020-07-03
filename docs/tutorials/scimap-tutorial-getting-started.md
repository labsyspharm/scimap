# Getting Started with Scimap


```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 23:11:32 2020
@author: Ajit Johnson Nirmal
Scimap Getting Started tutorial
"""
```




    '\nCreated on Fri Jun 26 23:11:32 2020\n@author: Ajit Johnson Nirmal\nScimap Getting Started tutorial\n'




```python
# Before you start make sure you have installed the following packages
# pip install scimap
# pip install scanpy
# pip install leidenalg
# pip install PyQt5
```

## Tutorial material

You can download the material for this tutorial from the following [link:](https://www.dropbox.com/s/rra13zir52o9hio/getting_started%20and%20phenotyping.zip?dl=0)  
The presentation files are available [here:](https://github.com/ajitjohnson/Jupyter-Notebooks/blob/master/tutorials/scimap_tutorial/getting_started%20and%20phenotyping/scimap_tutorial.pdf)

## Tutorial video


```python
from IPython.display import HTML
HTML('<iframe width="871" height="490" src="https://www.youtube.com/embed/knh5elRksUk" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>')
```




<iframe width="871" height="490" src="https://www.youtube.com/embed/knh5elRksUk" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>




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

## Load data using AnnData


```python
# Load data
data = pd.read_csv ('counts_table.csv') # Counts matrix
meta = pd.read_csv ('meta_data.csv') # Meta data like x and y coordinates 

# combine the data and metadata file to generate the AnnData object
adata = ad.AnnData (data)
adata.obs = meta
```

Print adata to check for it's content


```python
adata
```




    AnnData object with n_obs × n_vars = 4825 × 48
        obs: 'X_centroid', 'Y_centroid', 'Area', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', 'Solidity', 'Extent', 'Orientation'




```python
adata.obs # prints the meta data
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>X_centroid</th>
      <th>Y_centroid</th>
      <th>Area</th>
      <th>MajorAxisLength</th>
      <th>MinorAxisLength</th>
      <th>Eccentricity</th>
      <th>Solidity</th>
      <th>Extent</th>
      <th>Orientation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>511.555556</td>
      <td>9.846154</td>
      <td>117</td>
      <td>14.532270</td>
      <td>10.273628</td>
      <td>0.707261</td>
      <td>0.959016</td>
      <td>0.750000</td>
      <td>-0.695369</td>
    </tr>
    <tr>
      <th>1</th>
      <td>579.330097</td>
      <td>9.398058</td>
      <td>103</td>
      <td>16.056286</td>
      <td>8.776323</td>
      <td>0.837396</td>
      <td>0.903509</td>
      <td>0.613095</td>
      <td>1.115707</td>
    </tr>
    <tr>
      <th>2</th>
      <td>630.958333</td>
      <td>12.883333</td>
      <td>120</td>
      <td>15.222005</td>
      <td>10.310756</td>
      <td>0.735653</td>
      <td>0.975610</td>
      <td>0.681818</td>
      <td>0.151616</td>
    </tr>
    <tr>
      <th>3</th>
      <td>745.194631</td>
      <td>16.275168</td>
      <td>149</td>
      <td>14.380200</td>
      <td>13.404759</td>
      <td>0.362027</td>
      <td>0.967532</td>
      <td>0.662222</td>
      <td>-0.270451</td>
    </tr>
    <tr>
      <th>4</th>
      <td>657.173653</td>
      <td>18.035928</td>
      <td>167</td>
      <td>17.675831</td>
      <td>12.110106</td>
      <td>0.728428</td>
      <td>0.943503</td>
      <td>0.695833</td>
      <td>-0.810890</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>4820</th>
      <td>559.597403</td>
      <td>1091.577922</td>
      <td>154</td>
      <td>18.150307</td>
      <td>11.683288</td>
      <td>0.765281</td>
      <td>0.900585</td>
      <td>0.570370</td>
      <td>-0.342315</td>
    </tr>
    <tr>
      <th>4821</th>
      <td>619.983871</td>
      <td>1092.959677</td>
      <td>248</td>
      <td>21.734414</td>
      <td>15.565820</td>
      <td>0.697912</td>
      <td>0.864111</td>
      <td>0.551111</td>
      <td>1.432242</td>
    </tr>
    <tr>
      <th>4822</th>
      <td>583.317073</td>
      <td>1093.573171</td>
      <td>82</td>
      <td>12.060039</td>
      <td>9.539789</td>
      <td>0.611784</td>
      <td>0.964706</td>
      <td>0.630769</td>
      <td>0.203023</td>
    </tr>
    <tr>
      <th>4823</th>
      <td>607.064394</td>
      <td>1101.583333</td>
      <td>264</td>
      <td>22.549494</td>
      <td>15.905321</td>
      <td>0.708858</td>
      <td>0.882943</td>
      <td>0.661654</td>
      <td>0.691838</td>
    </tr>
    <tr>
      <th>4824</th>
      <td>641.592486</td>
      <td>1100.132948</td>
      <td>346</td>
      <td>23.149806</td>
      <td>19.375564</td>
      <td>0.547257</td>
      <td>0.945355</td>
      <td>0.791762</td>
      <td>-1.390516</td>
    </tr>
  </tbody>
</table>
<p>4825 rows × 9 columns</p>
</div>




```python
adata.X # prints the counts table
```




    array([[16640.564  ,   719.6325 ,   527.7094 , ...,  1085.735  ,
              218.54701,  3170.47   ],
           [16938.3    ,   686.5534 ,   469.30096, ...,  1075.6407 ,
              164.48544,  3116.767  ],
           [16243.542  ,   819.4167 ,   604.39166, ...,  1164.3917 ,
              227.74167,  3156.1084 ],
           ...,
           [28656.256  ,   878.2561 ,   585.3293 , ...,  1233.183  ,
             1243.5488 ,  3194.195  ],
           [22054.818  ,   685.8485 ,   424.85226, ...,  1031.2424 ,
              313.32574,  3038.8105 ],
           [23992.854  ,   850.25146,   529.89886, ...,  1000.5578 ,
              285.98267,  3087.3005 ]], dtype=float32)




```python
adata.var[0:5] # prints the first 5 channel or marker names
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>DNA1</th>
    </tr>
    <tr>
      <th>BG1</th>
    </tr>
    <tr>
      <th>BG2</th>
    </tr>
    <tr>
      <th>BG3</th>
    </tr>
    <tr>
      <th>DNA2</th>
    </tr>
  </tbody>
</table>
</div>



You would have noticed that
- the data is not in log scale
- All the DNA channels are there
- The background channels are there
If we diretly perform clustering or any other type of analysis, the above mentioned factors may affect the results and so it is recommended to remove them.

## Load data using scimap's helper function

Use this if the single-cell data was generated using **mcmicro pipeline**. With this function though many of the above limitations can be imediately addressed. By default it removes DNA channels and you can pass any channel name into `drop_markers` parameter inorder to not import them.


```python
image_path = ['/Users/aj/Desktop/scimap_tutorial/mcmicro_output.csv']
adata = sm.pp.mcmicro_to_scimap (image_path, drop_markers = ["PERK", "NOS2","BG1","BG2","BG3","ACTIN"])
```

    Loading mcmicro_output.csv


Check adata contents now as we did previously


```python
adata
```




    AnnData object with n_obs × n_vars = 4825 × 30
        obs: 'X_centroid', 'Y_centroid', 'Area', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', 'Solidity', 'Extent', 'Orientation', 'imageid'
        uns: 'all_markers'




```python
adata.X # Will now contain log normalized data
```




    array([[6.3674684, 6.4287267, 7.3826084, ..., 6.990933 , 5.3915663,
            8.061951 ],
           [6.340171 , 6.094227 , 7.339796 , ..., 6.981601 , 5.1088834,
            8.044872 ],
           [6.503502 , 6.3549495, 7.4734573, ..., 7.0608125, 5.4325933,
            8.057412 ],
           ...,
           [6.5583014, 6.660794 , 7.4199724, ..., 7.1181645, 7.1265283,
            8.069404 ],
           [6.3370404, 6.281594 , 7.2397914, ..., 6.939489 , 5.7504296,
            8.01955  ],
           [6.3805585, 6.180567 , 7.2547846, ..., 6.909312 , 5.659422 ,
            8.035377 ]], dtype=float32)




```python
adata.raw.X # contains the raw data
```




    array([[ 581.5812 ,  618.38464, 1606.7778 , ..., 1085.735  ,  218.54701,
            3170.47   ],
           [ 565.8932 ,  442.29126, 1539.3981 , ..., 1075.6407 ,  164.48544,
            3116.767  ],
           [ 666.475  ,  574.3333 , 1759.6833 , ..., 1164.3917 ,  227.74167,
            3156.1084 ],
           ...,
           [ 704.0732 ,  780.1707 , 1667.9878 , ..., 1233.183  , 1243.5488 ,
            3194.195  ],
           [ 564.1212 ,  533.64014, 1392.803  , ..., 1031.2424 ,  313.32574,
            3038.8105 ],
           [ 589.2572 ,  482.2659 , 1413.8584 , ..., 1000.5578 ,  285.98267,
            3087.3005 ]], dtype=float32)




```python
adata.obs # prints the meta data
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>X_centroid</th>
      <th>Y_centroid</th>
      <th>Area</th>
      <th>MajorAxisLength</th>
      <th>MinorAxisLength</th>
      <th>Eccentricity</th>
      <th>Solidity</th>
      <th>Extent</th>
      <th>Orientation</th>
      <th>imageid</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>mcmicro_output_1</th>
      <td>511.555556</td>
      <td>9.846154</td>
      <td>117</td>
      <td>14.532270</td>
      <td>10.273628</td>
      <td>0.707261</td>
      <td>0.959016</td>
      <td>0.750000</td>
      <td>-0.695369</td>
      <td>mcmicro_output</td>
    </tr>
    <tr>
      <th>mcmicro_output_2</th>
      <td>579.330097</td>
      <td>9.398058</td>
      <td>103</td>
      <td>16.056286</td>
      <td>8.776323</td>
      <td>0.837396</td>
      <td>0.903509</td>
      <td>0.613095</td>
      <td>1.115707</td>
      <td>mcmicro_output</td>
    </tr>
    <tr>
      <th>mcmicro_output_3</th>
      <td>630.958333</td>
      <td>12.883333</td>
      <td>120</td>
      <td>15.222005</td>
      <td>10.310756</td>
      <td>0.735653</td>
      <td>0.975610</td>
      <td>0.681818</td>
      <td>0.151616</td>
      <td>mcmicro_output</td>
    </tr>
    <tr>
      <th>mcmicro_output_4</th>
      <td>745.194631</td>
      <td>16.275168</td>
      <td>149</td>
      <td>14.380200</td>
      <td>13.404759</td>
      <td>0.362027</td>
      <td>0.967532</td>
      <td>0.662222</td>
      <td>-0.270451</td>
      <td>mcmicro_output</td>
    </tr>
    <tr>
      <th>mcmicro_output_5</th>
      <td>657.173653</td>
      <td>18.035928</td>
      <td>167</td>
      <td>17.675831</td>
      <td>12.110106</td>
      <td>0.728428</td>
      <td>0.943503</td>
      <td>0.695833</td>
      <td>-0.810890</td>
      <td>mcmicro_output</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>mcmicro_output_4821</th>
      <td>559.597403</td>
      <td>1091.577922</td>
      <td>154</td>
      <td>18.150307</td>
      <td>11.683288</td>
      <td>0.765281</td>
      <td>0.900585</td>
      <td>0.570370</td>
      <td>-0.342315</td>
      <td>mcmicro_output</td>
    </tr>
    <tr>
      <th>mcmicro_output_4822</th>
      <td>619.983871</td>
      <td>1092.959677</td>
      <td>248</td>
      <td>21.734414</td>
      <td>15.565820</td>
      <td>0.697912</td>
      <td>0.864111</td>
      <td>0.551111</td>
      <td>1.432242</td>
      <td>mcmicro_output</td>
    </tr>
    <tr>
      <th>mcmicro_output_4823</th>
      <td>583.317073</td>
      <td>1093.573171</td>
      <td>82</td>
      <td>12.060039</td>
      <td>9.539789</td>
      <td>0.611784</td>
      <td>0.964706</td>
      <td>0.630769</td>
      <td>0.203023</td>
      <td>mcmicro_output</td>
    </tr>
    <tr>
      <th>mcmicro_output_4824</th>
      <td>607.064394</td>
      <td>1101.583333</td>
      <td>264</td>
      <td>22.549494</td>
      <td>15.905321</td>
      <td>0.708858</td>
      <td>0.882943</td>
      <td>0.661654</td>
      <td>0.691838</td>
      <td>mcmicro_output</td>
    </tr>
    <tr>
      <th>mcmicro_output_4825</th>
      <td>641.592486</td>
      <td>1100.132948</td>
      <td>346</td>
      <td>23.149806</td>
      <td>19.375564</td>
      <td>0.547257</td>
      <td>0.945355</td>
      <td>0.791762</td>
      <td>-1.390516</td>
      <td>mcmicro_output</td>
    </tr>
  </tbody>
</table>
<p>4825 rows × 10 columns</p>
</div>



### We can use scanpy package to explore the data


```python
sc.pl.highest_expr_genes(adata, n_top=20, ) # Most expressing proteins
```


![png](scimap-tutorial-getting-started_files/scimap-tutorial-getting-started_26_0.png)



```python
sc.tl.pca(adata, svd_solver='arpack') # peform PCA
sc.pl.pca(adata, color='KI67') # scatter plot in the PCA coordinates
```


![png](scimap-tutorial-getting-started_files/scimap-tutorial-getting-started_27_0.png)



```python
sc.pl.pca_variance_ratio(adata) # PCs to the total variance in the data
```


![png](scimap-tutorial-getting-started_files/scimap-tutorial-getting-started_28_0.png)



```python
# Save the results
adata.write('tutorial_data.h5ad')
```

**This concludes the `getting started` tutorial, continue with the `phenotyping` tutorial.**
