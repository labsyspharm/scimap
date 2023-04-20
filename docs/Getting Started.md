---
title: Get Started
description: Spatial Single-Cell Analysis Toolkit
hide:
  - toc        # Hide table of contents
  - navigation
---

# Getting Started

Import Scimap as:

```python
pip install scimap
import scimap as sm
```

### Load Data

In order to make data analysis tools interoperable, `scimap` has adopted the the `AnnData` data structure. This allows users to use a wealth of single-cell analysis tools that works with AnnData structuring.

At the most basic level, an `AnnData` object `adata` stores a data matrix `adata.X`, annotation of observations `adata.obs` and variables `adata.var` as `pd.DataFrame` and unstructured annotation `adata.uns` as dict. Names of observations and variables can be accessed via `adata.obs_names` and `adata.var_names`, respectively. AnnData objects can be sliced like dataframes, for example, `adata_subset = adata[:, list_of_gene_names]`. For more, see the `AnnData` [page](https://anndata.readthedocs.io/en/stable/anndata.AnnData.html#anndata.AnnData).

To initialize an *AnnData object*, the following can be performed.

```python
import anndata as ad
import pandas as pd

# Load the data
data = pd.read_csv('counts_matrix.csv') # Single-Cell counts matrix
meta = pd.read_csv('meta_data.csv') # MetaData

# Create the AnnData object
adata = ad.AnnData(data)
adata.obs = meta
```

!!! note
    If you used [mcmicro](https://github.com/labsyspharm/mcmicro-nf) pipeline to process your images, `scimap` provides a handy function to convert `mcmicro` output to `AnnData` object.


``` python
filepath = ['/path/to/file.csv']
adata = sm.pp.mcmicro_to_scimap (filepath)

```

### Work Flow

The typical workflow then consists of subsequent calls of `scimap` tools:

> Pre-Processing Tools under `sm.pp.<tool>`  
> Analysis Tools under `sm.tl.<tool>`  
> Plotting Tools under `sm.pl.<tool>`  
> Helper Tools under `sm.hl.<tool>`  
