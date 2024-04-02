---
title: Get Started
description: Embark on Your Single-Cell Analysis Journey with Scimap
hide:
  - toc        # Hide table of contents
  - navigation
---

# Welcome to Scimap! 🚀

Begin your journey into the fascinating world of spatial single-cell analysis with Scimap, a comprehensive toolkit designed to empower your single-cell data exploration.

## Installation 📦

Kick off by installing Scimap through these simple commands:  
We highly advise creating a separate environment for installing Scimap to ensure a smooth and conflict-free setup. For comprehensive guidance on this process, please refer to our tutorials.

```
# If you have conda installed
conda create --name scimap python=3.10
conda activate scimap
```

Install `scimap` directly into an activated virtual environment:

```python
pip install scimap

# install scimap with napari for visualization
# Please note that napari requires a GUI toolkit like PyQt. If you encounter any errors due to your operating system, you should install scimap and napari separately, following their own documentation, and make sure to add them into the same virtual environment, instead of using this command.
pip install scimap[napari]

>>> import scimap as sm
```

This setup provides you with the foundational tools needed for your single-cell analysis endeavors.

## Data Loading and Integration 🔄

Scimap champions the interoperability of single-cell analysis tools by embracing the `AnnData` data structure. This strategic choice allows seamless use of numerous single-cell analysis utilities alongside `AnnData`.

### The AnnData Framework 🧬

An `AnnData` object, `adata`, encapsulates a data matrix `adata.X`, with annotations of observations `adata.obs` and variables `adata.var` as `pd.DataFrame`, and unstructured annotation `adata.uns` as a dictionary. Observation and variable names are accessible via `adata.obs_names` and `adata.var_names`, respectively. AnnData objects support slicing, similar to dataframes: `adata_subset = adata[:, list_of_gene_names]`. Explore `AnnData` further in the [official documentation](https://anndata.readthedocs.io/en/stable/anndata.AnnData.html#anndata.AnnData).

### Initializing AnnData Objects 🔄

To begin with an *AnnData object*, proceed as follows:

```python
import anndata as ad
import pandas as pd

# Load your data
data = pd.read_csv('counts_matrix.csv')  # Your single-cell counts matrix
meta = pd.read_csv('meta_data.csv')      # Associated MetaData

# Craft the AnnData object
adata = ad.AnnData(data)
adata.obs = meta
```

!!! Note 📝  
    Leveraging the [mcmicro](https://github.com/labsyspharm/mcmicro-nf) pipeline? Scimap simplifies converting `mcmicro` outputs into an `AnnData` object:

```python
filepath = ['/path/to/file.csv']
adata = sm.pp.mcmicro_to_scimap(filepath)
```

## Navigating Non-mcmicro Data 🧐

What to do if your dataset wasn't processed with mcmicro? Ensuring your data aligns with Scimap's expectations is vital for smooth analysis:

- **Spatial Assumptions:** Spatial functions expect XY coordinates in 'X_centroid' and 'Y_centroid' columns. If your data differs, specify your columns when using these functions.
- **Manual Data Integration:** Here's how to manually prepare your data for Scimap analysis:

```python
# Import necessary packages
import anndata as ad
import pandas as pd

# Load data
data = pd.read_csv('path/to/counts_table.csv')  # Counts matrix
meta = pd.read_csv('path/to/meta_data.csv')     # Meta data with coordinates

# Merge data and metadata for the AnnData object
adata = ad.AnnData(data)
adata.obs = meta

# preserve raw data
adata.raw = adata

# log transform data
adata = sm.pp.log1p(adata)

# Add marker annotation for visualization with Napari
adata.uns['all_markers'] = ['list', 'of', 'all', 'markers', 'in', 'image']

```

### Key Steps for Data Preparation 🗝️

1. **Unique Image Identification:** Include a `imageid` column in your metadata for easy data retrieval.
2. **Preserving Raw Data:** Store unprocessed data in `adata.raw` for reference.
3. **Log Transformation Layer:** Create a `log` layer for log-transformed data normalization.
4. **Marker Annotation:** Keep a record of image markers in `adata.uns['all_markers']` for clarity during analysis.

### Saving Your AnnData Object 💾

An AnnData object centralizes your data and analysis, making it simple to share and collaborate. To save your work:

```python
# Save your AnnData object
adata.write('/path/to/scimapExampleData.h5ad')
```

This streamlined approach facilitates comprehensive analyses, enabling you to leverage Scimap's full suite of tools and integrate with other analysis frameworks seamlessly.

## Your Workflow Journey 🛤️

With Scimap, navigate through a suite of tools designed to enrich your analysis:

- **Pre-Processing Tools:** `sm.pp.<tool>` for data preparation.
- **Analysis Tools:** `sm.tp.<tool>` for in-depth insights.
- **Plotting Tools:** `sm.pl.<tool>` for impactful visualizations.
- **Helper Tools:** `sm.hl.<tool>` for additional functionalities.

Embark on your spatial single-cell analysis journey with Scimap today and unlock the potential within your data. 🌟
