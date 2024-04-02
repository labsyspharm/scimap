# Single-Cell Image Analysis Package
<br>

[![build: Unix-Mac-Win](https://github.com/ajitjohnson/scimap/actions/workflows/build-unix-mac-win.yml/badge.svg)](https://github.com/ajitjohnson/scimap/actions/workflows/build-unix-mac-win.yml)
[![Docs](https://github.com/ajitjohnson/scimap/actions/workflows/docs.yml/badge.svg)](https://github.com/ajitjohnson/scimap/actions/workflows/docs.yml)
[![Downloads](https://pepy.tech/badge/scimap)](https://pepy.tech/project/scimap)
[![PyPI Version](https://img.shields.io/pypi/v/scimap.svg)](https://pypi.org/project/scimap)
[![PyPI License](https://img.shields.io/pypi/l/scimap.svg)](https://pypi.org/project/scimap)
[![DOI](https://zenodo.org/badge/271099296.svg)](https://zenodo.org/badge/latestdoi/271099296)

<br>

<img src="./docs/assets/scimap_logo.jpg" style="max-width:700px;width:100%" >

<br> 

Scimap is a scalable toolkit for analyzing spatial molecular data. The underlying framework is generalizable to spatial datasets mapped to XY coordinates. The package uses the [anndata](https://anndata.readthedocs.io/en/stable/anndata.AnnData.html) framework making it easy to integrate with other popular single-cell analysis toolkits. It includes preprocessing, phenotyping, visualization, clustering, spatial analysis and differential spatial testing. The Python-based implementation efficiently deals with large datasets of millions of cells.

## Installation

We strongly recommend installing `scimap` in a fresh virtual environment.

```
# If you have conda installed
conda create --name scimap python=3.10
conda activate scimap
```

Install `scimap` directly into an activated virtual environment:

```python
pip install scimap

# Install scimap with napari for visualization
# Please note that napari requires a GUI toolkit like PyQt. If you encounter any errors due to your operating system, you should install scimap and napari separately, following their own documentation, and make sure to add them into the same virtual environment, instead of using this command.
pip install scimap[napari]
```

After installation, the package can be imported as:

```python
$ python
>>> import scimap as sm
```


## Get Started

#### Detailed documentation of `scimap` functions and tutorials are available [here](http://scimap.xyz/).

*SCIMAP* development was led by [Ajit Johnson Nirmal](https://ajitjohnson.com/), Harvard Medical School.  
Check out other tools from the [Nirmal Lab](https://nirmallab.com/tools/). 


## Contibute
Interested in contributing to the package? Check out our guidelines at [https://scimap.xyz/contribute/](https://scimap.xyz/contribute/) for detailed instructions.


## Funding
This work was supported by the following NIH grant K99-CA256497

