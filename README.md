# Single-Cell Image Analysis Package
<br>

[![build: Unix-Mac-Win](https://github.com/ajitjohnson/scimap/actions/workflows/build-unix-mac-win.yml/badge.svg)](https://github.com/ajitjohnson/scimap/actions/workflows/build-unix-mac-win.yml)
[![Docs](https://github.com/ajitjohnson/scimap/actions/workflows/docs.yml/badge.svg)](https://github.com/ajitjohnson/scimap/actions/workflows/docs.yml)
[![Downloads](https://pepy.tech/badge/scimap)](https://pepy.tech/project/scimap)
[![PyPI Version](https://img.shields.io/pypi/v/scimap.svg)](https://pypi.org/project/scimap)
[![PyPI License](https://img.shields.io/pypi/l/scimap.svg)](https://pypi.org/project/scimap)
[![Gitter chat](https://badges.gitter.im/scimap_io/community.png)](https://gitter.im/scimap_io/community)

<br>

<img src="./docs/assets/scimap_logo.png" style="max-width:700px;width:100%" >

<br> 

Scimap is a scalable toolkit for analyzing spatial molecular data. The underlying framework is generalizable to spatial datasets mapped to XY coordinates. The package uses the [anndata](https://anndata.readthedocs.io/en/stable/anndata.AnnData.html) framework making it easy to integrate with other popular single-cell analysis toolkits. It includes preprocessing, phenotyping, visualization, clustering, spatial analysis and differential spatial testing. The Python-based implementation efficiently deals with large datasets of millions of cells.

## Installation

We strongly recommend installing `scimap` in a fresh virtual environment.

```
# If you have conda installed
conda create --name scimap python=3.7
conda activate scimap
```

Install `scimap` directly into an activated virtual environment:

```python
$ pip install scimap
```

After installation, the package can be imported as:

```python
$ python
>>> import scimap as sm
```

## Get Started


#### Detailed documentation of `scimap` functions and tutorials are available [here](http://scimap.xyz/).

*SCIMAP* development is led by [Ajit Johnson Nirmal](https://ajitjohnson.com/) at the Laboratory of Systems Pharmacology, Harvard Medical School.

## Funding
This work is supported by the following NIH grant K99-CA256497

