# Single-Cell Image Analysis Package
<br>

[![build: Unix-Mac-Win](https://github.com/ajitjohnson/scimap/actions/workflows/build-unix-mac-win.yml/badge.svg)](https://github.com/ajitjohnson/scimap/actions/workflows/build-unix-mac-win.yml)
[![Docs](https://github.com/ajitjohnson/scimap/actions/workflows/docs.yml/badge.svg)](https://github.com/ajitjohnson/scimap/actions/workflows/docs.yml)
[![Downloads](https://pepy.tech/badge/scimap)](https://pepy.tech/project/scimap)
[![PyPI Version](https://img.shields.io/pypi/v/scimap.svg)](https://pypi.org/project/scimap)
[![PyPI License](https://img.shields.io/pypi/l/scimap.svg)](https://pypi.org/project/scimap)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.06604/status.svg)](https://doi.org/10.21105/joss.06604)

<br>

<img src="./docs/assets/scimap_logo.jpg" style="max-width:700px;width:100%" >

<br> 

Scimap is a scalable toolkit for analyzing spatial molecular data. The underlying framework is generalizable to spatial datasets mapped to XY coordinates. The package uses the [anndata](https://anndata.readthedocs.io/en/stable/anndata.AnnData.html) framework making it easy to integrate with other popular single-cell analysis toolkits. It includes preprocessing, phenotyping, visualization, clustering, spatial analysis and differential spatial testing. The Python-based implementation efficiently deals with large datasets of millions of cells.

## Citing SCIMAP
Nirmal et al., (2024). SCIMAP: A Python Toolkit for Integrated Spatial Analysis of Multiplexed Imaging Data. Journal of Open Source Software, 9(97), 6604, [https://doi.org/10.21105/joss.06604](https://joss.theoj.org/papers/10.21105/joss.06604#)

## Installation

We strongly recommend installing `scimap` in a fresh virtual environment.

```
# If you have conda installed
conda create --name scimap python=3.10
conda activate scimap
```

Install `scimap` directly into an activated virtual environment:
  
**Firstly, we suggest installing `scimap` and `napari` together to enable visualization out of the box. Keep in mind, `napari` needs a GUI toolkit, such as PyQt. If you run into any issues because of your computer's operating system, install `scimap` and `napari` separately by following the guidance in `napari's` documentation.**

Here's how you can install both using pip:

```python
pip install "scimap[napari]"
```

**If you encounter a problem with PyQt6 during the installation, you can install `scimap` alone first. Later on, if you find you need `napari`, you can go ahead and install it by itself.**

To install just `scimap`:

```python
pip install scimap
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

