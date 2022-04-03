# Single-Cell Image Analysis Package
<br>

[![build: Unix-Mac-Win](https://github.com/ajitjohnson/scimap/actions/workflows/build-unix-mac-win.yml/badge.svg)](https://github.com/ajitjohnson/scimap/actions/workflows/build-unix-mac-win.yml)
[![Docs](https://github.com/ajitjohnson/scimap/actions/workflows/docs.yml/badge.svg)](https://github.com/ajitjohnson/scimap/actions/workflows/docs.yml)
[![Downloads](https://pepy.tech/badge/scimap)](https://pepy.tech/project/scimap)
[![PyPI Version](https://img.shields.io/pypi/v/scimap.svg)](https://pypi.org/project/scimap)
[![PyPI License](https://img.shields.io/pypi/l/scimap.svg)](https://pypi.org/project/scimap)
[![Gitter chat](https://badges.gitter.im/scimap_io/community.png)](https://gitter.im/scimap_io/community)
[![DOI](https://zenodo.org/badge/271099296.svg)](https://zenodo.org/badge/latestdoi/271099296)

<br>

<img src="./docs/assets/scimap_logo.jpg" style="max-width:700px;width:100%" >

<br> 

Scimap is a scalable toolkit for analyzing spatial molecular data. The underlying framework is generalizable to spatial datasets mapped to XY coordinates. The package uses the [anndata](https://anndata.readthedocs.io/en/stable/anndata.AnnData.html) framework making it easy to integrate with other popular single-cell analysis toolkits. It includes preprocessing, phenotyping, visualization, clustering, spatial analysis and differential spatial testing. The Python-based implementation efficiently deals with large datasets of millions of cells.

## Installation

We strongly recommend installing `scimap` in a fresh virtual environment.

```
# If you have conda installed
conda create --name scimap python=3.8
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

### Notice for Apple M1 users
Please note that multiple python packages have not yet extended support for M1 users. 
Below is a solution to install scimap in `Apple M1` machines

```
# reate and load a new environment
conda create --name scimap python=3.8 -y
conda activate scimap

# if you do not have xcode please install it
xcode-select --install

# if you do not have homebrew please install it
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# if you do not have cmake install it
brew install cmake

# install h5py
brew install hdf5@1.12
export HDF5_DIR=/opt/homebrew/Cellar/hdf5/1.12.1_1/
pip install --no-binary=h5py h5py

# install llvmlite
conda install llvmlite -y

# install leidenalg
pip install git+https://github.com/vtraag/leidenalg.git

# install scimap
pip install -U scimap

# uninstall 
conda remove llvmlite -y
pip uninstall numba -y
pip uninstall numpy -y

# reinstall this specific version of llvmlite (ignore errors/warning)
pip install -i https://pypi.anaconda.org/numba/label/wheels_experimental_m1/simple llvmlite

# reinstall this specific version of numpy (ignore errors/warning)
pip install numpy==1.22.3

# reinstall this specific version of numba (ignore errors/warning)
pip install -i https://pypi.anaconda.org/numba/label/wheels_experimental_m1/simple numba

```

## Get Started


#### Detailed documentation of `scimap` functions and tutorials are available [here](http://scimap.xyz/).

*SCIMAP* development is led by [Ajit Johnson Nirmal](https://ajitjohnson.com/) at the Laboratory of Systems Pharmacology, Harvard Medical School.

## Funding
This work is supported by the following NIH grant K99-CA256497

