---
title: SCIMAP
description: Spatial Single-Cell Analysis Toolkit
hide:
  - toc        # Hide table of contents
---

# Home

[![build: Unix-Mac-Win](https://github.com/ajitjohnson/scimap/actions/workflows/build-unix-mac-win.yml/badge.svg)](https://github.com/ajitjohnson/scimap/actions/workflows/build-unix-mac-win.yml)
[![Docs](https://github.com/ajitjohnson/scimap/actions/workflows/docs.yml/badge.svg)](https://github.com/ajitjohnson/scimap/actions/workflows/docs.yml)
[![Downloads](https://pepy.tech/badge/scimap)](https://pepy.tech/project/scimap)
[![PyPI Version](https://img.shields.io/pypi/v/scimap.svg)](https://pypi.org/project/scimap)
[![PyPI License](https://img.shields.io/pypi/l/scimap.svg)](https://pypi.org/project/scimap)
[![Gitter chat](https://badges.gitter.im/scimap_io/community.png)](https://gitter.im/scimap_io/community)
<br>
<br>
<img src="./assets/scimap_logo.jpg" style="" >
<br>

Scimap is a scalable toolkit for analyzing spatial molecular data. The underlying framework is generalizable to spatial datasets mapped to XY coordinates. The package uses the [anndata](https://anndata.readthedocs.io/en/stable/anndata.AnnData.html) framework making it easy to integrate with other popular single-cell analysis toolkits. It includes preprocessing, phenotyping, visualization, clustering, spatial analysis and differential spatial testing. The Python-based implementation efficiently deals with large datasets of millions of cells.
<br>

#### Notice for Apple M1 users
Please note that multiple python packages have not yet extended support for M1 users. 
Below is a temporary solution to install scimap in `Apple M1` machines. 
Please follow the instructions in the given order.

```
# create and load a new environment
conda create -y -n scimap -c andfoy python=3.9 pyqt
conda activate scimap

# if you do not have xcode please install it
xcode-select --install

# if you do not have homebrew please install it
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# if you do not have cmake install it
brew install cmake

# install h5py
brew install hdf5@1.12
HDF5_DIR=/opt/homebrew/Cellar/hdf5/ pip install --no-build-isolation h5py

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


*SCIMAP* development is led by [Ajit Johnson Nirmal](https://ajitjohnson.com/) at the Laboratory of Systems Pharmacology, Harvard Medical School.

### Funding
This work is supported by the following NIH grant K99-CA256497
