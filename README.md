# Single-Cell Image Analysis Package

Scimap is a scalable toolkit for analyzing single-cell multiplex imaging data. The package uses the [anndata](https://anndata.readthedocs.io/en/stable/anndata.AnnData.html) framework making it easy to integrate with other popular single-cell analysis toolkits such as [scanpy](https://scanpy.readthedocs.io/en/latest/#). It includes preprocessing, phenotyping, visualization, clustering, spatial analysis and differential spatial testing. The Python-based implementation efficiently deals with datasets of more than one million cells.


[![Unix Build Status](https://img.shields.io/travis/ajitjohnson/scimap/master.svg?label=unix)](https://travis-ci.org/ajitjohnson/scimap)
[![Windows Build Status](https://img.shields.io/appveyor/ci/ajitjohnson/scimap/master.svg?label=windows)](https://ci.appveyor.com/project/ajitjohnson/scimap)
[![Documentation Status](https://readthedocs.org/projects/scimap-doc/badge/?version=latest)](https://scimap-doc.readthedocs.io/en/latest/?badge=latest)
[![Downloads](https://pepy.tech/badge/scimap)](https://pepy.tech/project/scimap)
[![PyPI Version](https://img.shields.io/pypi/v/scimap.svg)](https://pypi.org/project/scimap)
[![PyPI License](https://img.shields.io/pypi/l/scimap.svg)](https://pypi.org/project/scimap)
<!--[![Scrutinizer Code Quality](https://img.shields.io/scrutinizer/g/ajitjohnson/scimap.svg)](https://scrutinizer-ci.com/g/ajitjohnson/scimap/?branch=master)-->
<!--[![Coverage Status](https://img.shields.io/coveralls/ajitjohnson/scimap/master.svg)](https://coveralls.io/r/ajitjohnson/scimap) -->

## Installation

Install it directly into an activated virtual environment:

```python
$ pip install scimap
```

After installation, the package can imported:

```python
$ python
>>> import scimap as sm
```

### To get started, read the [documentation](https://scimap-doc.readthedocs.io/en/latest/).
