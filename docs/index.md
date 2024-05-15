---
title: SCIMAP
description: Spatial Single-Cell Analysis Toolkit
hide:
  - toc        # Hide table of contents
  - navigation
---

# Home

[![build: Unix-Mac-Win](https://github.com/ajitjohnson/scimap/actions/workflows/build-unix-mac-win.yml/badge.svg)](https://github.com/ajitjohnson/scimap/actions/workflows/build-unix-mac-win.yml)
[![Docs](https://github.com/ajitjohnson/scimap/actions/workflows/docs.yml/badge.svg)](https://github.com/ajitjohnson/scimap/actions/workflows/docs.yml)
[![Downloads](https://pepy.tech/badge/scimap)](https://pepy.tech/project/scimap)
[![PyPI Version](https://img.shields.io/pypi/v/scimap.svg)](https://pypi.org/project/scimap)
[![PyPI License](https://img.shields.io/pypi/l/scimap.svg)](https://pypi.org/project/scimap)
<br>
<br>
<img src="./assets/scimap_logo.jpg" style="" >
<br>

Scimap is a scalable toolkit for analyzing spatial molecular data. The underlying framework is generalizable to spatial datasets mapped to XY coordinates. The package uses the [anndata](https://anndata.readthedocs.io/en/stable/anndata.AnnData.html) framework making it easy to integrate with other popular single-cell analysis toolkits. It includes preprocessing, phenotyping, visualization, clustering, spatial analysis and differential spatial testing. The Python-based implementation efficiently deals with large datasets of millions of cells.
  
SCIMAP operates on segmented single-cell data derived from imaging data using tools such as cellpose or MCMICRO. The essential inputs for SCIMAP are: (a) a single-cell expression matrix and (b) the X and Y coordinates for each cell. Additionally, multi-stack OME-TIFF or TIFF images can be optionally provided to enable visualization of the data analysis on the original raw images.
<br>


*SCIMAP* development was led by [Ajit Johnson Nirmal](https://ajitjohnson.com/), Harvard Medical School.  
Check out other tools from the [Nirmal Lab](https://nirmallab.com/tools/). 


### Funding
This work was supported by the following NIH grant K99-CA256497


