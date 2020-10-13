#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 17:03:56 2020
@author: Ajit Johnson Nirmal
Function to sub-cluster a cluster of interest. Particularly useful to check if 
there are sub-phenotypes after performing the gating based phenotyping.
"""

# Import library
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import parc
from sklearn.cluster import MiniBatchKMeans
