#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 10:18:22 2023
@author: aj
"""

# print out scimap version
import argparse
import toml

# Read the version number from pyproject.toml
with open('../../pyproject.toml', 'r') as f:
    pyproject = toml.load(f)
    version = pyproject['tool']['poetry']['version']

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-V', '--version', action='store_true', help='Print version number')
args = parser.parse_args()

# Print the version number if requested
if args.version:
    print(version)