---
title: Contributing Guide
description: Spatial Single-Cell Analysis Toolkit
hide:
  - navigation
---

# Contributing guide

We are thrilled you're interested in contributing to our tool! To ensure a smooth integration of your contributions, we have outlined a simple process. Contributions can be in the form of new functionalities, bug fixes, or enhancements to existing features. Please see the provided steps below and never hesitate to contact us. Follow these steps to get started:

If you are a new user, we recommend checking out the detailed Docs.  
  
## Setting up a development installation
In order to make changes to napari, you will need to [fork](https://github.com/labsyspharm/scimap) the repository. If you are not familiar with git, we recommend reading up on this [guide](https://docs.github.com/en/get-started/using-git/about-git#basic-git-commands).
  
Before we set up `SCIMAP`, we highly recommend using a environment manager like Conda. Using an environment manager like Conda allows you to create and manage isolated environments with specific package versions and dependencies.

```
# use the terminal (mac/linux) and anaconda promt (windows) to run the following
conda create --name scimap -y python=3.9
conda activate scimap
```

We use [poetry](https://python-poetry.org/docs/) to manage scimap dependencies. Intall poetry and scimap into the environment.

```
# install poetry within the environment
pipx install poetry

# install scimap within the environment
poetry install

```

## Set up the contributions

We invite contributions aimed at enhancing the performance, functionality, and reliability of the existing codebase within `scimap`. Additionally, we are open to the integration of innovative tools into scimap to facilitate the seamless analysis of multiplexed imaging data.

If you are interested in contributing a new tool to `scimap`, please encapsulate all requisite functions within a single Python script. This script should be comprehensively designed to ensure full functionality of the proposed tool. Once prepared, place this script in the designated directory path: `scimap/scimap/external`.

Your function must adhere to the structure outlined below.


```

# Required libraries
import [library name] as [alias]
...
...

# Your function
def functionName (adata,
                    ...
                    ...
                    ... # other necessary parameters
                    ...
                    verbose=True,
                    outputDir=None):
    
    # CODE BLOCK: Function implementation
    
    
    # OUTPUT
    # not needed if the output is a plot

    if outputDir:
        adata.write(outputDir / name the file)
    else:    
        return adata

```

If your function requires dependencies not present in `scimap`, ascertain this by examining the `pyproject.toml` file. In such cases, it is advisable to package your tool independently. Consequently, `scimap` could function merely as an API interface to your package, enabling it to be easily installed via pip when a user wishes to utilize your tool. This strategy is primarily aimed at reducing the maintenance effort associated with managing numerous dependencies.

## Add Documentation

All contributions must be documented directly within the code to maintain clarity and usability. This includes a short description at the top of your script, followed by a detailed comment block that explains the functionality of your code, parameters, return values, and provides an example usage.

Your code contributions should start with a script header followed by an abstract in a docstring, giving a brief overview of its functionality. Below is a template based on the provided script:

```
#!/usr/bin/env python3
# Created on [Date] by [Your Name]
"""


!!! abstract "Short Description"
    [Briefly describe the function's purpose and its utility in the tool. Explain how it contributes to the tool's functionality and any specific features it offers.]
    
    Results are stored in the [specify location, e.g., `.uns` section of the Anndata object] for easy access and further analysis.



## Function
"""

# Required libraries
import [library name] as [alias]

# Your function
def functionName (parameters):

    """
Parameters:
    parameter1 (type): 
        Description of parameter1.
    parameter2 (type, optional): 
        Description of parameter2. 
    
Returns:
    return_type (type): 
        Description of what is returned.
    
Example:

    ```python
    
    # Example usage of your function
    adata = sm.ex.yourFunction(parameter1, parameter2, ...)
    
    ```
    """
    
    # Function implementation

```


## Add Unit Tests

Contributing unit tests is as vital as contributing code! Well-crafted unit tests help maintain the tool's integrity and facilitate future enhancements. Here's how you can add unit tests following our specified format:

1. Unit Test Structure:  
Unit tests for our tool must be added to the test_external.py file, adhering to a specific template for consistency and effectiveness. Ensure your tests are comprehensive, covering various scenarios and edge cases.
  
2. Test Template and Example:  
Below is a template for creating a unit test, including the setup for test data and an example test case. This template uses pytest, a powerful testing framework for Python. Ensure you have pytest installed before proceeding.

```
import pytest
import sys, os
import anndata as ad

@pytest.fixture
def adata():
    # Adjust the path to where the example data is located in your environment
    image_path = os.getcwd() + '/scimap/tests/data/example_data.h5ad'
    adata = ad.read(image_path)
    return adata

# Example unit test for your function
def test_yourFunction(adata):

    from scimap.external.yourFunction import yourFunction
    
    # Add your testing logic here
    # Example: Assert the expected outcome from running your function on `adata`

```

Please ensure your test cases are well-documented, explaining the purpose of the test and the expected outcomes. This will help other contributors understand and maintain the test cases over time.

## Add tutorial (optional)

Develop a Jupyter notebook to demonstrate the capabilities of your tool, utilizing the example dataset available in /scimap/tests/data/ whenever feasible, even if the data does not precisely reflect biological accuracy. Please store the completed notebook in the `scimap/docs/tutorials/nb` directory.































































