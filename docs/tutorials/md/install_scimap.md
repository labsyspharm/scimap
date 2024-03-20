# 📁 Setting up SCIMAP

Before we set up *SCIMAP*, we highly recommend using an environment manager like Conda. Using an environment manager like Conda allows you to create and manage isolated environments with specific package versions and dependencies.

Download and Install the right [conda](https://docs.anaconda.com/free/miniconda/) based on the opertating system that you are using

<hr>

## Let's create a new conda environment and install SCIMAP

Use the terminal (mac/linux) and anaconda promt (windows) to run the following command

```
conda create --name scimap -y python=3.10
```

Install SCIMAP within the conda environment.

```
conda activate scimap
pip install scimap
```

<hr>

## Set up Jupyter Notebook / Spyder or any interactive IDE

Install jupyter notebook within the conda environment


```
pip install notebook ipywidgets
```

After installation, open Jupyter Notebook by typing the following command in the terminal, ensuring that the cspot environment is activated and you are within the environment before executing the jupyter notebook command.

```
jupyter notebook
```

We will talk about how to run *SCIMAP* in the next tutorial.
