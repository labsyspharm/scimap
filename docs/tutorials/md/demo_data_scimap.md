# 📥 SCIMAP Demo Data

A short guide to downloading our demo data for trying out the **scimap toolbox**.

**You can download the data in two methods:**

1. Utilize the built-in download function of scimap.
2. Visit the website through your browser and manually download the data.

### 1. Built in function


```python
import scimap as sm
```

    Running SCIMAP  1.3.16



```python
# Provide a download directory
download_directory = '/Users/aj/Downloads'
sm.hl.downloadDemoData (download_directory)
```

    Downloading scimapExampleData.zip...


    scimapExampleData.zip: 100%|███████████████████████████████████████████████████████████| 287M/287M [04:35<00:00, 1.09MiB/s]

    Downloaded scimapExampleData.zip to /Users/aj/Downloads/scimapExampleData.zip


    


A zipped data folder will appear in the specified directory. Simply unzip it, and you'll be all set to start.

### 2. Download from Zonodo

- Go to [Zonodo](https://zenodo.org/records/10845625)
- Click on the download button
- unzip the downloaded file

DOI: [10.5281/zenodo.10845624](https://zenodo.org/records/10845625)

If you have any questions or encounter any issues while following the tutorial, please don't hesitate to contact us for assistance.
