#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 15:04:25 2024
@author: aj
"""
      
import requests
from tqdm import tqdm
import os

def downloadDemoData(directory, api_url=None):
    """
    Downloads all files from a Zenodo record into the specified directory,
    showing a progress bar for each file.

    Parameters:
    - directory: The directory where the files will be saved.

    Returns:
    - None
    """

    # Ensure the target directory exists
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"Created directory: {directory}")

    # Get record details from Zenodo API
    if api_url is None:
        api_url = "https://zenodo.org/api/records/10845625"
    
    response = requests.get(api_url)
    if response.status_code != 200:
        print(f"Failed to retrieve record details. HTTP status code: {response.status_code}")
        return

    # Extract file links from the record information
    record_data = response.json()
    files = record_data.get('files', [])

    # Download each file
    for file_info in files:
        file_url = file_info['links']['self']
        filename = file_info['key']
        save_path = os.path.join(directory, filename)

        # Download the file with progress bar
        print(f"Downloading {filename}...")
        response = requests.get(file_url, stream=True)
        total_size_in_bytes = int(response.headers.get('content-length', 0))
        block_size = 1024  # 1 Kibibyte

        with open(save_path, 'wb') as file, tqdm(
                desc=filename,
                total=total_size_in_bytes,
                unit='iB',
                unit_scale=True,
                unit_divisor=1024,
            ) as progress_bar:
            for data in response.iter_content(block_size):
                progress_bar.update(len(data))
                file.write(data)

        print(f"Downloaded {filename} to {save_path}")
