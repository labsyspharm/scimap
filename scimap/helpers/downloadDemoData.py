#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 15:04:25 2024
@author: aj
"""

import requests
import os

def downloadDemoData(directory):

    # Ensure the target directory exists
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"Created directory: {directory}")

    # Get record details from Zenodo API
    api_url = f"https://zenodo.org/api/records/10845625"
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

        # Download the file
        print(f"Downloading {filename}...")
        with requests.get(file_url, stream=True) as file_response:
            with open(save_path, 'wb') as file:
                for chunk in file_response.iter_content(chunk_size=8192):
                    file.write(chunk)
        print(f"Downloaded {filename} to {save_path}")