#!/bin/bash

# Define relative destination directory
DEST_DIR="./raw_files/"

# Create the directory if it doesn't exist
mkdir -p "$DEST_DIR"

# GEO dataset download URL
FILE_URL="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE122960&format=file"

# Name of the file to save
FILE_NAME="GSE122960_RAW.tar"

# Download the file to the destination directory
echo "Downloading $FILE_NAME..."
wget -O $DEST_DIR$FILE_NAME $FILE_URL

# Extract the downloaded tar file
echo "Extracting $FILE_NAME..."
tar -xvf $DEST_DIR$FILE_NAME -C $DEST_DIR

echo "Download and extraction complete."
