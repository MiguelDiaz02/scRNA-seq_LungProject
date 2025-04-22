#!/bin/bash

# Exit on error
set -e

# Set input and output directories
INPUT_DIR="raw_files"
OUTPUT_DIR="cleaned_h5_files"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run CellBender on all *_raw_gene_bc_matrices_*.h5 files using conda run
echo "Running CellBender on raw input files..."
for infile in "$INPUT_DIR"/*_raw_gene_bc_matrices_*.h5; do
    filename=$(basename "$infile")
    outfile="${OUTPUT_DIR}/cb_${filename}"
    echo "Processing $filename..."
    conda run -n cellbender cellbender remove-background \
        --input "$infile" \
        --output "$outfile" \
        --total-droplets-included 10000 \
        --cuda
done

echo "All files processed. Outputs saved in $OUTPUT_DIR"
