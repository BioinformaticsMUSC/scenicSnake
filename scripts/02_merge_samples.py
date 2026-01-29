#!/usr/bin/env python3
"""
Merge metadata from multiple h5ad files for SCENIC analysis
"""

import scanpy as sc
import pandas as pd
import numpy as np

def main():
    # Get parameters from snakemake
    input_files = snakemake.input.h5ad_files
    output_file = snakemake.output[0]
    
    # Load and concatenate multiple h5ad files
    print(f"Loading {len(input_files)} h5ad files for metadata merging...")
    adata_list = []
    for i, input_file in enumerate(input_files):
        print(f"Loading file {i+1}/{len(input_files)}: {input_file}")
        adata = sc.read_h5ad(input_file)
        # Add sample information to obs
        adata.obs['sample_id'] = f"sample_{i+1}"
        adata.obs['source_file'] = input_file
        adata_list.append(adata)
    
    # Concatenate all samples
    print("Concatenating samples...")
    merged_adata = sc.concat(adata_list, axis=0, join="outer", fill_value=0)
    print(f"Merged data shape: {merged_adata.shape}")
    
    # Save merged data
    print(f"Saving merged metadata to: {output_file}")
    merged_adata.write_h5ad(output_file)
    print("Metadata merge complete!")

if __name__ == "__main__":
    main()