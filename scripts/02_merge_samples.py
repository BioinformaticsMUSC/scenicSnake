#!/usr/bin/env python3
"""
Merge multiple filtered samples into a single dataset
"""

import scanpy as sc
import pandas as pd
import numpy as np

def main():
    # Get parameters from snakemake
    input_files = snakemake.input
    output_file = snakemake.output.merged
    stats_file = snakemake.output.stats
    
    print(f"Merging {len(input_files)} samples...")
    
    # Load all samples
    adatas = []
    sample_stats = []
    
    for i, file_path in enumerate(input_files):
        print(f"Loading sample {i+1}: {file_path}")
        adata = sc.read(file_path)
        
        # Add sample information
        sample_name = file_path.split('/')[-1].replace('_filtered.h5ad', '')
        adata.obs['sample'] = sample_name
        adata.obs['batch'] = f"batch_{i+1}"
        
        adatas.append(adata)
        
        # Collect stats
        sample_stats.append({
            'sample': sample_name,
            'n_cells': adata.n_obs,
            'n_genes': adata.n_vars
        })
    
    # Concatenate all samples
    print("Concatenating samples...")
    merged_adata = sc.concat(adatas, join='outer', index_unique='-')
    
    # Remove genes that are not expressed in any cell after merging
    sc.pp.filter_genes(merged_adata, min_cells=1)
    
    print(f"Merged dataset: {merged_adata.n_obs} cells, {merged_adata.n_vars} genes")
    
    # Save merged data
    merged_adata.write(output_file)
    
    # Save statistics
    stats_df = pd.DataFrame(sample_stats)
    total_stats = {
        'total_samples': len(input_files),
        'total_cells': merged_adata.n_obs,
        'total_genes': merged_adata.n_vars,
        'samples_breakdown': stats_df.to_dict('records')
    }
    
    with open(stats_file, 'w') as f:
        f.write(f"Merged Dataset Statistics\n")
        f.write(f"========================\n\n")
        f.write(f"Total samples: {stats_df.shape[0]}\n")
        f.write(f"Total cells: {merged_adata.n_obs}\n")
        f.write(f"Total genes: {merged_adata.n_vars}\n\n")
        f.write("Per-sample breakdown:\n")
        for _, row in stats_df.iterrows():
            f.write(f"  {row['sample']}: {row['n_cells']} cells, {row['n_genes']} genes\n")
    
    print("Sample merging complete!")

if __name__ == "__main__":
    main()
