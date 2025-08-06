#!/usr/bin/env python3
"""
Load and filter single-cell RNA-seq data
"""

import scanpy as sc
import pandas as pd
import anndata as ann

# Set up scanpy
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

def main():
    # Get parameters from snakemake
    input_file = snakemake.input.data
    output_file = snakemake.output.filtered
    stats_file = snakemake.output.stats
    
    min_genes = snakemake.params.min_genes
    min_cells = snakemake.params.min_cells
    max_genes = snakemake.params.max_genes
    max_mito = snakemake.params.max_mito
    
    # Load data
    print(f"Loading data from {input_file}")
    adata = sc.read(input_file)
    
    # Store initial statistics
    n_cells_initial = adata.n_obs
    n_genes_initial = adata.n_vars
    
    # Calculate mitochondrial gene percentage
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    
    # Filter cells and genes
    print("Applying filters...")
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    # Filter cells with too many genes (potential doublets)
    adata = adata[adata.obs.n_genes_by_counts < max_genes, :]
    
    # Filter cells with high mitochondrial gene expression
    adata = adata[adata.obs.pct_counts_mt < max_mito, :]
    
    # Store final statistics
    n_cells_final = adata.n_obs
    n_genes_final = adata.n_vars
    
    # Save filtered data
    print(f"Saving filtered data to {output_file}")
    adata.write(output_file)
    
    # Save statistics
    stats = {
        'initial_cells': n_cells_initial,
        'initial_genes': n_genes_initial,
        'final_cells': n_cells_final,
        'final_genes': n_genes_final,
        'cells_removed': n_cells_initial - n_cells_final,
        'genes_removed': n_genes_initial - n_genes_final
    }
    
    with open(stats_file, 'w') as f:
        for key, value in stats.items():
            f.write(f"{key}: {value}\n")
    
    print(f"Filtering complete. Cells: {n_cells_initial} -> {n_cells_final}")
    print(f"Genes: {n_genes_initial} -> {n_genes_final}")

if __name__ == "__main__":
    main()
