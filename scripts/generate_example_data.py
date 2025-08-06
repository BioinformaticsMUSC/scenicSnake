#!/usr/bin/env python3
"""
Generate example single-cell data for testing the SCENIC workflow
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import argparse
import os

def generate_example_data(n_cells=1000, n_genes=2000, output_path="data/example_data.h5ad"):
    """
    Generate synthetic single-cell RNA-seq data for testing
    
    Parameters:
    - n_cells: Number of cells to generate
    - n_genes: Number of genes to generate  
    - output_path: Path to save the generated data
    """
    
    print(f"Generating example data with {n_cells} cells and {n_genes} genes...")
    
    # Generate random gene names
    gene_names = [f"GENE_{i:04d}" for i in range(n_genes)]
    
    # Add some common gene patterns
    gene_names[0:10] = ['MT-CO1', 'MT-CO2', 'MT-ATP6', 'MT-ATP8', 'MT-ND1', 
                        'MT-ND2', 'MT-ND3', 'MT-ND4', 'MT-ND5', 'MT-CYB']  # Mitochondrial genes
    gene_names[10:20] = ['ACTB', 'GAPDH', 'RPL13A', 'RPLP0', 'B2M',
                        'HPRT1', 'TBP', 'GUSB', 'HMBS', 'YWHAZ']  # Housekeeping genes
    
    # Generate cell barcodes
    cell_names = [f"CELL_{i:04d}" for i in range(n_cells)]
    
    # Generate expression matrix with realistic properties
    np.random.seed(42)
    
    # Create a realistic expression pattern
    # Most genes have low expression, few genes have high expression
    base_expression = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes))
    
    # Add some cell type specific patterns
    n_cell_types = 3
    cells_per_type = n_cells // n_cell_types
    cell_types = []
    
    for i in range(n_cell_types):
        start_idx = i * cells_per_type
        end_idx = (i + 1) * cells_per_type if i < n_cell_types - 1 else n_cells
        cell_types.extend([f"CellType_{i+1}"] * (end_idx - start_idx))
        
        # Add cell type specific expression
        type_specific_genes = range(i * 50, (i + 1) * 50)  # 50 genes per type
        base_expression[start_idx:end_idx, type_specific_genes] += np.random.negative_binomial(10, 0.2, 
                                                                                              size=(end_idx - start_idx, 50))
    
    # Convert to sparse matrix (more realistic for scRNA-seq)
    X = sparse.csr_matrix(base_expression.astype(np.float32))
    
    # Create AnnData object
    adata = sc.AnnData(X=X)
    adata.obs_names = cell_names
    adata.var_names = gene_names
    
    # Add metadata
    adata.obs['cell_type'] = cell_types
    adata.obs['condition'] = np.random.choice(['control', 'treatment'], size=n_cells)
    adata.obs['batch'] = np.random.choice(['batch1', 'batch2'], size=n_cells)
    
    # Add gene metadata
    adata.var['gene_type'] = 'protein_coding'
    adata.var.loc[gene_names[0:10], 'gene_type'] = 'mitochondrial'
    adata.var.loc[gene_names[10:20], 'gene_type'] = 'housekeeping'
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Save the data
    print(f"Saving example data to {output_path}")
    adata.write(output_path)
    
    # Also create a summary
    summary_path = output_path.replace('.h5ad', '_summary.txt')
    with open(summary_path, 'w') as f:
        f.write(f"Example Single-Cell Dataset Summary\n")
        f.write(f"===================================\n\n")
        f.write(f"Cells: {n_cells}\n")
        f.write(f"Genes: {n_genes}\n")
        f.write(f"Cell types: {', '.join(set(cell_types))}\n")
        f.write(f"Conditions: control, treatment\n")
        f.write(f"Batches: batch1, batch2\n")
        f.write(f"\nExpression matrix statistics:\n")
        f.write(f"Mean expression: {X.mean():.2f}\n")
        f.write(f"Sparsity: {(1 - X.nnz / (X.shape[0] * X.shape[1])) * 100:.1f}%\n")
    
    print(f"Example data generated successfully!")
    print(f"Summary saved to {summary_path}")
    
    return adata

def main():
    parser = argparse.ArgumentParser(description='Generate example single-cell data for SCENIC workflow')
    parser.add_argument('--cells', type=int, default=1000, help='Number of cells (default: 1000)')
    parser.add_argument('--genes', type=int, default=2000, help='Number of genes (default: 2000)')
    parser.add_argument('--output', type=str, default='data/example_data.h5ad', help='Output file path')
    
    args = parser.parse_args()
    
    generate_example_data(
        n_cells=args.cells,
        n_genes=args.genes,
        output_path=args.output
    )

if __name__ == "__main__":
    main()
