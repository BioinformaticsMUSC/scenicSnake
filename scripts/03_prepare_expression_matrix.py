#!/usr/bin/env python3
"""
Prepare expression matrix for SCENIC analysis
"""

import scanpy as sc
import pandas as pd
import numpy as np
import loompy as lp

def main():
    # Get parameters from snakemake
    input_file = snakemake.input.h5ad_file
    output_file = snakemake.output[0]
    sample_id = snakemake.params.sample_id
    
    # Load single h5ad file
    print(f"Loading h5ad file for sample {sample_id}: {input_file}")
    adata = sc.read_h5ad(input_file)
    print(f"Sample {sample_id} data shape: {adata.shape}")

    # Use snakemake.params for loom_preparation parameters
    loom_params = snakemake.config['loom_preparation']
    celltype_column = loom_params.get('celltype_column', None)
    
    # Create loom file
    row_attrs = {
        "Gene": np.array(adata.var_names),
    }
    col_attrs = {
        "CellID": np.array(adata.obs_names),
        "nGene": np.array(np.sum(adata.X.transpose() > 0, axis=0)).flatten(),
        "nUMI": np.array(np.sum(adata.X.transpose(), axis=0)).flatten(),
    }
    if celltype_column and celltype_column in adata.obs.columns:
        col_attrs["CellAnno"] = np.array(adata.obs[celltype_column])
    
    print(f"Creating loom file: {output_file}")
    lp.create(
        output_file,
        adata.X.transpose(),
        row_attrs=row_attrs,
        col_attrs=col_attrs,
    )
    print(f"Expression matrix preparation complete for sample {sample_id}!")

if __name__ == "__main__":
    main()
