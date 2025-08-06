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
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]

    adata = sc.read_h5ad(input_file)

    # Use snakemake.params for loom_preparation parameters
    loom_params = snakemake.params['loom_preparation']

    if loom_params['split']:
        split_condition = loom_params['split_condition']
        split_values = loom_params['split_values']
        celltype_column = loom_params.get('celltype_column', None)
        assert split_condition in adata.obs.columns
        for split_value in split_values:
            assert split_value in adata.obs[split_condition].values
            tmp_adata = adata[adata.obs[split_condition] == split_value].copy()
            row_attrs = {
                "Gene": np.array(tmp_adata.var_names),
            }
            col_attrs = {
                "CellID": np.array(tmp_adata.obs_names),
                "nGene": np.array(np.sum(tmp_adata.X.transpose() > 0, axis=0)).flatten(),
                "nUMI": np.array(np.sum(tmp_adata.X.transpose(), axis=0)).flatten(),
            }
            if celltype_column:
                col_attrs["CellAnno"] = np.array(tmp_adata.obs[celltype_column])
            save_name = f"{output_file}_{split_value}.loom"
            lp.create(
                save_name,
                tmp_adata.X.transpose(),
                row_attrs=row_attrs,
                col_attrs=col_attrs,
            )
    else:
        celltype_column = loom_params.get('celltype_column', None)
        row_attrs = {
            "Gene": np.array(adata.var_names),
        }
        col_attrs = {
            "CellID": np.array(adata.obs_names),
            "nGene": np.array(np.sum(adata.X.transpose() > 0, axis=0)).flatten(),
            "nUMI": np.array(np.sum(adata.X.transpose(), axis=0)).flatten(),
        }
        if celltype_column:
            col_attrs["CellAnno"] = np.array(adata.obs[celltype_column])
        save_name = f"{output_file}"
        lp.create(
            save_name,
            adata.X.transpose(),
            row_attrs=row_attrs,
            col_attrs=col_attrs,
        )
    print("Expression matrix preparation complete!")

if __name__ == "__main__":
    main()
