#!/usr/bin/env python3
"""
Infer gene regulatory network using GENIE3
"""

import pandas as pd
import numpy as np
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
import pickle

def main():
    # Get parameters from snakemake
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]
    n_jobs = snakemake.params.n_jobs
    seed = snakemake.params.seed
    
    print(f"Loading expression matrix from {input_file}")
    expression_df = pd.read_csv(input_file, index_col=0)
    
    print(f"Running GENIE3 with {n_jobs} cores...")
    
    # Import GENIE3 from arboreto
    from arboreto.algo import genie3
    
    # Run GENIE3
    adjacencies = genie3(
        expression_data=expression_df,
        gene_names=expression_df.index.tolist(),
        regulators='all',
        client_or_address=None,
        seed=seed
    )
    
    print(f"GENIE3 inference complete. Found {len(adjacencies)} gene-gene interactions.")
    print(f"Saving adjacencies to {output_file}")
    
    # Save adjacencies
    adjacencies.to_csv(output_file, sep='\t', index=False)
    
    print("Gene regulatory network inference complete!")

if __name__ == "__main__":
    main()
