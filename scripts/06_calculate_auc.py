#!/usr/bin/env python3
"""
Calculate AUC (Area Under the Curve) scores for regulon activity
"""

import pandas as pd
import numpy as np
import json
from pyscenic.aucell import aucell

def main():
    # Get parameters from snakemake
    regulons_file = snakemake.input.regulons
    expression_file = snakemake.input.expression
    auc_file = snakemake.output.auc_matrix
    binary_file = snakemake.output.binary_activity
    
    auc_threshold = snakemake.params.auc_threshold
    
    print("Calculating regulon activity scores (AUC)...")
    
    # Load regulons
    print(f"Loading regulons from {regulons_file}")
    with open(regulons_file, 'r') as f:
        regulons_dict = json.load(f)
    
    # Convert to pySCENIC regulon format
    regulons = []
    for regulon_name, regulon_data in regulons_dict.items():
        # Create a simple regulon object
        regulon = SimpleRegulon(
            name=regulon_name,
            targets=regulon_data['targets']
        )
        regulons.append(regulon)
    
    print(f"Loaded {len(regulons)} regulons")
    
    # Load expression matrix
    print(f"Loading expression matrix from {expression_file}")
    expression_df = pd.read_csv(expression_file, index_col=0)
    
    print(f"Expression matrix shape: {expression_df.shape}")
    
    # Calculate AUC scores
    print("Calculating AUC scores...")
    try:
        auc_mtx = aucell(
            expression_df,
            regulons,
            num_workers=1,
            normalize=False
        )
        print(f"AUC matrix shape: {auc_mtx.shape}")
        
    except Exception as e:
        print(f"Error with pySCENIC aucell: {e}")
        print("Using simplified AUC calculation...")
        auc_mtx = calculate_simple_auc(expression_df, regulons)
    
    # Save AUC matrix
    auc_mtx.to_csv(auc_file)
    print(f"AUC matrix saved to {auc_file}")
    
    # Calculate binary activity matrix
    print(f"Creating binary activity matrix with threshold {auc_threshold}")
    
    # Calculate percentile thresholds for each regulon
    binary_mtx = auc_mtx.copy()
    for regulon in auc_mtx.index:
        threshold = auc_mtx.loc[regulon].quantile(1 - auc_threshold)
        binary_mtx.loc[regulon] = (auc_mtx.loc[regulon] > threshold).astype(int)
    
    # Save binary matrix
    binary_mtx.to_csv(binary_file)
    print(f"Binary activity matrix saved to {binary_file}")
    
    # Print summary statistics
    print("\nAUC calculation summary:")
    print(f"  Regulons: {auc_mtx.shape[0]}")
    print(f"  Cells: {auc_mtx.shape[1]}")
    print(f"  Mean AUC: {auc_mtx.values.mean():.3f}")
    print(f"  Active cells per regulon (mean): {binary_mtx.sum(axis=1).mean():.1f}")

class SimpleRegulon:
    """Simple regulon class for compatibility with pySCENIC"""
    def __init__(self, name, targets):
        self.name = name
        self.gene2weight = {gene: 1.0 for gene in targets}
        
    @property
    def genes(self):
        return list(self.gene2weight.keys())

def calculate_simple_auc(expression_df, regulons):
    """
    Simple AUC calculation when pySCENIC fails
    """
    print("Using simplified AUC calculation...")
    
    n_cells = expression_df.shape[1]
    auc_data = []
    
    for regulon in regulons:
        regulon_genes = [g for g in regulon.genes if g in expression_df.index]
        
        if len(regulon_genes) == 0:
            # No genes found, create zero scores
            auc_scores = pd.Series(0.0, index=expression_df.columns)
        else:
            # Calculate mean expression of regulon genes
            regulon_expression = expression_df.loc[regulon_genes].mean(axis=0)
            
            # Rank transform to get AUC-like scores
            auc_scores = regulon_expression.rank(pct=True)
        
        auc_data.append(auc_scores)
    
    # Create AUC matrix
    auc_mtx = pd.DataFrame(
        auc_data,
        index=[r.name for r in regulons],
        columns=expression_df.columns
    )
    
    return auc_mtx

if __name__ == "__main__":
    main()
