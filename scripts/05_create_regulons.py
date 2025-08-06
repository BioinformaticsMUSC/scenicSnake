#!/usr/bin/env python3
"""
Create regulons from adjacencies using cisTarget motif enrichment
"""

import pandas as pd
import numpy as np
import json
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from pyscenic.cisTarget import cisTarget
import os

def main():
    # Get parameters from snakemake
    adjacencies_file = snakemake.input.adjacencies
    expression_file = snakemake.input.expression
    regulons_file = snakemake.output.regulons
    motif_file = snakemake.output.motif_enrichment
    
    species = snakemake.params.species
    motif_collection = snakemake.params.motif_collection
    min_genes = snakemake.params.min_genes
    
    print(f"Creating regulons for species: {species}")
    print(f"Using motif collection: {motif_collection}")
    
    # Load adjacencies
    print(f"Loading adjacencies from {adjacencies_file}")
    adjacencies = pd.read_csv(adjacencies_file, sep='\t')
    
    # Load expression matrix
    print(f"Loading expression matrix from {expression_file}")
    expression_df = pd.read_csv(expression_file, index_col=0)
    
    # Create modules from adjacencies
    print("Creating co-expression modules...")
    modules = list(modules_from_adjacencies(
        adjacencies, 
        expression_df,
        thresholds=[0.75, 0.80, 0.85, 0.90, 0.95],  # Multiple thresholds
        top_n_targets=[50, 100, 200],  # Different numbers of top targets
        top_n_regulators='all',
        min_genes=min_genes
    ))
    
    print(f"Found {len(modules)} co-expression modules")
    
    # Define motif databases based on species
    motif_databases = get_motif_databases(species, motif_collection)
    
    if not motif_databases:
        print(f"Warning: No motif databases found for {species}")
        # Create empty regulons as fallback
        regulons = []
    else:
        # Perform cisTarget motif enrichment
        print("Performing cisTarget motif enrichment...")
        try:
            df = prune2df(
                modules,
                motif_databases,
                motif_threshold=0.001,
                auc_threshold=0.05,
                nes_threshold=3.0,
                rank_threshold=5000,
                num_workers=1
            )
            
            # Save motif enrichment results
            df.to_csv(motif_file, index=False)
            print(f"Motif enrichment results saved to {motif_file}")
            
            # Create regulons
            regulons = df2regulons(df)
            print(f"Created {len(regulons)} regulons")
            
        except Exception as e:
            print(f"Warning: cisTarget failed with error: {e}")
            print("Creating regulons from modules without motif pruning...")
            # Fallback: create regulons directly from modules
            regulons = []
            for module in modules:
                if len(module) >= min_genes:
                    # Create a simple regulon structure
                    regulon = {
                        'name': f"{module[0]}_(+)",  # First gene as TF
                        'targets': list(module[1:]),  # Rest as targets
                        'score': 1.0
                    }
                    regulons.append(regulon)
    
    # Convert regulons to JSON-serializable format
    regulons_dict = {}
    for i, regulon in enumerate(regulons):
        if hasattr(regulon, 'name'):
            # pySCENIC regulon object
            regulon_name = regulon.name
            targets = list(regulon.gene2weight.keys())
        else:
            # Dictionary format
            regulon_name = regulon.get('name', f'Regulon_{i}')
            targets = regulon.get('targets', [])
        
        regulons_dict[regulon_name] = {
            'targets': targets,
            'n_targets': len(targets)
        }
    
    # Save regulons
    with open(regulons_file, 'w') as f:
        json.dump(regulons_dict, f, indent=2)
    
    print(f"Regulons saved to {regulons_file}")
    print(f"Final regulon count: {len(regulons_dict)}")
    
    # Print summary
    if regulons_dict:
        target_counts = [reg['n_targets'] for reg in regulons_dict.values()]
        print(f"Regulon size statistics:")
        print(f"  Mean: {np.mean(target_counts):.1f} targets")
        print(f"  Median: {np.median(target_counts):.1f} targets")
        print(f"  Range: {min(target_counts)}-{max(target_counts)} targets")

def get_motif_databases(species, collection):
    """
    Get motif database paths for different species
    Note: In a real implementation, you would download these databases
    """
    
    # Common database locations (these need to be downloaded separately)
    db_paths = {
        'homo_sapiens': {
            'v10nr_clust': [
                'databases/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather',
                'databases/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
            ]
        },
        'mus_musculus': {
            'v10nr_clust': [
                'databases/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather',
                'databases/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
            ]
        },
        'drosophila_melanogaster': {
            'v10nr_clust': [
                'databases/dm6_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather',
                'databases/dm6_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
            ]
        }
    }
    
    if species not in db_paths:
        print(f"Warning: Species {species} not supported")
        return []
    
    if collection not in db_paths[species]:
        print(f"Warning: Collection {collection} not available for {species}")
        return []
    
    # Check if databases exist
    databases = []
    for db_path in db_paths[species][collection]:
        if os.path.exists(db_path):
            databases.append(db_path)
        else:
            print(f"Warning: Database not found: {db_path}")
            print(f"Please download SCENIC databases from: https://resources.aertslab.org/cistarget/")
    
    return databases

if __name__ == "__main__":
    main()
