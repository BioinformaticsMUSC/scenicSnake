#!/usr/bin/env python3
"""
Create heatmap visualization of regulon activity
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from matplotlib.backends.backend_pdf import PdfPages

def main():
    # Get parameters from snakemake
    auc_file = snakemake.input.auc_matrix
    metadata_file = snakemake.input.metadata
    output_file = snakemake.output[0]
    
    top_regulons = snakemake.params.top_regulons
    
    print("Creating regulon activity heatmap...")
    
    # Load AUC matrix
    print(f"Loading AUC matrix from {auc_file}")
    auc_df = pd.read_csv(auc_file, index_col=0)
    
    # Load metadata
    print(f"Loading metadata from {metadata_file}")
    adata = sc.read_h5ad(metadata_file)
    split_column = snakemake.params.split_column
    adata = adata[adata.obs[split_column] == snakemake.params.split_value].copy()
    
    # Ensure cell order matches
    common_cells = list(set(auc_df.columns) & set(adata.obs_names))
    auc_df = auc_df[common_cells]
    adata = adata[common_cells]
    
    print(f"Plotting data for {len(common_cells)} cells and {auc_df.shape[0]} regulons")
    
    # Select top variable regulons
    if top_regulons and top_regulons < auc_df.shape[0]:
        regulon_variance = auc_df.var(axis=1).sort_values(ascending=False)
        top_regulon_names = regulon_variance.head(top_regulons).index
        auc_subset = auc_df.loc[top_regulon_names]
        print(f"Selected top {top_regulons} most variable regulons")
    else:
        auc_subset = auc_df
        print("Using all regulons")
    
    # Create the plot
    with PdfPages(output_file) as pdf:
        # Main heatmap
        fig, axes = plt.subplots(2, 1, figsize=(12, 10), 
                                gridspec_kw={'height_ratios': [1, 4]})
        
        # Cell type annotation bar (if available)
        if 'cell_type' in adata.obs.columns:
            cell_types = adata.obs['cell_type'].astype('category')
            cell_type_colors = dict(zip(
                cell_types.cat.categories,
                plt.cm.Set3(np.linspace(0, 1, len(cell_types.cat.categories)))
            ))
            
            # Create color bar
            colors = [cell_type_colors[ct] for ct in cell_types]
            axes[0].imshow([colors], aspect='auto')
            axes[0].set_xlim(0, len(colors))
            axes[0].set_yticks([])
            axes[0].set_title('Cell Types')
            axes[0].set_xticks([])
        else:
            axes[0].axis('off')
        
        # Regulon activity heatmap
        sns.heatmap(
            auc_subset,
            ax=axes[1],
            cmap='viridis',
            cbar_kws={'label': 'AUC Score'},
            xticklabels=False,
            yticklabels=True
        )
        
        axes[1].set_xlabel('Cells')
        axes[1].set_ylabel('Regulons')
        axes[1].set_title(f'Regulon Activity Heatmap (Top {auc_subset.shape[0]} Regulons)')
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
        
        # Summary statistics plot
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        
        # Regulon activity distribution
        axes[0, 0].hist(auc_subset.values.flatten(), bins=50, alpha=0.7)
        axes[0, 0].set_xlabel('AUC Score')
        axes[0, 0].set_ylabel('Frequency')
        axes[0, 0].set_title('Distribution of AUC Scores')
        
        # Number of active regulons per cell
        active_regulons_per_cell = (auc_subset > auc_subset.quantile(0.95, axis=1).values[:, np.newaxis]).sum(axis=0)
        axes[0, 1].hist(active_regulons_per_cell, bins=30, alpha=0.7)
        axes[0, 1].set_xlabel('Number of Active Regulons')
        axes[0, 1].set_ylabel('Number of Cells')
        axes[0, 1].set_title('Active Regulons per Cell')
        
        # Mean activity per regulon
        mean_activity = auc_subset.mean(axis=1).sort_values(ascending=True)
        y_pos = np.arange(len(mean_activity))
        axes[1, 0].barh(y_pos[-20:], mean_activity.values[-20:])  # Top 20
        axes[1, 0].set_yticks(y_pos[-20:])
        axes[1, 0].set_yticklabels(mean_activity.index[-20:])
        axes[1, 0].set_xlabel('Mean AUC Score')
        axes[1, 0].set_title('Top 20 Most Active Regulons')
        
        # Regulon activity variance
        variance_activity = auc_subset.var(axis=1).sort_values(ascending=True)
        y_pos = np.arange(len(variance_activity))
        axes[1, 1].barh(y_pos[-20:], variance_activity.values[-20:])  # Top 20
        axes[1, 1].set_yticks(y_pos[-20:])
        axes[1, 1].set_yticklabels(variance_activity.index[-20:])
        axes[1, 1].set_xlabel('AUC Score Variance')
        axes[1, 1].set_title('Top 20 Most Variable Regulons')
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    print(f"Heatmap saved to {output_file}")

if __name__ == "__main__":
    main()
