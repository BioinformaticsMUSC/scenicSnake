#!/usr/bin/env python3
"""
Create UMAP visualization colored by regulon activity
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
    cell_type_column = snakemake.params.cell_type_column
    
    selected_regulons = snakemake.params.selected_regulons
    
    print("Creating UMAP plots colored by regulon activity...")
    
    # Load AUC matrix
    print(f"Loading AUC matrix from {auc_file}")
    auc_df = pd.read_csv(auc_file, index_col=0)
    
    # Load metadata with UMAP coordinates
    print(f"Loading metadata from {metadata_file}")
    adata = sc.read_h5ad(metadata_file)
    split_column = snakemake.config['loom_preparation']['split_condition']
    # Only filter by split condition if it's specified and not empty
    if split_column and split_column.strip() and snakemake.params.split_value != "all":
        adata = adata[adata.obs[split_column] == snakemake.params.split_value].copy()
    
    # Compute UMAP if not present
    if 'X_umap' not in adata.obsm:
        print("Computing UMAP coordinates...")
        sc.tl.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
    
    # Ensure cell order matches
    common_cells = list(set(auc_df.index) & set(adata.obs_names))
    auc_df = auc_df.loc[common_cells, :]
    adata = adata[common_cells]
    
    print(f"Plotting data for {len(common_cells)} cells")
    
    # Get UMAP coordinates
    umap_coords = adata.obsm['X_umap']
    
    # Select regulons to plot
    if selected_regulons and len(selected_regulons) > 0:
        # Use user-specified regulons
        available_regulons = [r for r in selected_regulons if r in auc_df.columns]
        if len(available_regulons) == 0:
            print("Warning: None of the selected regulons found in data")
            # Fall back to top variable regulons
            regulon_variance = auc_df.var(axis=0).sort_values(ascending=False)
            plot_regulons = regulon_variance.head(6).index.tolist()
        else:
            plot_regulons = available_regulons[:6]  # Limit to 6 for visualization
    else:
        # Auto-select most variable regulons
        regulon_variance = auc_df.var(axis=0).sort_values(ascending=False)
        plot_regulons = regulon_variance.head(6).index.tolist()
    
    print(f"Plotting {len(plot_regulons)} regulons: {plot_regulons}")
    
    # Create the plots
    with PdfPages(output_file) as pdf:
        # Overview plot: cell types (if available)
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        
        if cell_type_column in adata.obs.columns:
            cell_types = adata.obs[cell_type_column].astype('category')
            scatter = ax.scatter(
                umap_coords[:, 0], 
                umap_coords[:, 1],
                c=cell_types.cat.codes,
                cmap='tab10',
                s=1,
                alpha=0.7
            )
            
            # Add legend
            handles = []
            for i, ct in enumerate(cell_types.cat.categories):
                handles.append(plt.Line2D([0], [0], marker='o', color='w', 
                                        markerfacecolor=plt.cm.tab10(i), markersize=8, label=ct))
            ax.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left')
            ax.set_title('UMAP - Cell Types')
        else:
            ax.scatter(umap_coords[:, 0], umap_coords[:, 1], s=1, alpha=0.7, c='gray')
            ax.set_title('UMAP - All Cells')
        
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
        
        # Regulon activity plots
        n_regulons = len(plot_regulons)
        n_cols = 3
        n_rows = (n_regulons + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5 * n_rows))
        if n_rows == 1:
            axes = axes.reshape(1, -1)
        
        for i, regulon in enumerate(plot_regulons):
            row = i // n_cols
            col = i % n_cols
            ax = axes[row, col]
            
            # Get regulon activity scores
            activity_scores = auc_df.loc[common_cells, regulon]
            
            # Create scatter plot
            scatter = ax.scatter(
                umap_coords[:, 0],
                umap_coords[:, 1],
                c=activity_scores,
                cmap='viridis',
                s=1,
                alpha=0.7
            )
            
            # Add colorbar
            cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
            cbar.set_label('AUC Score')
            
            ax.set_title(f'{regulon}')
            ax.set_xlabel('UMAP 1')
            ax.set_ylabel('UMAP 2')
        
        # Hide empty subplots
        for i in range(n_regulons, n_rows * n_cols):
            row = i // n_cols
            col = i % n_cols
            axes[row, col].axis('off')
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
        
        # Summary plot: distribution of regulon activities
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()
        
        for i, regulon in enumerate(plot_regulons[:6]):
            activity_scores = auc_df.loc[common_cells, regulon]

            axes[i].hist(activity_scores, bins=50, alpha=0.7, edgecolor='black')
            axes[i].set_xlabel('AUC Score')
            axes[i].set_ylabel('Number of Cells')
            axes[i].set_title(f'{regulon}\n(mean: {activity_scores.mean():.3f})')
            axes[i].grid(True, alpha=0.3)
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    print(f"UMAP plots saved to {output_file}")

if __name__ == "__main__":
    main()
