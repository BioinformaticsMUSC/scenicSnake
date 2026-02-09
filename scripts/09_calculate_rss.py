#!/usr/bin/env python3
"""
Calculate Regulon Specificity Score (RSS) for cell type specificity
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from matplotlib.backends.backend_pdf import PdfPages
from pyscenic.rss import regulon_specificity_scores
from scipy.stats import zscore

def main():
    # Get parameters from snakemake
    auc_file = snakemake.input.auc_matrix
    metadata_file = snakemake.input.metadata
    rss_file = snakemake.output.rss_scores
    plot_file = snakemake.output.rss_plot
    
    cell_type_column = snakemake.params.cell_type_column
    
    print("Calculating Regulon Specificity Scores (RSS)...")
    
    # Load AUC matrix
    print(f"Loading AUC matrix from {auc_file}")
    auc_df = pd.read_csv(auc_file, index_col=0)
    
    # Load metadata
    print(f"Loading metadata from {metadata_file}")
    adata = sc.read_h5ad(metadata_file)
    split_column = snakemake.config['loom_preparation']['split_condition']
    # Only filter by split condition if it's specified and not empty
    if split_column and split_column.strip() and snakemake.params.split_value != "all":
        adata = adata[adata.obs[split_column] == snakemake.params.split_value].copy()

    # Ensure cell order matches
    common_cells = list(set(auc_df.index) & set(adata.obs_names))
    auc_df = auc_df.loc[common_cells, :]
    adata = adata[common_cells]
    
    # Check if cell type column exists
    if cell_type_column not in adata.obs.columns:
        print(f"Warning: Cell type column '{cell_type_column}' not found in metadata")
        print("Available columns:", list(adata.obs.columns))
        
        # Try to find alternative column names
        potential_columns = ['cell_type', 'celltype', 'cluster', 'leiden', 'seurat_clusters']
        cell_type_column = None
        for col in potential_columns:
            if col in adata.obs.columns:
                cell_type_column = col
                print(f"Using '{col}' as cell type column")
                break
        
        if cell_type_column is None:
            print("No suitable cell type column found. Creating dummy clusters...")
            # Create dummy clusters based on data structure
            from sklearn.cluster import KMeans
            
            # Use first 2 PCs for clustering if available
            if 'X_pca' in adata.obsm:
                features = adata.obsm['X_pca'][:, :2]
            else:
                # Use highly variable genes
                sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
                features = adata[:, adata.var.highly_variable].X.toarray()
            
            # Perform clustering
            n_clusters = min(10, len(common_cells) // 50)  # Reasonable number of clusters
            kmeans = KMeans(n_clusters=n_clusters, random_state=42)
            clusters = kmeans.fit_predict(features)
            
            adata.obs['dummy_clusters'] = [f'Cluster_{i}' for i in clusters]
            cell_type_column = 'dummy_clusters'
            print(f"Created {n_clusters} dummy clusters")
    
    # Get cell type information
    print(cell_type_column)
    cell_types = adata.obs[cell_type_column].astype(str)
    unique_cell_types = cell_types.unique()
    
    print(f"Found {len(unique_cell_types)} cell types: {unique_cell_types}")
    
    # Calculate RSS scores
    rss_scores = regulon_specificity_scores(auc_df, cell_types)
    
    # Save RSS scores
    rss_scores.to_csv(rss_file)
    print(f"RSS scores saved to {rss_file}")
    
    # Create visualizations
    create_rss_plots(rss_scores, cell_types, plot_file)
    print(f"RSS plots saved to {plot_file}")

def calculate_rss_alternative(auc_df, cell_types):
    """
    Calculate Regulon Specificity Score (RSS)
    RSS measures how specific a regulon is to each cell type
    """
    
    unique_cell_types = cell_types.unique()
    rss_data = []
    
    for cell_type in unique_cell_types:
        # Get cells of this type
        type_mask = cell_types == cell_type
        type_cells = cell_types[type_mask].index
        other_cells = cell_types[~type_mask].index
        
        if len(type_cells) == 0:
            continue
        
        # Calculate mean AUC for this cell type vs others
        type_mean = auc_df.loc[type_cells,:].mean(axis=0)
        other_mean = auc_df.loc[other_cells,:].mean(axis=0) if len(other_cells) > 0 else pd.Series(0, index=auc_df.columns)
        
        # RSS = log2(fold change) with pseudocount
        pseudocount = 0.001
        rss = np.log2((type_mean + pseudocount) / (other_mean + pseudocount))
        
        rss_data.append(pd.Series(rss, name=cell_type))
    
    # Combine into DataFrame
    rss_df = pd.concat(rss_data, axis=1)
    rss_df = rss_df.fillna(0)
    
    return rss_df

def overall_auc_heatmap(auc_df, plot_file, top_n=15):
    """Create overall AUC heatmap"""

    #zscore the auc_df for each regulon
    auc_mean = auc_df.mean(axis=0)
    auc_mean_z = (auc_mean - auc_mean.mean()) / auc_mean.std()
    top_auc = auc_mean_z.sort_values(ascending=False).head(top_n)
    fig, ax = plt.subplots(1, len(unique_conditions), figsize=(12, 8))
    for i, cond in enumerate(unique_conditions):
        sns.heatmap(
            top_auc,
            ax=ax[i],
            annot=True,
            cmap='Reds',
            center=0,
            cbar_kws={'label': 'AUC Score, z-scored'},
            yticklabels=True,
            xticklabels=True
        )
        ax.set_title(f'Overall Mean AUC scores - {cond}')
    ax.set_ylabel('Regulons')
    plt.tight_layout()
    plt.savefig(plot_file)
    plt.close()

def create_rss_plots(rss_scores, cell_types, plot_file):
    """Create RSS visualization plots"""
    
    with PdfPages(plot_file) as pdf:
        # RSS heatmap
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        
        # Clip extreme values for better visualization
        rss_clipped = rss_scores.clip(-3, 3)
        
        sns.heatmap(
            rss_clipped,
            ax=ax,
            cmap='RdBu_r',
            center=0,
            cbar_kws={'label': 'RSS Score (log2 FC)'},
            yticklabels=True,
            xticklabels=True
        )
        
        ax.set_title('Regulon Specificity Scores (RSS)')
        ax.set_xlabel('Cell Types')
        ax.set_ylabel('Regulons')
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
        
        # Top specific regulons per cell type
        n_cell_types = len(rss_scores.columns)
        n_cols = min(3, n_cell_types)
        n_rows = (n_cell_types + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5 * n_rows))
        if n_rows == 1:
            axes = axes.reshape(1, -1)
        elif n_cols == 1:
            axes = axes.reshape(-1, 1)
        
        for i, cell_type in enumerate(rss_scores.columns):
            row = i // n_cols
            col = i % n_cols
            ax = axes[row, col]
            
            # Get top 10 specific regulons for this cell type
            top_regulons = rss_scores[cell_type].sort_values(ascending=False).head(10)
            
            if len(top_regulons) > 0:
                y_pos = np.arange(len(top_regulons))
                ax.barh(y_pos, top_regulons.values, alpha=0.7)
                ax.set_yticks(y_pos)
                ax.set_yticklabels(top_regulons.index)
                ax.set_xlabel('RSS Score')
                ax.set_title(f'{cell_type} - Top Specific Regulons')
                ax.grid(True, alpha=0.3)
        
        # Hide empty subplots
        for i in range(n_cell_types, n_rows * n_cols):
            row = i // n_cols
            col = i % n_cols
            axes[row, col].axis('off')
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
        
        # RSS distribution
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # Overall RSS distribution
        all_rss = rss_scores.values.flatten()
        axes[0].hist(all_rss, bins=50, alpha=0.7, edgecolor='black')
        axes[0].set_xlabel('RSS Score')
        axes[0].set_ylabel('Frequency')
        axes[0].set_title('Distribution of RSS Scores')
        axes[0].axvline(0, color='red', linestyle='--', alpha=0.7, label='No specificity')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        
        # Max RSS per regulon (most specific cell type)
        max_rss = rss_scores.max(axis=1).sort_values(ascending=False)
        axes[1].hist(max_rss, bins=30, alpha=0.7, edgecolor='black')
        axes[1].set_xlabel('Max RSS Score')
        axes[1].set_ylabel('Number of Regulons')
        axes[1].set_title('Maximum Specificity per Regulon')
        axes[1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()

if __name__ == "__main__":
    main()
