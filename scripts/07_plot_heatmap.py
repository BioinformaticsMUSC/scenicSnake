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
    rss_file = snakemake.input.rss_scores
    output_file = snakemake.output[0]
    
    top_regulons = snakemake.params.top_regulons
    cell_type_column = snakemake.params.cell_type_column

    print("Creating regulon activity heatmap...")
    
    # Load AUC matrix
    print(f"Loading AUC matrix from {auc_file}")
    auc_df = pd.read_csv(auc_file, index_col=0)
    rss_df = pd.read_csv(rss_file, index_col=0)
    # Load metadata
    print(f"Loading metadata from {metadata_file}")
    adata = sc.read_h5ad(metadata_file)
    # split_column = snakemake.config['loom_preparation']['split_condition']
    # adata = adata[adata.obs[split_column] == snakemake.params.split_value].copy()
    
    # Ensure cell order matches
    common_cells = list(set(auc_df.index) & set(adata.obs_names))
    auc_df = auc_df.loc[common_cells, :]
    adata = adata[common_cells]
    auc_df_ct = pd.concat([auc_df, adata.obs[cell_type_column]], axis=1)

    print(f"Plotting data for {len(common_cells)} cells and {auc_df.shape[1]} regulons")
    
    # Select top variable regulons
    if top_regulons and top_regulons < auc_df.shape[0]:
        regulon_variance = auc_df.var(axis=0).sort_values(ascending=False)
        top_regulon_names = regulon_variance.head(top_regulons).index
        auc_subset = auc_df.loc[:,top_regulon_names]
        print(f"Selected top {top_regulons} most variable regulons")
    else:
        auc_subset = auc_df
        print("Using all regulons")

    # get top regulons per RSS
    reg_dict = {}
    for ct in rss_df.index:
        reg_dict[ct] = list(rss_df.T[ct].sort_values(ascending=False)[:5].index)

    # Create the plot
    with PdfPages(output_file) as pdf:
        # Main heatmap
        fig, ax = plt.subplots(1, 1, figsize=(12, 10), 
                                )
        
        # Cell type annotation bar (if available)
        # if snakemake.config['loom_preparation']['celltype_column'] in adata.obs.columns:
        #     cell_types = adata.obs[snakemake.config['loom_preparation']['celltype_column']].astype('category')
        #     cell_type_colors = dict(zip(
        #         cell_types.cat.categories,
        #         plt.cm.Set3(np.linspace(0, 1, len(cell_types.cat.categories)))
        #     ))
        #     print(cell_types.cat.categories)
        #     print(cell_type_colors)
            
        #     # Create color bar
        #     colors = [cell_type_colors[ct] for ct in cell_types]
        #     axes[1].imshow([colors], aspect='auto')
        #     axes[1].set_xlim(0, len(colors))
        #     axes[1].set_yticks([])
        #     axes[1].set_title('Cell Types')
        #     axes[1].set_xticks([])
        # else:
        #     axes[1].axis('off')

        # Regulon activity heatmap
        auc_means = auc_df_ct.groupby(cell_type_column).mean().T
        top_regs_list = []
        for c in reversed(auc_means.columns):
            top_cc = reg_dict[c]
            top_cc_sorted = auc_means.loc[top_cc,c].sort_values(ascending=True).index
            for r in top_cc_sorted:
                if r not in top_regs_list:
                    top_regs_list.append(r)
        hm_data = auc_means.loc[top_regs_list,:]
        sns.heatmap(data=hm_data, annot=False, cmap="Reds", cbar_kws={"shrink": 0.35}, ax=ax)
        
        ax.set_ylabel('Regulons')
        ax.set_xlabel('Cells')
        ax.set_title(f'Regulon Activity Heatmap for Top 5 RSS Regulons per Cell Type')

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
        mean_activity = auc_subset.mean(axis=0).sort_values(ascending=True)
        y_pos = np.arange(len(mean_activity))
        axes[1, 0].barh(y_pos[-20:], mean_activity.values[-20:])  # Top 20
        axes[1, 0].set_yticks(y_pos[-20:])
        axes[1, 0].set_yticklabels(mean_activity.index[-20:])
        axes[1, 0].set_xlabel('Mean AUC Score')
        axes[1, 0].set_title('Top 20 Most Active Regulons')
        
        # Regulon activity variance
        variance_activity = auc_subset.var(axis=0).sort_values(ascending=True)
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

if __name__ == "__main__":
    main()
