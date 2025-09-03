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

def get_heatmap_data(adata, auc_df, rss_data, cell_type_column, n_regs_per_celltype=5):
    reg_dict = {}
    for ct in rss_data.index:
        reg_dict[ct] = list(rss_data.T[ct].sort_values(ascending=False)[:n_regs_per_celltype].index)
    auc_means = auc_df.groupby(cell_type_column).mean().T
    top_regs_list = []
    for c in reversed(auc_means.columns):
        top_cc = reg_dict[c]
        top_cc_sorted = auc_means.loc[top_cc,c].sort_values(ascending=True).index
        for r in top_cc_sorted:
            if r not in top_regs_list:
                top_regs_list.append(r)
    hm_data = auc_means.loc[top_regs_list,:]
    return hm_data

def main():
    # Load the data
    auc_file_template = snakemake.params.auc_path
    rss_file_template = snakemake.params.rss_path
    adata_file = snakemake.input.adata_file
    split_values = snakemake.params.split_values
    cell_type_column = snakemake.params.cell_type_column

    hm_data_dict = {}
    max_values = []
    for sv in split_values:
        auc_path = auc_file_template.replace("??GROUP??", sv)
        rss_path = rss_file_template.replace("??GROUP??", sv)
        adata = sc.read(adata_file)
        adata_split = adata[adata.obs[snakemake.config['loom_preparation']['split_condition']] == sv].copy()
        auc_df = pd.read_csv(auc_path, index_col=0)
        rss_data = pd.read_csv(rss_path, index_col=0)

        auc_df = pd.concat([auc_df, adata_split.obs[cell_type_column]], axis=1)

        hm_data = get_heatmap_data(adata, auc_df, rss_data, cell_type_column, n_regs_per_celltype=5)
        hm_data_dict[sv] = hm_data
        max_values.append(hm_data.values.max())

    # Plot the heatmaps in one figure
    plt.figure(figsize=(12, 10))
    for i, (sv, hm_data) in enumerate(hm_data_dict.items()):
        plt.subplot(1, len(hm_data_dict), i + 1)
        sns.heatmap(hm_data, cmap="viridis", vmax=max(max_values), yticklabels=1, cbar=i==len(hm_data_dict)-1)
        plt.title(f"Regulon Activity Heatmap - {sv}")
        plt.xlabel("Cell Type")
        plt.xticks(rotation=90, ha='center')
        plt.ylabel("Regulon")
    plt.tight_layout()
    plt.savefig(f"results/plots/OVERALL_group_heatmap.pdf")
    plt.close()

main()