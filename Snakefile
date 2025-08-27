"""
SCENIC (Single-Cell rEgulatory Network Inference and Clustering) Snakemake Workflow
================================================================================

This workflow implements the SCENIC algorithm for single-cell RNA-seq data analysis.
SCENIC infers gene regulatory networks and identifies cell states.

Author: Bryan Granger
Date: August 4, 2025
"""

import pandas as pd
from snakemake.utils import min_version

# Minimum Snakemake version
min_version("7.0")

# Configuration
configfile: "config/config.yaml"

# Load sample information
split_column = config['loom_preparation'].get("split_condition", None)
split_values = config['loom_preparation'].get("split_values", ["all"])
print(split_values)
def get_regulatory_feature_dbs():
    """Get regulatory feature databases from config"""
    return " ".join(config["scenic"].get("regulatory_feature_dbs", []))

def get_split_filenames(pattern):
    """
    Generate file names based on split_values.
    If split_values is empty, return pattern with split='all'.
    Otherwise, expand pattern for each split value.
    """
    if split_values:
        return expand(pattern, split=split_values)
    else:
        return [pattern.format(split="all")]

# Example patterns:
# "results/scenic/{split}_adjacencies.tsv"
# "results/plots/{split}_regulon_heatmap.pdf"
# "results/reports/{split}_scenic_report.html"


# Define rule all (final outputs)
rule all:
    input:
        # SCENIC core outputs
        expand("results/scenic/{split_value}_adjacencies.csv", split_value=split_values),
        expand("results/scenic/{split_value}_regulons.csv", split_value=split_values),
        expand("results/scenic/{split_value}_auc_matrix.csv", split_value=split_values),

        # Visualization outputs
        expand("results/plots/{split_value}_regulon_heatmap.pdf", split_value=split_values),
        expand("results/plots/{split_value}_umap_regulon_activity.pdf", split_value=split_values),
        expand("results/plots/{split_value}_rss_plot.pdf", split_value=split_values),

        # Final report
        expand("results/reports/{split_value}_scenic_report.html", split_value=split_values)
        


# SCENIC core rules
rule create_expression_matrix:
    """Prepare expression matrix for SCENIC"""
    input:
        config["loom_preparation"]["adata_file_path"]
    output:
        "results/scenic/{split_value}_expression_matrix.loom"
    params:
        split_value = "{split_value}"
    conda:
        "envs/scenic.yaml"
    script:
        "scripts/03_prepare_expression_matrix.py"

rule grn_inference:
    """Infer gene regulatory network using GRNBOOST2"""
    input:
        "results/scenic/{split_value}_expression_matrix.loom"
    output:
        "results/scenic/{split_value}_adjacencies.csv"
    params:
        n_jobs = config["scenic"]["n_jobs"],
        seed = config["scenic"]["seed"],
        path_to_TF_list = config["scenic"]["path_to_TF_list"]
    threads: config["scenic"]["n_jobs"]
    conda:
        "envs/scenic.yaml"
    shell:
        """
        pyscenic grn \
            {input} \
            {params.path_to_TF_list} \
            --num_workers {params.n_jobs} \
            --seed {params.seed} \
            -o {output}
        """


rule create_regulons:
    """Create regulons from adjacencies using cisTarget"""
    input:
        adjacencies = "results/scenic/{split_value}_adjacencies.csv",
        expression = "results/scenic/{split_value}_expression_matrix.loom"
    output:
        regulons = "results/scenic/{split_value}_regulons.csv",
    params:
        database_fname = get_regulatory_feature_dbs(),
        annotations_fname = config["scenic"]["annotations_fname"],
        min_genes = config["scenic"]["min_genes_per_regulon"],
        auc_threshold = config["scenic"]["auc_threshold"]
    conda:
        "envs/scenic.yaml"
    shell:
        """
        pyscenic ctx \
            {input.adjacencies} \
            {params.database_fname} \
            --annotations_fname {params.annotations_fname} \
            --expression_mtx_fname {input.expression} \
            --output {output.regulons} \
            --auc_threshold {params.auc_threshold} \
            --cell_id_attribute CellAnno \
            --min_genes {params.min_genes}
        """

rule calculate_auc:
    """Calculate AUC scores for regulons"""
    input:
        regulons = "results/scenic/{split_value}_regulons.csv",
        expression = "results/scenic/{split_value}_expression_matrix.loom"
    output:
        auc_matrix = "results/scenic/{split_value}_auc_matrix.csv",
    params:
        auc_threshold = config["scenic"]["auc_threshold"]
    conda:
        "envs/scenic.yaml"
    shell:
        """
        pyscenic aucell \
            {input.expression} \
            {input.regulons} \
            -o {output.auc_matrix} \
            --auc_threshold {params.auc_threshold}
        """

# Visualization rules
rule plot_regulon_heatmap:
    """Create heatmap of regulon activity"""
    input:
        auc_matrix = "results/scenic/{split_value}_auc_matrix.csv",
        metadata = config["loom_preparation"]["adata_file_path"],
        rss_scores = "results/scenic/{split_value}_rss_scores.csv",
    output:
        "results/plots/{split_value}_regulon_heatmap.pdf"
    params:
        top_regulons = config["visualization"]["top_regulons"],
        split_value = "{split_value}",
        split_column = config["loom_preparation"]["split_condition"],
        cell_type_column = config["visualization"]["cell_type_column"]
    conda:
        "envs/scenic.yaml"
    script:
        "scripts/07_plot_heatmap.py"

rule plot_umap_regulon_activity:
    """Plot UMAP colored by regulon activity"""
    input:
        auc_matrix = "results/scenic/{split_value}_auc_matrix.csv",
        metadata = config["loom_preparation"]["adata_file_path"]
    output:
        "results/plots/{split_value}_umap_regulon_activity.pdf"
    params:
        selected_regulons = config["visualization"]["selected_regulons"],
        split_value = "{split_value}",
        split_column = config["loom_preparation"]["split_condition"],
        cell_type_column = config["visualization"]["cell_type_column"]
    conda:
        "envs/scenic.yaml"
    script:
        "scripts/08_plot_umap.py"

rule calculate_rss:
    """Calculate Regulon Specificity Score (RSS)"""
    input:
        auc_matrix = "results/scenic/{split_value}_auc_matrix.csv",
        metadata = config["loom_preparation"]["adata_file_path"]
    output:
        rss_scores = "results/scenic/{split_value}_rss_scores.csv",
        rss_plot = "results/plots/{split_value}_rss_plot.pdf"
    params:
        cell_type_column = config["visualization"]["cell_type_column"],
        split_value = "{split_value}",
        split_column = config["loom_preparation"]["split_condition"]
    conda:
        "envs/scenic.yaml"
    script:
        "scripts/09_calculate_rss.py"

# Report generation
rule generate_report:
    """Generate final HTML report"""
    input:
        auc_matrix = "results/scenic/{split_value}_auc_matrix.csv",
        regulons = "results/scenic/{split_value}_regulons.csv",
        rss_scores = "results/scenic/{split_value}_rss_scores.csv",
        metadata = config["loom_preparation"]["adata_file_path"],
        plots = [
            "results/plots/{split_value}_regulon_heatmap.pdf",
            "results/plots/{split_value}_umap_regulon_activity.pdf",
            "results/plots/{split_value}_rss_plot.pdf"
        ]
    output:
        "results/reports/{split_value}_scenic_report.html"
    params:
        split= "{split_value}",
        cell_type_column = config["visualization"]["cell_type_column"]
    conda:
        "envs/scenic.yaml"
    script:
        "scripts/10_generate_report.py"

# Utility rules
rule clean:
    """Clean intermediate files"""
    shell:
        """
        rm -rf results/preprocessing/
        rm -rf results/scenic/expression_matrix.csv
        """

rule clean_all:
    """Clean all output files"""
    shell:
        """
        rm -rf results/
        """
