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
split_column = config.get("split_column", None)
split_values = config.get("split_values", [])

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
        get_split_filenames("results/scenic/{split}_adjacencies.tsv"),
        get_split_filenames("results/scenic/{split}_regulons.json"),
        get_split_filenames("results/scenic/{split}_auc_matrix.csv"),
        get_split_filenames("results/scenic/{split}_binary_regulon_activity.csv"),

        # Visualization outputs
        get_split_filenames("results/plots/{split}_regulon_heatmap.pdf"),
        get_split_filenames("results/plots/{split}_umap_regulon_activity.pdf"),
        get_split_filenames("results/plots/{split}_rss_plot.pdf"),

        # Final report
        get_split_filenames("results/reports/{split}_scenic_report.html")


# SCENIC core rules
rule create_expression_matrix:
    """Prepare expression matrix for SCENIC"""
    input:
        lambda wildcards: config["loom_preparation"]["adata_file_path"]
    output:
        get_split_filenames("results/scenic/{split}_expression_matrix.loom")
    conda:
        "envs/scenic.yaml"
    script:
        "scripts/03_prepare_expression_matrix.py"

rule genie3_inference:
    """Infer gene regulatory network using GENIE3"""
    input:
        get_split_filenames("results/scenic/{split}_expression_matrix.loom")
    output:
        get_split_filenames("results/scenic/{split}_adjacencies.csv")
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
        adjacencies = get_split_filenames("results/scenic/{split}_adjacencies.csv"),
        expression = get_split_filenames("results/scenic/{split}_expression_matrix.loom")
    output:
        regulons = get_split_filenames("results/scenic/{split}_regulons.json"),
        motif_enrichment = get_split_filenames("results/scenic/{split}_motif_enrichment.csv")
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
            --output_motif_enrichment {output.motif_enrichment} \
            --auc_threshold {params.auc_threshold} \
            --cell_id_attribute CellAnno \
            --min_genes {params.min_genes}
        """

rule calculate_auc:
    """Calculate AUC scores for regulons"""
    input:
        regulons = get_split_filenames("results/scenic/{split}_regulons.json"),
        expression = get_split_filenames("results/scenic/{split}_expression_matrix.loom")
    output:
        auc_matrix = get_split_filenames("results/scenic/{split}_auc_matrix.csv"),
        binary_activity = get_split_filenames("results/scenic/{split}_binary_regulon_activity.csv")
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
        auc_matrix = get_split_filenames("results/scenic/{split}_auc_matrix.csv"),
        metadata = config["loom_preparation"]["adata_file_path"]
    output:
        get_split_filenames("results/plots/{split}_regulon_heatmap.pdf")
    params:
        top_regulons = config["visualization"]["top_regulons"],
        split_value = (lambda split_col: (lambda wildcards: wildcards.split if split_col else "all"))(split_column)
    conda:
        "envs/scenic.yaml"
    script:
        "scripts/07_plot_heatmap.py"

rule plot_umap_regulon_activity:
    """Plot UMAP colored by regulon activity"""
    input:
        auc_matrix = get_split_filenames("results/scenic/{split}_auc_matrix.csv"),
        metadata = config["loom_preparation"]["adata_file_path"]
    output:
        get_split_filenames("results/plots/{split}_umap_regulon_activity.pdf")
    params:
        selected_regulons = config["visualization"]["selected_regulons"]
    conda:
        "envs/scenic.yaml"
    script:
        "scripts/08_plot_umap.py"

rule calculate_rss:
    """Calculate Regulon Specificity Score (RSS)"""
    input:
        auc_matrix = get_split_filenames("results/scenic/{split}_auc_matrix.csv"),
        metadata = config["loom_preparation"]["adata_file_path"]
    output:
        rss_scores = get_split_filenames("results/scenic/{split}_rss_scores.csv"),
        rss_plot = get_split_filenames("results/plots/{split}_rss_plot.pdf")
    params:
        cell_type_column = config["visualization"]["cell_type_column"]
    conda:
        "envs/scenic.yaml"
    script:
        "scripts/09_calculate_rss.py"

# Report generation
rule generate_report:
    """Generate final HTML report"""
    input:
        auc_matrix = get_split_filenames("results/scenic/{split}_auc_matrix.csv"),
        regulons = get_split_filenames("results/scenic/{split}_regulons.json"),
        rss_scores = get_split_filenames("results/scenic/{split}_rss_scores.csv"),
        plots = [
            get_split_filenames("results/plots/{split}_regulon_heatmap.pdf"),
            get_split_filenames("results/plots/{split}_umap_regulon_activity.pdf"),
            get_split_filenames("results/plots/{split}_rss_plot.pdf")
        ]
    output:
        get_split_filenames("results/reports/{split}_scenic_report.html")
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
