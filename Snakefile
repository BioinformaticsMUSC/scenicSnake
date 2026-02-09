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

# Container configuration for Singularity/Docker compatibility
# You can modify this to use different container sources:
# - "docker://condaforge/mambaforge:latest" for conda-based execution (requires internet)
# - "docker://biocontainers/pyscenic:0.12.1--pyhdfd78af_0" for pre-built SCENIC
# - "scenic-workflow.sif" for local Singularity image (recommended for HPC)
# - "docker://yourusername/scenic-snakemake:latest" for custom Docker Hub image
CONTAINER_IMAGE = "scenic-workflow.sif"

# Load sample information
samples_df = pd.read_csv(config["samples"], sep="\t", comment="#")
samples = samples_df.set_index("sample_id", drop=False).to_dict("index")
sample_ids = list(samples.keys())

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
        # SCENIC core outputs for each sample
        expand("results/scenic/{sample_id}_adjacencies.csv", sample_id=sample_ids),
        expand("results/scenic/{sample_id}_regulons.csv", sample_id=sample_ids),
        expand("results/scenic/{sample_id}_auc_matrix.csv", sample_id=sample_ids),

        # Visualization outputs for each sample
        expand("results/plots/{sample_id}_regulon_heatmap.pdf", sample_id=sample_ids),
        expand("results/plots/{sample_id}_umap_regulon_activity.pdf", sample_id=sample_ids),
        expand("results/plots/{sample_id}_rss_plot.pdf", sample_id=sample_ids),

        # Final report for each sample
        expand("results/reports/{sample_id}_scenic_report.html", sample_id=sample_ids)
        


# SCENIC core rules
rule create_expression_matrix:
    """Prepare expression matrix for SCENIC"""
    input:
        h5ad_file = lambda wildcards: samples[wildcards.sample_id]["file_path"]
    output:
        "results/scenic/{sample_id}_expression_matrix.loom"
    params:
        sample_id = "{sample_id}"
    container:
        CONTAINER_IMAGE
    script:
        "scripts/03_prepare_expression_matrix.py"

rule grn_inference:
    """Infer gene regulatory network using GRNBOOST2"""
    input:
        "results/scenic/{sample_id}_expression_matrix.loom"
    output:
        "results/scenic/{sample_id}_adjacencies.csv"
    params:
        n_jobs = config["scenic"]["n_jobs"],
        seed = config["scenic"]["seed"],
        path_to_TF_list = config["scenic"]["path_to_TF_list"]
    threads: config["scenic"]["n_jobs"]
    container:
        CONTAINER_IMAGE
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
        adjacencies = "results/scenic/{sample_id}_adjacencies.csv",
        expression = "results/scenic/{sample_id}_expression_matrix.loom"
    output:
        regulons = "results/scenic/{sample_id}_regulons.csv",
    params:
        database_fname = get_regulatory_feature_dbs(),
        annotations_fname = config["scenic"]["annotations_fname"],
        min_genes = config["scenic"]["min_genes_per_regulon"],
        auc_threshold = config["scenic"]["auc_threshold"]
    container:
        CONTAINER_IMAGE
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
            --min_genes {params.min_genes} \
            --num_workers {config["scenic"]["n_jobs"]}
        """

rule calculate_auc:
    """Calculate AUC scores for regulons"""
    input:
        regulons = "results/scenic/{sample_id}_regulons.csv",
        expression = "results/scenic/{sample_id}_expression_matrix.loom"
    output:
        auc_matrix = "results/scenic/{sample_id}_auc_matrix.csv",
    params:
        auc_threshold = config["scenic"]["auc_threshold"]
    container:
        CONTAINER_IMAGE
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
        auc_matrix = "results/scenic/{sample_id}_auc_matrix.csv",
        metadata = lambda wildcards: samples[wildcards.sample_id]["file_path"],
        rss_scores = "results/scenic/{sample_id}_rss_scores.csv",
    output:
        "results/plots/{sample_id}_regulon_heatmap.pdf"
    params:
        top_regulons = config["visualization"]["top_regulons"],
        sample_id = "{sample_id}",
        cell_type_column = config["visualization"]["cell_type_column"]
    container:
        CONTAINER_IMAGE
    script:
        "scripts/07_plot_heatmap.py"

rule plot_umap_regulon_activity:
    """Plot UMAP colored by regulon activity"""
    input:
        auc_matrix = "results/scenic/{sample_id}_auc_matrix.csv",
        metadata = lambda wildcards: samples[wildcards.sample_id]["file_path"]
    output:
        "results/plots/{sample_id}_umap_regulon_activity.pdf"
    params:
        selected_regulons = config["visualization"]["selected_regulons"],
        sample_id = "{sample_id}",
        cell_type_column = config["visualization"]["cell_type_column"]
    container:
        CONTAINER_IMAGE
    script:
        "scripts/08_plot_umap.py"

rule calculate_rss:
    """Calculate Regulon Specificity Score (RSS)"""
    input:
        auc_matrix = "results/scenic/{sample_id}_auc_matrix.csv",
        metadata = lambda wildcards: samples[wildcards.sample_id]["file_path"]
    output:
        rss_scores = "results/scenic/{sample_id}_rss_scores.csv",
        rss_plot = "results/plots/{sample_id}_rss_plot.pdf"
    params:
        sample_id = "{sample_id}",
        cell_type_column = config["visualization"]["cell_type_column"]
    container:
        CONTAINER_IMAGE
    script:
        "scripts/09_calculate_rss.py"

# Report generation
rule generate_report:
    """Generate final HTML report"""
    input:
        auc_matrix = "results/scenic/{sample_id}_auc_matrix.csv",
        regulons = "results/scenic/{sample_id}_regulons.csv",
        rss_scores = "results/scenic/{sample_id}_rss_scores.csv",
        metadata = lambda wildcards: samples[wildcards.sample_id]["file_path"],
        plots = [
            "results/plots/{sample_id}_regulon_heatmap.pdf",
            "results/plots/{sample_id}_umap_regulon_activity.pdf",
            "results/plots/{sample_id}_rss_plot.pdf"
        ]
    output:
        "results/reports/{sample_id}_scenic_report.html"
    params:
        sample_id = "{sample_id}",
        cell_type_column = config["visualization"]["cell_type_column"]
    container:
        CONTAINER_IMAGE
    script:
        "scripts/10_generate_report.py"

# Utility rules
rule clean:
    """Clean intermediate files"""
    shell:
        """
        rm -rf results/preprocessing/
        rm -rf results/scenic/*_expression_matrix.loom
        """

rule clean_all:
    """Clean all output files"""
    shell:
        """
        rm -rf results/
        """
