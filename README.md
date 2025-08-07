# SCENIC Snakemake Workflow

A comprehensive Snakemake workflow for running SCENIC analysis on single-cell RNA-seq data.

## Overview

This workflow implements the complete SCENIC pipeline using Snakemake, including downstream analysis.

## Features

- **Data preprocessing**: Quality control and filtering of single-cell data
- **Network inference**: Gene regulatory network inference using GENIE3
- **Regulon discovery**: Identification of regulons using cisTarget
- **Activity scoring**: Calculation of regulon activity scores (AUC)
- **Visualization**: Comprehensive plots and reports
- **Scalable**: Parallel processing support
- **Reproducible**: Conda environments and version control

## Workflow Steps

1. **Data Loading & Filtering** (`01_load_and_filter.py`)
   - Load single-cell expression data
   - Filter low-quality cells and genes
   - Calculate QC metrics

2. **Expression Matrix Preparation** (`03_prepare_expression_matrix.py`)
   - Normalize expression data
   - Log-transform values
   - Format for SCENIC input

3. **Gene Regulatory Network Inference** (`04_genie3_inference.py`)
   - Run GENIE3 algorithm
   - Generate gene-gene adjacency matrix

4. **Regulon Creation** (`05_create_regulons.py`)
   - Use cisTarget for motif enrichment
   - Create transcription factor regulons

5. **Activity Scoring** (`06_calculate_auc.py`)
   - Calculate AUC scores for each regulon
   - Generate binary activity matrix

6. **Visualization & Analysis**
   - Regulon activity heatmaps
   - UMAP plots colored by regulon activity
   - Regulon Specificity Score (RSS) analysis

## Quick Start

### 1. Clone and Setup

```bash
git clone <repository-url>
cd scenicSnake
```

### 2. Install Dependencies

Create the conda environment:

```bash
conda env create -f envs/scenic.yaml
conda activate scenic
```

### 3. Configure Your Analysis

Edit the configuration files:

- `config/config.yaml`: Main workflow parameters

### 4. Prepare Your Data

Create an AnnData object (.h5ad) of the preprocessed data. 


### 5. Run the Workflow

```bash
# Dry run to check the workflow
snakemake -n

# Run the complete workflow
snakemake --cores 8 --use-conda

# Run specific steps
snakemake results/scenic/regulons.json --cores 4 --use-conda
```

## Configuration

### Main Parameters (`config/config.yaml`)
Make sure the config.yaml file is updated prior to running.

## Output Files

### Core Results
Files will be saved for each split condition if applicable. 
- `results/scenic/adjacencies.tsv`: Gene-gene adjacency matrix
- `results/scenic/regulons.json`: Discovered regulons
- `results/scenic/auc_matrix.csv`: Regulon activity scores
- `results/scenic/binary_regulon_activity.csv`: Binary regulon activity

### Visualizations
- `results/plots/regulon_heatmap.pdf`: Regulon activity heatmap
- `results/plots/umap_regulon_activity.pdf`: UMAP with regulon overlay
- `results/plots/rss_plot.pdf`: Regulon specificity scores

### Reports
- `results/reports/scenic_report.html`: Comprehensive analysis report

### Performance Optimization

1. **Use cluster execution**:
```bash
snakemake --cluster "sbatch --time={resources.time} --mem={resources.mem}" --cores 32
```

2. **Adjust resource allocation** in `config/config.yaml`


## Citation

If you use this workflow, please cite:

- **SCENIC**: Aibar et al. Nature Methods (2017)
- **Snakemake**: Köster & Rahmann, Bioinformatics (2012)
- **scanpy**: Wolf et al. Genome Biology (2018)

## License

This workflow is released under the MIT License.
