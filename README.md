# SCENIC Snakemake Workflow

A comprehensive Snakemake workflow for running SCENIC (Single-Cell rEgulatory Network Inference and Clustering) analysis on single-cell RNA-seq data.

## Overview

SCENIC is a computational method to infer gene regulatory networks (GRNs) from single-cell RNA-seq data and to identify cell states. This workflow implements the complete SCENIC pipeline using Snakemake for reproducible and scalable analysis.

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
- `config/samples.tsv`: Sample information and file paths

### 4. Prepare Your Data

Ensure your single-cell data is in one of these formats:
- AnnData (.h5ad)
- CSV/TSV files
- 10X Genomics format

Update `config/samples.tsv` with your file paths:

```tsv
sample_id	file_path	condition
sample1	data/sample1.h5ad	control
sample2	data/sample2.h5ad	treatment
```

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

```yaml
# Data filtering
filtering:
  min_genes_per_cell: 200
  min_cells_per_gene: 3
  max_genes_per_cell: 5000
  max_mitochondrial_percent: 20

# SCENIC parameters
scenic:
  species: "homo_sapiens"
  n_jobs: 4
  min_genes_per_regulon: 10
  auc_threshold: 0.05
```

### Species Support

Currently supported species:
- `homo_sapiens` (Human)
- `mus_musculus` (Mouse)
- `drosophila_melanogaster` (Fly)

## Output Files

### Core Results
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

## Resource Requirements

### Computational Requirements
- **CPU**: 4-16 cores recommended
- **Memory**: 16-64 GB RAM (depends on dataset size)
- **Storage**: 10-100 GB (depends on dataset size)

### Time Estimates
- Small dataset (1K cells): 1-2 hours
- Medium dataset (10K cells): 4-8 hours
- Large dataset (100K cells): 12-24 hours

## Troubleshooting

### Common Issues

1. **Memory errors during GENIE3**
   - Reduce the number of genes by increasing filtering thresholds
   - Increase memory allocation in cluster configuration

2. **Long runtime**
   - Increase the number of cores (`scenic.n_jobs`)
   - Consider using a smaller gene set for initial testing

3. **Empty regulons**
   - Check species parameter matches your data
   - Verify motif collection is appropriate
   - Lower `min_genes_per_regulon` threshold

### Performance Optimization

1. **Use cluster execution**:
```bash
snakemake --cluster "sbatch --time={resources.time} --mem={resources.mem}" --cores 100
```

2. **Adjust resource allocation** in `config/config.yaml`

3. **Use solid-state storage** for better I/O performance

## Citation

If you use this workflow, please cite:

- **SCENIC**: Aibar et al. Nature Methods (2017)
- **Snakemake**: Köster & Rahmann, Bioinformatics (2012)
- **scanpy**: Wolf et al. Genome Biology (2018)

## License

This workflow is released under the MIT License.

## Support

For questions and issues:
1. Check the [troubleshooting section](#troubleshooting)
2. Open an issue on GitHub
3. Consult the SCENIC documentation: https://pyscenic.readthedocs.io/

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request with detailed description
