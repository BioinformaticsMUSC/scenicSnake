# HPC Setup Guide (No Docker Required)

This guide is for running SCENIC on HPC systems without Docker access.

## Quick Start (HPC Environment)

### 1. Build Singularity Image

```bash
# Method 1: Build from definition file (recommended)
singularity build scenic-workflow.sif Singularity.def

# Method 2: Pull from Docker Hub (requires internet)
singularity build scenic-workflow.sif docker://condaforge/mambaforge:latest
```

### 2. Test the Image

```bash
# Test that the environment works
singularity exec scenic-workflow.sif conda run -n scenic python -c "import scanpy; print('Success!')"

# Test Snakemake
singularity exec scenic-workflow.sif conda run -n scenic snakemake --help
```

### 3. Prepare Your Data

```bash
# Create directory structure
mkdir -p data results config logs

# Copy your h5ad files to data/
cp /path/to/your/files/*.h5ad data/

# Update config/samples.tsv with your file paths
# Update config/config.yaml with your settings
```

### 4. Run the Workflow

```bash
# Using the helper script
./singularity-run.sh dry-run    # Test first
./singularity-run.sh run        # Full execution

# Or manually
singularity exec \
    --bind $(pwd)/data:/data \
    --bind $(pwd)/results:/opt/scenic/results \
    --bind $(pwd)/config:/opt/scenic/config \
    scenic-workflow.sif \
    conda run -n scenic snakemake --cores 8 --use-conda
```

## SLURM Job Script Example

```bash
#!/bin/bash
#SBATCH --job-name=scenic
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --partition=compute

# Load Singularity module (adjust for your HPC)
module load singularity

# Change to workflow directory
cd $SLURM_SUBMIT_DIR

# Run SCENIC workflow
singularity exec \
    --bind $(pwd)/data:/data \
    --bind $(pwd)/results:/opt/scenic/results \
    --bind $(pwd)/config:/opt/scenic/config \
    --bind $(pwd)/logs:/opt/scenic/logs \
    --bind /tmp:/tmp \
    scenic-workflow.sif \
    conda run -n scenic snakemake \
    --cores $SLURM_NTASKS \
    --use-conda \
    --printshellcmds \
    --keep-going
```

## Alternative: Using Conda on HPC

If you prefer not to use containers:

```bash
# Load conda/mamba module
module load miniconda3

# Create environment from file
mamba env create -f envs/scenic.yaml

# Activate environment
conda activate scenic

# Run workflow
snakemake --cores 8 --use-conda
```

## Pre-built Image Locations

If building is slow, you can use pre-built images:

```bash
# From Docker Hub (if available)
singularity pull docker://biocontainers/pyscenic:0.12.1--pyhdfd78af_0

# Or request your HPC admin to build centrally:
# /shared/containers/scenic-workflow.sif
```

## Troubleshooting

### Build Issues
- **No internet**: Copy the definition file and build on a machine with internet, then transfer
- **Permissions**: Use `--fakeroot` flag if available: `singularity build --fakeroot scenic-workflow.sif Singularity.def`
- **Space**: Use `--tmpdir /path/to/large/tmp` for build temporary files

### Runtime Issues
- **Bind mount errors**: Ensure all paths exist before running
- **Permission denied**: Check file permissions in bind-mounted directories
- **Memory errors**: Increase memory allocation in SLURM script

### Module Loading
```bash
# Common HPC module names
module load singularity/3.8.0
module load apptainer  # newer name for Singularity
```

## File Transfer

For large images, consider:

```bash
# Compress for transfer
gzip scenic-workflow.sif

# Transfer to HPC
scp scenic-workflow.sif.gz user@hpc.institution.edu:~/

# Decompress on HPC
gunzip scenic-workflow.sif.gz
```