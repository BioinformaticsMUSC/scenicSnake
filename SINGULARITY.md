# Singularity Usage Guide for SCENIC Workflow

This guide explains how to run the SCENIC workflow with Singularity for HPC environments and better container management.

## Singularity vs Docker

**Singularity Benefits:**
- Better suited for HPC/cluster environments
- No root privileges required after image build
- Native support for MPI and GPU
- Better integration with schedulers (SLURM, PBS, etc.)
- More secure for multi-user environments

**Docker Benefits:**
- Easier development and testing
- Better ecosystem and tooling
- More familiar to most developers

## Singularity Usage Options

### Option 1: Native Singularity (Recommended for HPC)

Use the dedicated Singularity runner script:

```bash
# Build Singularity image from Docker
./singularity-run.sh build

# Test the workflow
./singularity-run.sh test

# Dry run
./singularity-run.sh dry-run

# Full run
./singularity-run.sh run

# Interactive shell
./singularity-run.sh shell
```

### Option 2: Docker-in-Docker with Singularity

Use the existing Docker script but with Singularity orchestration:

```bash
# This uses Docker to run Snakemake, which then uses Singularity for rules
./docker-run.sh run
```

### Option 3: Direct Singularity Commands

Manual Singularity execution:

```bash
# Build image
singularity build scenic-workflow.sif docker-daemon://scenic-snakemake:latest

# Run workflow
singularity exec \
    --bind $(pwd)/data:/data \
    --bind $(pwd)/results:/opt/scenic/results \
    --bind $(pwd)/config:/opt/scenic/config \
    scenic-workflow.sif \
    snakemake --cores 8 --use-singularity
```

## Container Image Options

The workflow supports multiple container sources (configured in Snakefile):

```python
# Option 1: Conda-based image (requires environment installation)
CONTAINER_IMAGE = "docker://condaforge/mambaforge:latest"

# Option 2: Pre-built SCENIC image
CONTAINER_IMAGE = "docker://biocontainers/pyscenic:0.12.1--pyhdfd78af_0"

# Option 3: Local Singularity image
CONTAINER_IMAGE = "scenic-workflow.sif"

# Option 4: Custom Docker Hub image
CONTAINER_IMAGE = "docker://yourusername/scenic-snakemake:latest"
```

## Building Custom Singularity Image

### From Singularity Definition File

```bash
# Build from definition file
singularity build scenic-workflow.sif Singularity.def
```

### From Docker Image

```bash
# Build local Docker image first
docker build -t scenic-snakemake:latest .

# Convert to Singularity
singularity build scenic-workflow.sif docker-daemon://scenic-snakemake:latest
```

### From Docker Hub

```bash
# Direct from Docker Hub (if you push your image)
singularity build scenic-workflow.sif docker://yourusername/scenic-snakemake:latest
```

## HPC Cluster Integration

### SLURM Example

```bash
#!/bin/bash
#SBATCH --job-name=scenic
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=24:00:00
#SBATCH --mem=32G

module load singularity

# Run SCENIC workflow
singularity exec \
    --bind $PWD/data:/data \
    --bind $PWD/results:/results \
    scenic-workflow.sif \
    snakemake --cores $SLURM_NTASKS --use-singularity
```

### PBS Example

```bash
#!/bin/bash
#PBS -N scenic
#PBS -l nodes=1:ppn=8
#PBS -l walltime=24:00:00
#PBS -l mem=32gb

cd $PBS_O_WORKDIR
module load singularity

singularity exec \
    --bind $PWD/data:/data \
    --bind $PWD/results:/results \
    scenic-workflow.sif \
    snakemake --cores $PBS_NP --use-singularity
```

## Troubleshooting

### Common Issues

1. **Permission Errors**: Ensure proper bind mounts for data directories
2. **Missing Dependencies**: Check container image has all required software
3. **Path Issues**: Use absolute paths for bind mounts
4. **Memory Limits**: Increase memory allocation for large datasets

### Debugging Commands

```bash
# Check container contents
singularity shell scenic-workflow.sif

# Test individual steps
singularity exec scenic-workflow.sif which snakemake
singularity exec scenic-workflow.sif python -c "import scanpy; print(scanpy.__version__)"

# Verbose execution
singularity exec -v scenic-workflow.sif snakemake --cores 1 -n
```

## Performance Considerations

- **Image Size**: Use minimal base images when possible
- **Bind Mounts**: Only bind necessary directories
- **Shared Storage**: Store images on shared filesystem for clusters
- **Caching**: Use `--singularity-prefix` to cache downloaded images
- **Parallel Jobs**: Adjust `--cores` based on available resources