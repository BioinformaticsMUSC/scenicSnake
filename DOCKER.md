# Docker Usage Guide for SCENIC Snakemake Workflow

This guide explains how to run the SCENIC workflow using Docker for maximum portability and reproducibility.

## Quick Start with Docker

### 1. Build the Docker Image

```bash
# Using the helper script
./docker-run.sh help

# Or using make
make docker-build

# Or manually
docker build -t scenic-snakemake:latest .
```

### 2. Prepare Your Data

Create the necessary directories and place your data:

```bash
mkdir -p data results config logs
# Place your .h5ad file in the data/ directory
# Update config/config.yaml with your settings
```

### 3. Run the Workflow

```bash
# Test the setup
./docker-run.sh test

# Dry run to check workflow
./docker-run.sh dry-run

# Run the complete workflow
./docker-run.sh run
```

## Docker Commands

### Using the Helper Script

The `docker-run.sh` script provides convenient commands:

```bash
./docker-run.sh help          # Show available commands
./docker-run.sh test          # Test workflow setup
./docker-run.sh dry-run       # Perform dry run
./docker-run.sh run           # Run complete workflow
./docker-run.sh shell         # Interactive shell
./docker-run.sh clean         # Clean intermediate files
./docker-run.sh jupyter       # Start Jupyter notebook
```

### Using Make Commands

```bash
make docker-build     # Build Docker image
make docker-run       # Run workflow in Docker
make docker-shell     # Start interactive shell
make docker-test      # Test workflow
make docker-clean     # Clean Docker resources
```

### Using Docker Compose

```bash
# Start the main workflow container
docker-compose up scenic

# Start with Jupyter notebook
docker-compose --profile jupyter up

# Run in background
docker-compose up -d scenic
```

## Volume Mounts

The Docker setup automatically mounts these directories:

| Host Directory | Container Path | Purpose |
|---------------|----------------|---------|
| `./data` | `/data` | Input data files |
| `./results` | `/opt/scenic/results` | Output files |
| `./config` | `/opt/scenic/config` | Configuration (read-only) |
| `./logs` | `/opt/scenic/logs` | Log files |

## Manual Docker Commands

If you prefer to run Docker commands manually:

```bash
# Build the image
docker build -t scenic-snakemake:latest .

# Run workflow
docker run --rm \
  -v "$(pwd)/data:/data" \
  -v "$(pwd)/results:/opt/scenic/results" \
  -v "$(pwd)/config:/opt/scenic/config:ro" \
  scenic-snakemake:latest \
  conda run -n scenic snakemake --cores 8 --use-conda

# Interactive shell
docker run --rm -it \
  -v "$(pwd)/data:/data" \
  -v "$(pwd)/results:/opt/scenic/results" \
  -v "$(pwd)/config:/opt/scenic/config:ro" \
  scenic-snakemake:latest \
  conda run -n scenic bash
```

## Jupyter Notebook in Docker

Start a Jupyter notebook server for interactive analysis:

```bash
# Using helper script
./docker-run.sh jupyter

# Using docker-compose
docker-compose --profile jupyter up

# Access at: http://localhost:8888
```

## Cluster Usage

### Singularity (for HPC clusters)

Convert Docker image to Singularity:

```bash
# Build Singularity image from Docker
singularity build scenic.sif docker://scenic-snakemake:latest

# Run on cluster
singularity exec scenic.sif conda run -n scenic snakemake --cores 8
```

### SLURM Integration

Create a SLURM job script:

```bash
#!/bin/bash
#SBATCH --job-name=scenic
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

module load singularity

singularity exec scenic.sif \
  conda run -n scenic snakemake \
  --cores 8 --use-conda
```

## Troubleshooting

### Common Issues

1. **Permission Errors**
   ```bash
   # Fix file permissions
   sudo chown -R $USER:$USER results/
   ```

2. **Port Already in Use** (Jupyter)
   ```bash
   # Use different port
   docker run -p 8889:8888 ...
   ```

3. **Out of Disk Space**
   ```bash
   # Clean Docker system
   make docker-clean
   docker system prune -a
   ```

### Debug Mode

Run container in debug mode:

```bash
./docker-run.sh shell
# Inside container:
conda run -n scenic snakemake -n --verbose
```

## Performance Tips

1. **Allocate More Memory**: Increase Docker Desktop memory allocation
2. **Use SSD Storage**: Mount data from SSD for better I/O performance
3. **Parallel Processing**: Adjust `--cores` parameter based on your system
4. **Cache Conda**: Reuse conda cache between runs

## Security Considerations

- The container runs as root for compatibility
- Mounted directories are accessible from within container
- Use read-only mounts for sensitive configuration files
- Consider using Docker secrets for credentials

## Building Custom Images

Customize the Dockerfile for your needs:

```dockerfile
# Add custom packages
RUN conda install -c conda-forge mypackage

# Install system dependencies
RUN apt-get install -y mytool

# Add custom scripts
COPY my_scripts/ /opt/scenic/scripts/
```

Then rebuild:

```bash
docker build -t my-scenic:latest .
```
