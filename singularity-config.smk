# Snakemake Configuration for Singularity
# This file provides container configurations that work well with Singularity

# For better Singularity compatibility, we can use these container options:

# Option 1: Use a public Docker Hub image with SCENIC dependencies
container_scenic = "docker://biocontainers/pyscenic:0.12.1--pyhdfd78af_0"

# Option 2: Use a conda-forge based image and rely on conda for dependencies
container_conda = "docker://condaforge/mambaforge:latest"

# Option 3: Use the local Singularity image file (built from Docker)
container_local_sif = "scenic-workflow.sif"

# Option 4: Use Docker Hub image after pushing local image
container_dockerhub = "docker://bryangr/scenic-snakemake:latest"  # Replace with your Docker Hub username