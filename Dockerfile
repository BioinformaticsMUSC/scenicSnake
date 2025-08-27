# Use mambaforge as base image for better conda performance
FROM condaforge/mambaforge:latest

# Set metadata
LABEL maintainer="Bryan Granger"
LABEL description="SCENIC Snakemake Workflow - Single-Cell Regulatory Network Analysis"
LABEL version="1.0"

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    wget \
    curl \
    build-essential \
    libssl-dev \
    libffi-dev \
    libhdf5-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Create working directory
WORKDIR /opt/scenic

# Copy environment file first (for better Docker layer caching)
COPY envs/scenic.yaml /opt/scenic/envs/

# Create conda environment
RUN mamba env create -f envs/scenic.yaml && \
    mamba clean -afy

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "scenic", "/bin/bash", "-c"]

# Copy the workflow files
COPY . /opt/scenic/

# Create necessary directories
RUN mkdir -p /opt/scenic/data \
    /opt/scenic/results \
    /opt/scenic/logs

# Set permissions
RUN chmod +x /opt/scenic/scripts/*.py
RUN chmod +x /opt/scenic/test_setup.py

# Add conda environment to PATH
ENV PATH="/opt/conda/envs/scenic/bin:$PATH"

# Expose port for potential web interfaces
EXPOSE 8080

# Set default working directory for user
WORKDIR /data

# Default command
CMD ["conda", "run", "-n", "scenic", "snakemake", "--help"]
