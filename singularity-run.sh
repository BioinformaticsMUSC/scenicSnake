#!/bin/bash

# Singularity-specific runner for SCENIC workflow
# This script runs the workflow directly with Singularity (not Docker-in-Singularity)

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}🔬 SCENIC Singularity Workflow Runner${NC}"
echo "====================================="

# Build Singularity image from definition file (no Docker required)
echo -e "${YELLOW}Checking for Singularity image...${NC}"
if [ ! -f "scenic-workflow.sif" ]; then
    if [ -f "Singularity.def" ]; then
        echo -e "${YELLOW}Building Singularity image from definition file...${NC}"
        singularity build scenic-workflow.sif Singularity.def
        echo -e "${GREEN}Singularity image built: scenic-workflow.sif${NC}"
    else
        echo -e "${YELLOW}No Singularity.def found. Trying to pull from Docker Hub...${NC}"
        # Try to pull from a public registry (requires internet)
        if singularity build scenic-workflow.sif docker://condaforge/mambaforge:latest; then
            echo -e "${GREEN}Base image pulled. You may need to install additional dependencies.${NC}"
        else
            echo -e "${RED}Failed to build image. Please ensure Singularity.def exists or internet connectivity.${NC}"
            exit 1
        fi
    fi
else
    echo -e "${GREEN}Using existing Singularity image: scenic-workflow.sif${NC}"
fi

# Create necessary directories
echo -e "${YELLOW}Creating directories...${NC}"
mkdir -p data results logs config

# Default command is help
COMMAND=${1:-"--help"}

case $COMMAND in
    "help"|"--help"|"-h")
        echo -e "${GREEN}Available commands:${NC}"
        echo "  ./singularity-run.sh help          - Show this help"
        echo "  ./singularity-run.sh build         - Build Singularity image"
        echo "  ./singularity-run.sh test          - Test the workflow setup"
        echo "  ./singularity-run.sh dry-run       - Perform a dry run"
        echo "  ./singularity-run.sh run           - Run the complete workflow"
        echo "  ./singularity-run.sh shell         - Start interactive shell"
        echo "  ./singularity-run.sh clean         - Clean intermediate files"
        echo ""
        echo -e "${GREEN}Singularity Configuration:${NC}"
        echo "  Uses native Singularity execution (no Docker required)"
        echo "  Each rule runs in Singularity containers"
        echo "  Bind mounts: ./data, ./results, ./config, ./logs"
        ;;
    
    "build")
        echo -e "${YELLOW}Building Singularity image...${NC}"
        if [ -f "Singularity.def" ]; then
            singularity build --force scenic-workflow.sif Singularity.def
            echo -e "${GREEN}Image built from Singularity.def${NC}"
        else
            echo -e "${YELLOW}No Singularity.def found. Building from Docker Hub...${NC}"
            singularity build --force scenic-workflow.sif docker://condaforge/mambaforge:latest
            echo -e "${YELLOW}Note: You may need to install additional dependencies manually.${NC}"
        fi
        ;;
    
    "test")
        echo -e "${YELLOW}Testing workflow setup...${NC}"
        singularity exec \
            --bind "$(pwd)/data:/data" \
            --bind "$(pwd)/config:/opt/scenic/config" \
            scenic-workflow.sif \
            conda run -n scenic python /opt/scenic/test_setup.py
        ;;
    
    "dry-run")
        echo -e "${YELLOW}Performing dry run...${NC}"
        singularity exec \
            --bind "$(pwd)/data:/data" \
            --bind "$(pwd)/results:/opt/scenic/results" \
            --bind "$(pwd)/config:/opt/scenic/config" \
            --bind "$(pwd)/logs:/opt/scenic/logs" \
            --bind /tmp:/tmp \
            --pwd /opt/scenic \
            scenic-workflow.sif \
            conda run -n scenic snakemake -n --cores 1 --use-singularity --singularity-prefix "$(pwd)/.snakemake/singularity"
        ;;
    
    "run")
        echo -e "${YELLOW}Running SCENIC workflow...${NC}"
        singularity exec \
            --bind "$(pwd)/data:/data" \
            --bind "$(pwd)/results:/opt/scenic/results" \
            --bind "$(pwd)/config:/opt/scenic/config" \
            --bind "$(pwd)/logs:/opt/scenic/logs" \
            --bind /tmp:/tmp \
            --pwd /opt/scenic \
            scenic-workflow.sif \
            conda run -n scenic snakemake --cores 8 --use-singularity --singularity-prefix "$(pwd)/.snakemake/singularity" --printshellcmds
        ;;
    
    "shell")
        echo -e "${YELLOW}Starting interactive shell...${NC}"
        singularity shell \
            --bind "$(pwd)/data:/data" \
            --bind "$(pwd)/results:/opt/scenic/results" \
            --bind "$(pwd)/config:/opt/scenic/config" \
            --bind "$(pwd)/logs:/opt/scenic/logs" \
            scenic-workflow.sif
        ;;
    
    "clean")
        echo -e "${YELLOW}Cleaning intermediate files...${NC}"
        singularity exec \
            --bind "$(pwd)/results:/opt/scenic/results" \
            scenic-workflow.sif \
            conda run -n scenic snakemake clean
        ;;
    
    *)
        echo -e "${RED}Unknown command: $COMMAND${NC}"
        echo "Use './singularity-run.sh help' for available commands"
        exit 1
        ;;
esac