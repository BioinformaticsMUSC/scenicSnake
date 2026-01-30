#!/bin/bash

# Build and run SCENIC workflow in Docker
# Usage: ./docker-run.sh [command]

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}🐳 SCENIC Docker Workflow Runner${NC}"
echo "=================================="

# Build the Docker image
echo -e "${YELLOW}Building Docker image...${NC}"
docker build -t scenic-snakemake:latest .

# Create necessary directories
echo -e "${YELLOW}Creating directories...${NC}"
mkdir -p data results logs config

# Default command is help
COMMAND=${1:-"--help"}

case $COMMAND in
    "help"|"--help"|"-h")
        echo -e "${GREEN}Available commands:${NC}"
        echo "  ./docker-run.sh help          - Show this help"
        echo "  ./docker-run.sh test          - Test the workflow setup"
        echo "  ./docker-run.sh dry-run       - Perform a dry run (Docker containers)"
        echo "  ./docker-run.sh run           - Run the complete workflow (Docker containers)"
        echo "  ./docker-run.sh shell         - Start interactive shell"
        echo "  ./docker-run.sh clean         - Clean intermediate files"
        echo "  ./docker-run.sh jupyter       - Start Jupyter notebook"
        echo ""
        echo -e "${GREEN}Docker Configuration:${NC}"
        echo "  Uses Docker containers instead of conda environments"
        echo "  Each rule runs in isolated Docker containers"
        echo "  Requires Docker socket access for container orchestration"
        echo ""
        echo -e "${GREEN}Volume mounts:${NC}"
        echo "  ./data    -> /data (input data)"
        echo "  ./results -> /opt/scenic/results (outputs)"
        echo "  ./config  -> /opt/scenic/config (configuration)"
        echo "  ./logs    -> /opt/scenic/logs (log files)"
        ;;
    
    "test")
        echo -e "${YELLOW}Testing workflow setup...${NC}"
        docker run --rm \
            -v "$(pwd)/data:/data" \
            -v "$(pwd)/config:/opt/scenic/config:ro" \
            scenic-snakemake:latest \
            conda run -n scenic python test_setup.py
        ;;
    
    "dry-run")
        echo -e "${YELLOW}Performing dry run...${NC}"
        docker run --rm \
            -v "$(pwd)/data:/data" \
            -v "$(pwd)/results:/opt/scenic/results" \
            -v "$(pwd)/config:/opt/scenic/config:ro" \
            -v "$(pwd)/logs:/opt/scenic/logs" \
            -v /var/run/docker.sock:/var/run/docker.sock \
            -v /tmp:/tmp \
            scenic-snakemake:latest \
            conda run -n scenic snakemake -n --cores 1 --use-singularity \
            --singularity-args "--bind /tmp:/tmp --bind /data:/data"
        ;;
    
    "run")
        echo -e "${YELLOW}Running SCENIC workflow...${NC}"
        docker run --rm \
            -v "$(pwd)/data:/data" \
            -v "$(pwd)/results:/opt/scenic/results" \
            -v "$(pwd)/config:/opt/scenic/config:ro" \
            -v "$(pwd)/logs:/opt/scenic/logs" \
            -v /var/run/docker.sock:/var/run/docker.sock \
            -v /tmp:/tmp \
            scenic-snakemake:latest \
            conda run -n scenic snakemake --cores 8 --use-singularity \
            --singularity-args "--bind /tmp:/tmp --bind /data:/data" \
            --printshellcmds
        ;;
    
    "shell")
        echo -e "${YELLOW}Starting interactive shell...${NC}"
        docker run --rm -it \
            -v "$(pwd)/data:/data" \
            -v "$(pwd)/results:/opt/scenic/results" \
            -v "$(pwd)/config:/opt/scenic/config:ro" \
            -v "$(pwd)/logs:/opt/scenic/logs" \
            scenic-snakemake:latest \
            conda run -n scenic bash
        ;;
    
    "clean")
        echo -e "${YELLOW}Cleaning intermediate files...${NC}"
        docker run --rm \
            -v "$(pwd)/results:/opt/scenic/results" \
            scenic-snakemake:latest \
            conda run -n scenic snakemake clean
        ;;
    
    "jupyter")
        echo -e "${YELLOW}Starting Jupyter notebook...${NC}"
        echo -e "${GREEN}Access at: http://localhost:8888${NC}"
        mkdir -p notebooks
        docker run --rm -it \
            -p 8888:8888 \
            -v "$(pwd)/data:/data" \
            -v "$(pwd)/results:/opt/scenic/results" \
            -v "$(pwd)/notebooks:/opt/scenic/notebooks" \
            scenic-snakemake:latest \
            conda run -n scenic jupyter lab \
            --ip=0.0.0.0 \
            --port=8888 \
            --no-browser \
            --allow-root \
            --NotebookApp.token=''
        ;;
    
    *)
        echo -e "${RED}Unknown command: $COMMAND${NC}"
        echo "Use './docker-run.sh help' for available commands"
        exit 1
        ;;
esac
