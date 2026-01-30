#!/bin/bash

# SCENIC Workflow Mode Selector
# Choose between Docker containers or Conda environments

set -e

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${GREEN}🐍 SCENIC Workflow Execution Mode${NC}"
echo "=================================="
echo ""
echo "Choose how to run the SCENIC workflow:"
echo ""
echo -e "${BLUE}1) Docker Container Mode${NC}"
echo "   • Each rule runs in isolated Docker containers"
echo "   • Uses Snakemake's --use-singularity flag"
echo "   • Good for development and local execution"
echo "   • Requires Docker daemon and socket access"
echo ""
echo -e "${BLUE}2) Singularity Mode (HPC Recommended)${NC}"
echo "   • Native Singularity execution"
echo "   • Better for HPC/cluster environments"  
echo "   • No root privileges needed after build"
echo "   • Maximum portability and security"
echo ""
echo -e "${BLUE}3) Conda Mode (Legacy)${NC}"
echo "   • Uses conda environments within Docker"
echo "   • Uses Snakemake's --use-conda flag"  
echo "   • Faster startup, less isolation"
echo "   • May have environment conflicts"
echo ""

read -p "Enter your choice (1, 2, or 3): " choice

case $choice in
    1)
        echo -e "${YELLOW}Setting up Docker Container Mode...${NC}"
        echo ""
        echo "Running workflow with Docker containers:"
        echo "• Each rule executes in its own container"
        echo "• Using configurable container images"
        echo ""
        echo "Commands to run:"
        echo "  ./docker-run.sh dry-run    # Test with containers"
        echo "  ./docker-run.sh run        # Execute with containers"
        echo ""
        echo -e "${GREEN}Docker container mode selected!${NC}"
        ;;
    2)
        echo -e "${YELLOW}Setting up Singularity Mode...${NC}"
        echo ""
        echo "Running workflow with Singularity:"
        echo "• Native Singularity execution"
        echo "• Ideal for HPC environments"
        echo ""
        echo "Commands to run:"
        echo "  ./singularity-run.sh build     # Build Singularity image"
        echo "  ./singularity-run.sh dry-run   # Test with Singularity"
        echo "  ./singularity-run.sh run       # Execute with Singularity"
        echo ""
        echo -e "${GREEN}Singularity mode selected!${NC}"
        ;;
    3)
        echo -e "${YELLOW}Setting up Conda Mode...${NC}"
        echo ""
        echo "To use Conda mode, modify docker-run.sh:"
        echo "• Change --use-singularity to --use-conda"
        echo "• Remove Docker socket mount"
        echo ""
        echo "Or run manually:"
        echo "  docker run --rm -v \$(pwd):/data scenic-snakemake:latest \\"
        echo "    conda run -n scenic snakemake --use-conda --cores 8"
        echo ""
        echo -e "${GREEN}Conda mode guidance provided!${NC}"
        ;;
    *)
        echo "Invalid choice. Please run again and select 1, 2, or 3."
        exit 1
        ;;
esac