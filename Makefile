# Makefile for SCENIC Snakemake Workflow
# Provides convenient shortcuts for common tasks

.PHONY: help setup test clean dry-run run report

# Default target
help:
	@echo "SCENIC Snakemake Workflow"
	@echo "========================"
	@echo ""
	@echo "Available targets:"
	@echo "  setup     - Create conda environment and setup workflow"
	@echo "  test      - Test workflow configuration and setup"
	@echo "  dry-run   - Perform a dry run to check workflow"
	@echo "  run       - Run the complete SCENIC workflow"
	@echo "  report    - Generate Snakemake HTML report"
	@echo "  clean     - Clean intermediate files"
	@echo "  clean-all - Remove all output files"
	@echo ""

# Setup conda environment
setup:
	@echo "Setting up SCENIC workflow..."
	conda env create -f envs/scenic.yaml
	@echo "Environment created. Activate with: conda activate scenic"

# Test workflow setup
test:
	@echo "Testing workflow setup..."
	python test_setup.py

# Dry run
dry-run:
	@echo "Performing dry run..."
	snakemake -n --cores 1

# Run workflow with reasonable defaults
run:
	@echo "Running SCENIC workflow..."
	snakemake --cores 8 --use-conda --printshellcmds

# Run with cluster support (SLURM example)
run-cluster:
	@echo "Running on cluster..."
	snakemake --cluster "sbatch --time=24:00:00 --mem=16G --cpus-per-task={threads}" \
		--cores 100 --use-conda --latency-wait 60

# Generate workflow report
report:
	@echo "Generating workflow report..."
	snakemake --report results/reports/workflow_report.html

# Clean intermediate files
clean:
	@echo "Cleaning intermediate files..."
	snakemake clean

# Clean all output files
clean-all:
	@echo "Cleaning all output files..."
	snakemake clean_all
	rm -rf results/
	rm -rf .snakemake/

# Create example data directory structure
create-data-dirs:
	@echo "Creating data directory structure..."
	mkdir -p data/{raw,processed}
	mkdir -p results/{preprocessing,scenic,plots,reports,qc}

# Quick start for new users
quickstart: setup create-data-dirs
	@echo "Quick start complete!"
	@echo ""
	@echo "Next steps:"
	@echo "1. Place your data files in data/raw/"
	@echo "2. Update config/samples.tsv with your file paths"
	@echo "3. Run 'make test' to verify setup"
	@echo "4. Run 'make dry-run' to check workflow"
	@echo "5. Run 'make run' to execute the workflow"
