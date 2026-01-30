#!/usr/bin/env python3
"""
Test script to validate SCENIC workflow setup
"""

import sys
import os
import subprocess
from pathlib import Path

def check_conda_env():
    """Check if scenic conda environment exists"""
    try:
        result = subprocess.run(['conda', 'env', 'list'], capture_output=True, text=True)
        if 'scenic' in result.stdout:
            print("✓ Conda environment 'scenic' found")
            return True
        else:
            print("✗ Conda environment 'scenic' not found")
            print("  Run: conda env create -f envs/scenic.yaml")
            return False
    except FileNotFoundError:
        print("✗ Conda not found in PATH")
        return False

def check_config_files():
    """Check if configuration files exist"""
    config_files = [
        'config/config.yaml',
        'config/samples.tsv'
    ]
    
    all_exist = True
    for file in config_files:
        if os.path.exists(file):
            print(f"✓ Configuration file {file} exists")
        else:
            print(f"✗ Configuration file {file} missing")
            all_exist = False
    
    return all_exist

def check_snakemake():
    """Check if Snakemake is available"""
    try:
        result = subprocess.run(['snakemake', '--version'], capture_output=True, text=True)
        version = result.stdout.strip()
        print(f"✓ Snakemake version {version} found")
        return True
    except FileNotFoundError:
        print("✗ Snakemake not found")
        print("  Install with: conda install -c bioconda snakemake")
        return False

def test_workflow_syntax():
    """Test Snakemake workflow syntax"""
    try:
        result = subprocess.run(['snakemake', '-n', '-q'], capture_output=True, text=True)
        if result.returncode == 0:
            print("✓ Snakemake workflow syntax is valid")
            return True
        else:
            print("✗ Snakemake workflow syntax error:")
            print(result.stderr)
            return False
    except Exception as e:
        print(f"✗ Error testing workflow: {e}")
        return False

def main():
    print("SCENIC Snakemake Workflow - Setup Test")
    print("=" * 50)
    
    checks = [
        check_conda_env,
        check_config_files,
        check_snakemake,
        test_workflow_syntax
    ]
    
    results = []
    for check in checks:
        results.append(check())
        print()
    
    if all(results):
        print("🎉 All checks passed! Your SCENIC workflow is ready to run.")
        print("\nNext steps:")
        print("1. Update config/samples.tsv with your data paths")
        print("2. Adjust parameters in config/config.yaml if needed")
        print("3. Run: ./singularity-run.sh run")
    else:
        print("❌ Some checks failed. Please fix the issues above.")
        sys.exit(1)

if __name__ == "__main__":
    main()
