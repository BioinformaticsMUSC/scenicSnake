#!/bin/bash

# Smart update script for HPC deployment
# This script updates the workflow code while preserving your local config files

set -e

echo "🔄 Updating SCENIC workflow from GitHub..."

# Option 1: Stash local config changes, pull, then restore
if [ "$1" = "stash" ]; then
    echo "Using git stash method..."
    git stash push config/ -m "Local config files"
    git pull origin main
    git stash pop
    echo "✅ Updated with config files restored"

# Option 2: Pull only specific directories/files (recommended)
elif [ "$1" = "selective" ]; then
    echo "Using selective update method..."
    git fetch origin main
    
    # Update only workflow files, not config
    git checkout origin/main -- Snakefile
    git checkout origin/main -- scripts/
    git checkout origin/main -- envs/
    git checkout origin/main -- Singularity.def
    git checkout origin/main -- singularity-run.sh
    git checkout origin/main -- docker-run.sh
    git checkout origin/main -- *.md
    git checkout origin/main -- *.yml
    
    echo "✅ Updated workflow files only, configs preserved"

# Option 3: Backup configs, pull everything, restore configs
elif [ "$1" = "backup" ]; then
    echo "Using backup method..."
    
    # Backup current configs
    cp -r config/ config_backup/
    
    # Pull updates
    git pull origin main
    
    # Restore configs
    cp -r config_backup/* config/
    rm -rf config_backup/
    
    echo "✅ Updated with configs backed up and restored"

else
    echo "Usage: ./hpc-update.sh [stash|selective|backup]"
    echo ""
    echo "Methods:"
    echo "  stash     - Use git stash to temporarily save config changes"
    echo "  selective - Update only workflow files, preserve all configs"
    echo "  backup    - Backup configs, pull everything, restore configs"
    echo ""
    echo "Recommended: ./hpc-update.sh selective"
fi