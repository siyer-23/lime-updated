#!/bin/bash

#SBATCH --job-name=lime_rmats_catalog
#SBATCH --account=cphg-millerlab-vip
#SBATCH --partition=standard
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

#SBATCH --output=slurm/lime_rmats_catalog.%j.out
#SBATCH --error=slurm/lime_rmats_catalog.%j.err

#SBATCH --mail-user=bhe3qa@virginia.edu
#SBATCH --mail-type=END,FAIL


set -euo pipefail


PROJECT_DIR="/sfs/gpfs/tardis/project/cphg-millerlab/ShayI/30_lime-refreshed_2026"


cd "${PROJECT_DIR}"


module load miniforge


export PYTHONPATH="${PROJECT_DIR}/src:${PYTHONPATH:-}"


python -m lime.rmats.catalog \
    --manifest config/rmats_runs.tsv \
    --output-dir results/rmats_catalog \
    --count-mode JC \
    --fdr 0.05 \
    --min-abs-dpsi 0