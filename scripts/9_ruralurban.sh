#!/bin/bash
#SBATCH -c 48
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH -p 'genoa'
eval "$(conda shell.bash hook)"
conda activate xgb
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name rural_urban \
    -path /projects/0/prjs0784/rodam/rural_urban \
    -x class \
    -n 200 \
    -t 48 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json