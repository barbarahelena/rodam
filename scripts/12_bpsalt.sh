#!/bin/bash
#SBATCH -c 48
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH -p 'genoa'
eval "$(conda shell.bash hook)"
conda activate xgb
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name sbp \
    -path /projects/0/prjs0784/rodam/sbp \
    -x reg \
    -n 200 \
    -t 48 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name dbp \
    -path /projects/0/prjs0784/rodam/dbp \
    -x reg \
    -n 200 \
    -t 48 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name salt \
    -path /projects/0/prjs0784/rodam/salt \
    -x reg \
    -n 200 \
    -t 48 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json