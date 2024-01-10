#!/bin/bash
#SBATCH -c 48
#SBATCH --mem=32G
#SBATCH --time=00:45:00
#SBATCH -p 'genoa'
eval "$(conda shell.bash hook)"
conda activate xgb
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name renin \
    -path /projects/0/prjs0784/rodam/renin \
    -x reg \
    -n 200 \
    -t 48 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name aldosterone \
    -path /projects/0/prjs0784/rodam/aldo \
    -x reg \
    -n 200 \
    -t 48 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name arr \
    -path /projects/0/prjs0784/rodam/ARR \
    -x reg \
    -n 200 \
    -t 48 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json