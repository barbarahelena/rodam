#!/bin/bash
#SBATCH -c 68
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH -p 'genoa'
eval "$(conda shell.bash hook)"
conda activate xgb
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name urban_ams \
    -path /projects/0/prjs0784/rodam/urban_ams_dich \
    -x class \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name urban_ams \
    -path /projects/0/prjs0784/rodam/urban_ams_dich \
    -x class \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json \
    -permute