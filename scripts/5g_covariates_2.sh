#!/bin/bash
#SBATCH -c 68
#SBATCH --mem=24G
#SBATCH --time=4:00:00
#SBATCH -p 'genoa'
eval "$(conda shell.bash hook)"
conda activate xgb
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name fibre \
    -path /projects/0/prjs0784/rodam/covariatemodels/fibre \
    -x reg \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name gfr \
    -path /projects/0/prjs0784/rodam/covariatemodels/gfr \
    -x reg \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name kcal \
    -path /projects/0/prjs0784/rodam/covariatemodels/kcal \
    -x reg \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name ldl \
    -path /projects/0/prjs0784/rodam/covariatemodels/ldl \
    -x reg \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name occupationmanual \
    -path /projects/0/prjs0784/rodam/covariatemodels/occupationmanual \
    -x class \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name physicalactivity \
    -path /projects/0/prjs0784/rodam/covariatemodels/physicalactivity \
    -x reg \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name probiotics \
    -path /projects/0/prjs0784/rodam/covariatemodels/probiotics \
    -x class \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name proteins \
    -path /projects/0/prjs0784/rodam/covariatemodels/proteins \
    -x reg \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name sodium \
    -path /projects/0/prjs0784/rodam/covariatemodels/sodium \
    -x reg \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name women \
    -path /projects/0/prjs0784/rodam/covariatemodels/women \
    -x class \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json