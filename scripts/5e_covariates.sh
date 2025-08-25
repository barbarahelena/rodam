#!/bin/bash
#SBATCH -c 68
#SBATCH --mem=24G
#SBATCH --time=4:00:00
#SBATCH -p 'genoa'
eval "$(conda shell.bash hook)"
conda activate xgb
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name age \
    -path /projects/0/prjs0784/rodam/covariatemodels/age \
    -x reg \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name bmi \
    -path /projects/0/prjs0784/rodam/covariatemodels/bmi \
    -x reg \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name bpmed \
    -path /projects/0/prjs0784/rodam/covariatemodels/bpmed \
    -x class \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name bristol \
    -path /projects/0/prjs0784/rodam/covariatemodels/bristol \
    -x reg \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name carbohydrates \
    -path /projects/0/prjs0784/rodam/covariatemodels/carbohydrates \
    -x reg \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name crp \
    -path /projects/0/prjs0784/rodam/covariatemodels/crp \
    -x reg \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name earlylifeurban \
    -path /projects/0/prjs0784/rodam/covariatemodels/earlylifeurban \
    -x class \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json
python /projects/0/prjs0784/rodam/scripts/XGBeast.py \
    -name fat \
    -path /projects/0/prjs0784/rodam/covariatemodels/fat \
    -x class \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/rodam/scripts/param_grid_mb.json