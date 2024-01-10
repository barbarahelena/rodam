#!/bin/bash
#SBATCH -c 24
#SBATCH --mem=10G
#SBATCH --time=12:00:00
#SBATCH -p 'genoa'
eval "$(conda shell.bash hook)"
conda activate nextflow
mkdir /projects/0/prjs0784/rodam/results
nextflow run /home/bverhaar/projects/nf-core-vsearchpipeline/main.nf \
        -profile singularity,snellius \
        --input /projects/0/prjs0784/rodam/data/Samplesheet.csv \
        --outdir /projects/0/prjs0784/rodam/results \
        --skip_primers \
        --email 'bjh.verhaar@gmail.com' \
        --treetool 'iqtree' \
        -resume
