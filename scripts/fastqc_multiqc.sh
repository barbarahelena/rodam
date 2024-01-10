#!/bin/bash
#SBATCH -c 24
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH -p 'genoa'
eval "$(conda shell.bash hook)"
conda activate qc
mkdir /home/bverhaar/projects/rodam/fastqc
fastqc \
        /home/bverhaar/projects/rodam/data/*.fastq.gz \
        --threads 24 \
        -o /home/bverhaar/projects/rodam/fastqc
mkdir /home/bverhaar/projects/rodam/multiqc
multiqc \
        /home/bverhaar/projects/rodam/fastqc \
        -o /home/bverhaar/projects/rodam/multiqc