#!/bin/bash
#SBATCH --job-name=scplus
#SBATCH --partition=hmem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=72
#SBATCH --mem=3000G
#SBATCH --output=scplus.%j.out
#SBATCH --error=scplus.%j.err
#SBATCH --time=2-00:00:00
#SBATCH --dependency=afterok:4248322

source /data/stemcell/jwhittle/mambaforge/etc/profile.d/conda.sh

conda activate scenic-plus

mkdir tmp

snakemake --cores 72

rm -r tmp

conda deactivate
