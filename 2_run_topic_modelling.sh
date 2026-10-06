#!/bin/bash
#SBATCH --job-name=cis_topic_model
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=244G
#SBATCH --output=cis_topic_model.%j.out
#SBATCH --error=cis_topic_model.%j.err
#SBATCH --time=5-00:00:00
#SBATCH --dependency=afterok:5659337

#$1=n_topics- --topics 5 10 15 20 30 40 50 60 70 80 90 100 150 200 300 \
source /data/stemcell/jwhittle/mambaforge/etc/profile.d/conda.sh

conda activate scenic-plus

mkdir ${SLURM_SUBMIT_DIR}/${1}_tmp/

pycistopic topic_modeling mallet \
	--input scATAC/cistopic_obj.pkl \
	--output scATAC/${1}_model.pkl \
	--temp_dir ${SLURM_SUBMIT_DIR}/${1}_tmp/ \
	--topics ${1} \
	--iterations 500 \
	--alpha 50 \
	--parallel 48 \
	--keep True \
	--seed 555 \
	--mallet_path /data/stemcell/jwhittle/scenicplus/Mallet-202108/bin/mallet

rm -r ${SLURM_SUBMIT_DIR}/${1}_tmp/

conda deactivate
