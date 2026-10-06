#!/bin/bash
#SBATCH --job-name=cis_db
#SBATCH --partition=hmem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=2000G
#SBATCH --output=cis_db.%j.out
#SBATCH --error=cis_db.%j.err
#SBATCH --time=2-00:00:00
#SBATCH --dependency=afterok:3951498

source /data/stemcell/jwhittle/mambaforge/etc/profile.d/conda.sh

conda activate scenic-plus

ml apps/bedtools/2.31.0

OUT_DIR="."
DATABASE_PREFIX="hg38_Lambo"
REGION_BED="../scATAC/consensus_peak_calling/consensus_regions.bed"
GENOME_FASTA="/data/stemcell/jwhittle/ref/refdata-cellranger-atac-GRCh38-1.1.0/fasta/genome.fa"
CHROMSIZES="/data/stemcell/jwhittle/ref/refdata-cellranger-atac-GRCh38-1.1.0/fasta/genome.chrom.sizes"
SCRIPT_DIR="/data/stemcell/jwhittle/scenicplus/cistarget_db/create_cisTarget_databases"
CBDIR="/data/stemcell/jwhittle/scenicplus/cistarget_db/aertslab_motif_colleciton/v10nr_clust_public/singletons"
MOTIF_LIST="/data/stemcell/jwhittle/scenicplus/cistarget_db/aertslab_motif_colleciton/v10nr_clust_public/motifs.txt"
CBUST_PATH="/data/stemcell/jwhittle/scenicplus/cistarget_db/cbust"

${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        ${OUT_DIR}/${DATABASE_PREFIX}.fa \
        1000 \
        yes

${SCRIPT_DIR}/create_cistarget_motif_databases.py \
    -f ${OUT_DIR}/${DATABASE_PREFIX}.fa \
    -M ${CBDIR} \
    -m ${MOTIF_LIST} \
    -o ${OUT_DIR}/${DATABASE_PREFIX} \
    -c ${CBUST_PATH} \
    --bgpadding 1000 \
    -t ${SLURM_CPUS_PER_TASK}


conda deactivate
