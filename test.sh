#!/bin/bash
#SBATCH --job-name test
#SBATCH --mail-user sdk18@duke.edu
#SBATCH --mail-type END,FAIL
#SBATCH --time 12:00:00
#SBATCH -c 2
#SBATCH --output logs/test.cl.out
#SBATCH --error logs/test.cl.err

### Run snippets for testing purposes

# MUT_DATASET=SKCM-US
# TFBS_DATASET=skcm

# MUT_FILE="./temp_data/sample.ssm.open.${MUT_DATASET}.tsv"
# MUT_FILE_PREFIX="./temp_data/sample.ssm.open.${MUT_DATASET}"
# MUT_PREP="${MUT_FILE_PREFIX}_prepped.bed"
# MUT_INTR="${MUT_FILE_PREFIX}_intersect.bed"
# MUT_CNTR="${MUT_FILE_PREFIX}_centered.bed"

# TFBS_FILE="../datasets/proximalTFBS-DHS_${TFBS_DATASET}.bed"
# TFBS_CNTR="./data/sample.proximalTFBS-DHS_${TFBS_DATASET}_center1000.bed"

# cut -f9-11,16,17,34 $MUT_FILE | # select cols
#   awk '$6=="WGS"' | # get only WGS
#   cut -f1-5 | # remove sequencing_strategy col
#   sort -V | # sort
#   sed -e $'s/\t/>/4' | # preprocess to BED format
#   sed -e 's/^/chr/' |
#   uniq > $MUT_PREP # remove duplicates

module load bedtools2

CDS_FILE="../datasets/coding_exons.bed"
GEN_FILE="../datasets/human.hg38.genome"
NONCODING="./data/supplementary/coding_exons_complement.bed"

grep -P '^chr(\d+|[MXY])\t' "${CDS_FILE}" |
  sort -V |
  bedtools complement -i - -g "${GEN_FILE}" > "${NONCODING}"