#!/bin/bash
#SBATCH --job-name test
#SBATCH --mail-user sdk18@duke.edu
#SBATCH --mail-type END,FAIL
#SBATCH --time 12:00:00
#SBATCH -c 2
#SBATCH --output logs/test.cl.out
#SBATCH --error logs/test.cl.err

### Run snippets for testing purposes

# module load bedtools2

# MUT_DATASET=BLCA
# TFBS_DATASET=blca

# MUT_FILE="../datasets/simple_somatic_mutation.open.${MUT_DATASET}.tsv"
# TSS_FILE="../datasets/refseq_TSS_hg19_170929.bed"
# TFBS_FILE="../datasets/distalTFBS-${DHS}_${TFBS_DATASET}.bed"

# TSS_PROC="./tss-processed.bed"
# MUT_INTR="./${MUT_DATASET}_intr_tss.bed"

# ## TSS_FILE:
# #  Transcription start sites on genome
# #  1. chromosome
# #  2. chromosome_start
# #  3. chromosome_end
# #  4. name
# #  5. score
# #  6. strand

# cut -f1-3 "${TSS_FILE}" |
#   sort -V |
#   uniq > "${TSS_PROC}"

# cut -f9-11,16-17 "${MUT_FILE}" | # select cols
#   sort -V | # sort
#   sed -e $'s/\t/>/4' | # preprocess to BED format
#   sed -e 's/^/chr/' |
#   uniq | # remove duplicates
#   bedtools closest -a - -b "${TSS_PROC}" -D ref |
#   sort -V > "${MUT_INTR}"

ENH_FILE="../datasets/permissive_enhancers.bed"
ENH_PROC="./test_permissive_enhancers_proc.bed"
cut -f1-3 "${ENH_FILE}" |
  tail -n +2 |
  sort -V > "${ENH_PROC}"