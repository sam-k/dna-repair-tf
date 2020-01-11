#!/bin/bash
#SBATCH --job-name mut-profile-SKCM
#SBATCH --mail-user sdk18@duke.edu
#SBATCH --mail-type ALL
#SBATCH --time 12:00:00
#SBATCH -c 2
#SBATCH --output logs/mut-profile-SKCM-noDHS.cl.out
#SBATCH --error logs/mut-profile-SKCM-noDHS.cl.err

### Run on full data using cluster

module load bedtools2

MUT_DATASET=SKCM-US
TFBS_DATASET=skcm

MUT_FILE="../datasets/simple_somatic_mutation.open.${MUT_DATASET}.tsv"
# TFBS_FILE="../datasets/proximalTFBS-DHS_${TFBS_DATASET}.bed"
TFBS_NODHS_FILE="../datasets/proximalTFBS-noDHS_${TFBS_DATASET}.bed"

# Transform TFBSs into TFBS centers ±1000 bp.
TFBS_CENTER="data/proximalTFBS-DHS_${TFBS_DATASET}_center1000.bed"
awk '{center=int(($2+$3)/2); $2=center-1000; $3=center+1000; print $1"\t"$2"\t"$3"\t"$4}' $TFBS_FILE |
  sort -V > $TFBS_CENTER

MUT_FILE_PREFIX="data/ssm.open.${MUT_DATASET}"
MUT_PREP="${MUT_FILE_PREFIX}_prepped.bed"
MUT_INTR="${MUT_FILE_PREFIX}_intersect.bed"
MUT_CENT="${MUT_FILE_PREFIX}_centered.bed"
cut -f9-11,16,17 $MUT_FILE | # select
  sed -e 1d | # sort
  sort -V |
  sed -e $'s/\t/>/4' | # preprocess
  sed -e 's/^/chr/' |
  uniq > $MUT_PREP # remove duplicates
bedtools intersect -a $MUT_PREP -b $TFBS_CENTER -wa -wb | # intersect with TFBS centers
  cut -f1-2,4,6,8 |
  awk '{print $1"\t"$2"\t"$2"\t"$3"\t"($2-$4-1000)"\t"$5}' > $MUT_INTR
cut -f1,5 $MUT_INTR |
  awk '{print $1"\t"$2"\t"$2}' > $MUT_CENT
