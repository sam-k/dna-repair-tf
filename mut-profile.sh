#!/bin/bash

### Run locally on small sample of data

MUT_DATASET=SKCM-US
TFBS_DATASET=skcm

MUT_FILE="../datasets/simple_somatic_mutation.open.${MUT_DATASET}.tsv"
TFBS_FILE="../datasets/proximalTFBS-DHS_${TFBS_DATASET}.bed"
# TFBS_DHS_FILE="../datasets/proximalTFBS-noDHS_${TFBS_DATASET}.bed"

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




MUT_FILE=../datasets/simple_somatic_mutation.open.SKCM-US.tsv
TFBS_FILE=../datasets/proximalTFBS-DHS_skcm.bed

# Transform TFBSs into TFBS centers ±1000 bp.
TFBS_CENTER=data/proximalTFBS-DHS_skcm_center1000.bed
awk '{center=int(($2+$3)/2); $2=center-1000; $3=center+1000; print $1"\t"$2"\t"$3"\t"$4}' $TFBS_FILE |
  sort -V > $TFBS_CENTER

# Sample dataset to work with.
MUT_SAMPLE=data/samples/sample.ssm.SKCM-US_0.tsv
head -n 1 $MUT_FILE > $MUT_SAMPLE
sed -e 1d $MUT_FILE | gshuf -n 1000 >> $MUT_SAMPLE

# Select, sort and preprocess columns.
MUT_SAMPLE_CUT=data/samples/sample.ssm.SKCM-US_1_mut_only.tsv #select
cut -f9-11,16,17 $MUT_SAMPLE > $MUT_SAMPLE_CUT
MUT_SAMPLE_SORT=data/samples/sample.ssm.SKCM-US_2_sorted.tsv #sort
(head -n 1 $MUT_SAMPLE_CUT && sed -e 1d $MUT_SAMPLE_CUT | sort -V) > $MUT_SAMPLE_SORT
MUT_SAMPLE_BED=data/samples/sample.ssm.SKCM-US_3_prepped.bed #preprocess
sed -e 1d $MUT_SAMPLE_SORT |
 sed -e $'s/\t/>/4' |
 sed -e 's/^/chr/' > $MUT_SAMPLE_BED
cut -f1-2,4 $MUT_SAMPLE_BED |
 uniq -d
MUT_SAMPLE_UNIQ=data/samples/sample.ssm.SKCM-US_4_unique.bed #remove duplicates
uniq $MUT_SAMPLE_BED > $MUT_SAMPLE_UNIQ

# Intersect with TFBS centers.
MUT_SAMPLE_INT=data/samples/sample.ssm.SKCM-US_5_intersect.bed #intersect
bedtools intersect -a $MUT_SAMPLE_UNIQ -b $TFBS_CENTER -wa -wb |
 cut -f1-2,4,6,8 |
 awk '{print $1"\t"$2"\t"$2"\t"$3"\t"($2-$4-1000)"\t"$5}' > $MUT_SAMPLE_INT
MUT_SAMPLE_CENTER=data/samples/sample.ssm.SKCM-US_6_centered.bed #center
cut -f1,5 $MUT_SAMPLE_INT |
 awk '{print $1"\t"$2"\t"$2}' > $MUT_SAMPLE_CENTER

# Perform operations on actual dataset.
# MUT_FILE_PREFIX=data/ssm.open.SKCM-US
# MUT_PREP="${MUT_FILE_PREFIX}_prepped.bed"
# MUT_DUPE="${MUT_FILE_PREFIX}_duplicated.bed"
# MUT_UNIQ="${MUT_FILE_PREFIX}_unique.bed"
# MUT_INTR="${MUT_FILE_PREFIX}_intersect.bed"
# MUT_CENT="${MUT_FILE_PREFIX}_centered.bed"
# cut -f9-11,16,17 $MUT_FILE |
#   sed -e 1d |
#   sort -V |
#   sed -e $'s/\t/>/4' |
#   sed -e 's/^/chr/' > $MUT_PREP
# cut -f1-2,4 $MUT_PREP |
#   uniq -d > $MUT_DUPE
# uniq $MUT_PREP > $MUT_UNIQ
# bedtools intersect -a $MUT_UNIQ -b $TFBS_CENTER -wa -wb |
#   cut -f1-2,4,6,8 |
#   awk '{print $1"\t"$2"\t"$2"\t"$3"\t"($2-$4-1000)"\t"$5}' > $MUT_INTR
# cut -f1,5 $MUT_INTR |
#   awk '{print $1"\t"$2"\t"$2}' > $MUT_CENT

# python test-mut-profile.py
