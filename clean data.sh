#!/bin/bash

DATAPATH=../datasets
MUT_FILE=$DATAPATH/simple_somatic_mutation.open.SKCM-US.tsv
TFBS_FILE=$DATAPATH/proximalTFBS-DHS_skcm.bed

# Transform TFBSs into TFBS centers ±1000 bp.
TFBS_CENTER=data/proximalTFBS-DHS_skcm_center1000.bed
awk '{center=int(($2+$3)/2); $2=center-1000; $3=center+1000; print $1"\t"$2"\t"$3"\t"$4}' $TFBS_FILE |
  sort -V > $TFBS_CENTER

# Sample dataset to work with.
MUT_SAMPLE=data/samples/sample.ssm.SKCM-US.tsv
gshuf -n 1000 $MUT_FILE > $MUT_SAMPLE

# Select, sort and preprocess columns.
MUT_SAMPLE_CUT=data/samples/sample.ssm.SKCM-US_mut_only.tsv #select
cut -f9-11,16,17 $MUT_SAMPLE > $MUT_SAMPLE_CUT
MUT_SAMPLE_SORT=data/samples/sample.ssm.SKCM-US_sorted.tsv #sort
(head -n 1 $MUT_SAMPLE_CUT && tail -n +2 $MUT_SAMPLE_CUT | sort -V) > $MUT_SAMPLE_SORT
MUT_SAMPLE_BED=data/samples/sample.ssm.SKCM-US_prepped.bed #preprocess
tail -n +2 $MUT_SAMPLE_SORT |
  sed -e $'s/\t/>/4' |
  sed -e 's/^/chr/' > $MUT_SAMPLE_BED

## Intersect with TFBS centers.
MUT_SAMPLE_INT=data/samples/sample.ssm.SKCM-US_intersect.bed #intersect
bedtools intersect -a $MUT_SAMPLE_BED -b $TFBS_CENTER -wa > $MUT_SAMPLE_INT

# Perform operations on actual dataset.
#MUT_PREP=$DATAPATH/ssm.open.SKCM-US_prepped.bed
#MUT_INT=$DATAPATH/ssm.open.SKCM-US_intersect.bed
#cut -f9-11,16,17 $MUT_FILE |
#  tail -n +2 |
#  sort -V |
#  sed -e $'s/\t/>/4' |
#  sed -e 's/^/chr/' > $MUT_PREP
#bedtools intersect -a $MUT_PREP -b $TFBS_CENTER -wa > $MUT_INT