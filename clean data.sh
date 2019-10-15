#!/bin/bash

DATAPATH=../datasets
MUT_FILE=$DATAPATH/simple_somatic_mutation.open.SKCM-US.tsv
TFBS_FILE=$DATAPATH/proximalTFBS-DHS_skcm.bed

# Transform TFBSs into TFBS centers ±1000 bp.
TFBS_CENTER=$DATAPATH/proximalTFBS-DHS_skcm_center1000.bed
awk '{center=int(($2+$3)/2); $2=center-1000; $3=center+1000; print $1"\t"$2"\t"$3"\t"$4}' $TFBS_FILE |
  sort -V > $TFBS_CENTER

# Sample dataset to work with.
MUT_SAMPLE=$DATAPATH/sample_ssm_SKCM-US_0.tsv
gshuf -n 1000 $MUT_FILE > $MUT_SAMPLE

# Select, sort and preprocess columns.
MUT_SAMPLE_CUT=$DATAPATH/sample_ssm_SKCM-US_1_mut_only.tsv #select
cut -f9-11,16,17 $MUT_SAMPLE > $MUT_SAMPLE_CUT
MUT_SAMPLE_SORT=$DATAPATH/sample_ssm_SKCM-US_2_sorted.tsv #sort
(head -n 1 $MUT_SAMPLE_CUT && tail -n +2 $MUT_SAMPLE_CUT | sort -V) > $MUT_SAMPLE_SORT
MUT_SAMPLE_BED=$DATAPATH/sample_ssm_SKCM-US_3_prepped.bed #preprocess
tail -n +2 $MUT_SAMPLE_SORT |
  sed -e $'s/\t/>/4' |
  sed -e 's/^/chr/' > $MUT_SAMPLE_BED

# Intersect with TFBS centers.
MUT_SAMPLE_INT=$DATAPATH/sample_ssm_SKCM-US_4_intersect.bed #intersect
bedtools intersect -a $MUT_SAMPLE_BED -b $TFBS_CENTER -wa > $MUT_SAMPLE_INT

# Perform operations on actual dataset.
cut -f9-11,16,17 $MUT_FILE |
