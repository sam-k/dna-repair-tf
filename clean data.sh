#!/bin/bash

DATAPATH=../datasets
MUT_FILE=$DATAPATH/simple_somatic_mutation.open.SKCM-US.tsv

# Sample dataset to work with.
MUT_SAMPLE=$DATAPATH/sample_ssm_SKCM-US.tsv
head -n 1000 $MUT_FILE > $MUT_SAMPLE

# Select and sort columns.
MUT_SAMPLE_CUT=$DATAPATH/mut_only_sample_ssm_SKCM-US.tsv
cut -f9-11,16,17 $MUT_SAMPLE > $MUT_SAMPLE_CUT
MUT_SAMPLE_SORT=$DATAPATH/sorted_sample_ssm_SKCM-US.tsv
(head -n 1 $MUT_SAMPLE_CUT && tail -n +2 $MUT_SAMPLE_CUT | sort -V) > $MUT_SAMPLE_SORT

# Perform operations on actual dataset.
#MUT_FILE_CUT=$DATAPATH/mut_only_ssm_SKCM-US.tsv
#cut -f9-11,16,17 $MUT_FILE | sort > $MUT_FILE_CUT
