#!/bin/bash

### Run snippets for testing purposes

MUT_DATASET=SKCM-US
TFBS_DATASET=skcm

MUT_FILE="./data/sample.ssm.open.${MUT_DATASET}.tsv"
MUT_FILE_PREFIX="./data/sample.ssm.open.${MUT_DATASET}"
MUT_PREP="${MUT_FILE_PREFIX}_prepped.bed"
MUT_INTR="${MUT_FILE_PREFIX}_intersect.bed"
MUT_CNTR="${MUT_FILE_PREFIX}_centered.bed"

TFBS_FILE="../datasets/proximalTFBS-DHS_${TFBS_DATASET}.bed"
TFBS_CNTR="./data/sample.proximalTFBS-DHS_${TFBS_DATASET}_center1000.bed"

cut -f1,6-7 $MUT_INTR |
  awk '{print $1"\t"$2"\t"$2"\t"$3}' > $MUT_CNTR