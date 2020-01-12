MUT_DATASET=SKCM-US
MUT_FILE="./data/sample.ssm.open.${MUT_DATASET}.tsv"
MUT_FILE_PREFIX="./data/sample.ssm.open.${MUT_DATASET}"
MUT_PREP="${MUT_FILE_PREFIX}_prepped.bed"
MUT_INTR="${MUT_FILE_PREFIX}_intersect.bed"

TFBS_DATASET=skcm
TFBS_FILE="../datasets/proximalTFBS-DHS_${TFBS_DATASET}.bed"
TFBS_CNTR="./data/sample.proximalTFBS-DHS_${TFBS_DATASET}_center1000.bed"

bedtools intersect -a $MUT_PREP -b $TFBS_CNTR -wa -wb > $MUT_INTR # intersect with TFBS centers
    # cut -f1-2,4,6,8 |
    # awk '{print $1"\t"$2"\t"$2"\t"$3"\t"($2-$4-1000)"\t"$5}' > $MUT_INTR