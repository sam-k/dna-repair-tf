#!/bin/bash

### Run locally on small sample of data

MUT_DATASET=SKCM-US
TFBS_DATASET=skcm

MUT_FILE="../datasets/simple_somatic_mutation.open.${MUT_DATASET}.tsv"
TFBS_FILE="../datasets/proximalTFBS-DHS_${TFBS_DATASET}.bed"
# TFBS_DHS_FILE="../datasets/proximalTFBS-noDHS_${TFBS_DATASET}.bed"

## MUT_FILE:
#  Mutation locations on patient genomes
#  1. icgc_mutation_id
#  2. icgc_donor_id
#  3. project_code
#  4. icgc_specimen_id
#  5. icgc_sample_id
#  6. matched_icgc_sample_id
#  7. submitted_sample_id
#  8. submitted_matched_sample_id
#  9. chromosome
# 10. chromosome_start
# 11. chromosome_end
# 12. chromosome_strand
# 13. assembly_version
# 14. mutation_type
# 15. reference_genome_allele
# 16. mutated_from_allele
# 17. mutated_to_allele
# 18. quality_score
# 19. probability
# 20. total_read_count
# 21. mutant_allele_read_count
# 22. verification_status
# 23. verification_platform
# 24. biological_validation_status
# 25. biological_validation_platform
# 26. consequence_type
# 27. aa_mutation
# 28. cds_mutation
# 29. gene_affected
# 30. transcript_affected
# 31. gene_build_version
# 32. platform
# 33. experimental_protocol
# 34. sequencing_strategy
# 35. base_calling_algorithm
# 36. alignment_algorithm
# 37. variation_calling_algorithm
# 38. other_analysis_algorithm
# 39. seq_coverage
# 40. raw_data_repository	
# 41. raw_data_accession
# 42. initial_data_release_date

## TFBS_FILE:
#  Transcription factor-binding sites on genome
#  1. chromosome
#  2. chromosome_start
#  3. chromosome_end
#  4. transcription_factor

# Transform TFBSs into TFBS centers ±1000 bp.
TFBS_CNTR="./data/proximalTFBS-DHS_${TFBS_DATASET}_center1000.bed"
awk '{center=int(($2+$3)/2); print $1"\t"(center-1000)"\t"(center+1000)"\t"$4}' $TFBS_FILE |
  sort -V > $TFBS_CNTR

## TFBS_CNTR:
#  Region of ±1000 bp around center of each TFBS
#  1. chromosome
#  2. region_start_neg1000
#  3. region_end_pos1000
#  4. transcription_factor

MUT_FILE_PREFIX="./data/ssm.open.${MUT_DATASET}"
MUT_PREP="${MUT_FILE_PREFIX}_prepped.bed"
cut -f9-11,16,17 $MUT_FILE | # select cols
  sed -e 1d | # remove header
  sort -V | # sort
  sed -e $'s/\t/>/4' | # preprocess to BED format
  sed -e 's/^/chr/' |
  uniq > $MUT_PREP # remove duplicates

## MUT_PREP:
#  Mut locations preprocessed to fit BED format
#  1. chromosome
#  2. mutation_location
#  3. mutation_location
#  4. mutation_allele

MUT_INTR="${MUT_FILE_PREFIX}_intersect.bed"
bedtools intersect -a $MUT_PREP -b $TFBS_CNTR -wa -wb | # intersect with TFBS ±1000bp regions
  cut -f1-4,6-8  > $MUT_INTR

## MUT_INTR:
#  Preprocessed mut locations that are in TFBS ±1000bp regions
#  1. chromosome
#  2. mutation_location
#  3. mutation_location
#  4. mutation_allele
#  5. region_start_neg1000
#  6. region_end_pos1000
#  7. transcription_factor

MUT_CNTR="${MUT_FILE_PREFIX}_centered.bed"
cut -f1,5 $MUT_INTR |
  awk '{print $1"\t"$2"\t"$2}' > $MUT_CNTR







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
