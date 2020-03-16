#!/bin/bash
#SBATCH --job-name mut-prof_NC
#SBATCH --mail-user sdk18@duke.edu
#SBATCH --mail-type END,FAIL
#SBATCH --time 12:00:00
#SBATCH -c 4
#SBATCH --output logs/mut-profile_NC.cl.out
#SBATCH --error logs/mut-profile_NC.cl.err

## Intersects somatic mutation coords w/ TFBS coords.
## Mutation file is filtered by noncoding coords.
#  Does not produce intermediate files.
#  Run on full data using cluster.

module load bedtools2

MUT_DATASET="$1"
TFBS_DATASET="$2"
TFBS_DHS="$3"
TFBS_TYPE="$4"

MUT_FILE="../datasets/simple_somatic_mutation.open.${MUT_DATASET}.tsv"
CDS_FILE="../datasets/coding_exons.bed" # NC1
# CDS_FILE="../datasets/cds.regions" # NC2
GEN_FILE="../datasets/human.hg38.genome"
TFBS_FILE="../datasets/${TFBS_TYPE}TFBS-${TFBS_DHS}_${TFBS_DATASET}.bed"

case "${TFBS_TYPE}" in
  distal )
    temp="distalTFBS_"
    ;;
  * )
    temp=""
    ;;
esac
MUT_CNTR="./data/ssm.open.${temp}NC1_${TFBS_DHS}_${MUT_DATASET}_centered.bed"

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

## CDS_FILE:
#  Coding regions on genome
#  1. chromosome
#  2. chromosome_start
#  3. chromosome_end
#  4. name
#  5. score
#  6. strand
#  NC1 : coding_exons.bed
#  NC2 : cds.regions

## GEN_FILE:
#  Lengths of each chromosome (1-26/M/X/Y) in genome
#  1. chromosome
#  2. length

## TFBS_FILE:
#  Transcription factor-binding sites on genome
#  1. chromosome
#  2. chromosome_start
#  3. chromosome_end
#  4. transcription_factor

# Complement coding regions to get noncoding regions.
NONCODING="./data/supplementary/cds_regions_complement.bed"
grep -P '^chr(\d+|[MXY])\t' "${CDS_FILE}" | # remove alt chr coords
  sort -V |
  bedtools complement -i - -g "${GEN_FILE}" > "${NONCODING}"

# Transform TFBSs into TFBS centers ±1000 bp.
TFBS_CNTR="./data/supplementary/distalTFBS-${DHS}_${TFBS_DATASET}_center1000.bed"
awk '{center=int(($2+$3)/2); print $1"\t"(center-1000)"\t"(center+1000)"\t"$4}' "${TFBS_FILE}" |
  sort -V > "${TFBS_CNTR}"

## TFBS_CNTR:
#  Region of ±1000 bp around center of each TFBS
#  1. chromosome
#  2. region_start_neg1000
#  3. region_end_pos1000
#  4. transcription_factor

cut -f9-11,16-17 "${MUT_FILE}" | # select cols
  sort -V | # sort
  sed -e $'s/\t/>/4' | # preprocess to BED format
  sed -e 's/^/chr/' |
  uniq | # remove duplicates
  bedtools intersect -a - -b "${NONCODING}" -wa -sorted | # intersect with noncoding regions
  bedtools intersect -a - -b "${TFBS_CNTR}" -wa -wb -sorted | # intersect with TFBS ±1000bp regions
  cut -f1-2,4,6,8 |
  awk '{dist=$2-$4-1000; print $1"\t"dist"\t"dist"\t"$3"\t"$5}' |
  sort -V > "${MUT_CNTR}"

## MUT_CNTR:
#  Mut locations as distances from centers of ±1000bp TFBS regions
#  1. chromosome
#  2. mutation_distance_from_center
#  3. mutation_distance_from_center
#  4. mutation_allele
#  5. transcription_factor