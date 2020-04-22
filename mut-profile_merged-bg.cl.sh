#!/usr/bin/env bash
#SBATCH --job-name mut-prof_merged
#SBATCH --mail-user sdk18@duke.edu
#SBATCH --mail-type END,FAIL
#SBATCH --time 12:00:00
#SBATCH -c 4
#SBATCH --output logs/OUT_mut-profile_merged-bg.txt
#SBATCH --error logs/ERR_mut-profile_merged-bg.txt

## Intersects somatic mutation coords w/ TFBS coords.
#  Does not produce intermediate files.
#  Run on full data using cluster.

module load bedtools2
module load bedops
module load python

# Definition of promoter regions
UPSTREAM="2000"
DOWNSTREAM="1000"

MUT_DATASET="$1"
TFBS_DATASET="$2"
RUN_ID="$3"
PACKAGE="$4"
_BENCHMARK="$5"

case "$TFBS_DATASET" in
  brca )      DHS_ID="E028";;
  crc )       DHS_ID="E084";;
  luad_lusc ) DHS_ID="E088";;
  skcm )      DHS_ID="E059";;
  * )
    echo "Unsupported TFBS_DATASET argument: ${TFBS_DATASET}"
    exit 1
    ;;
esac

# Import data
MUT_FILE="../datasets/ssm/simple_somatic_mutation.open.${MUT_DATASET}.tsv"
MERGED_TFBS_FILE="../datasets/tfbs/merged_ENCODE.tf.bound.union.bed"
DHS_FILE="../datasets/total_dhs/${DHS_ID}-DNase.hotspot.fdr0.01.peaks.v2.bed"
TSS_FILE="../datasets/refseq_TSS_hg19_170929.bed"
ENH_FILE="../datasets/permissive_enhancers.bed"
GEN_FILE="../datasets/hg19/bedtools_hg19_sorted.txt"
GEN_FA="../datasets/hg19/hg19.fa"

# Export data
BOUND_MUT_CNTR="./data/ssm.open.${RUN_ID}_${MUT_DATASET}_bound_centered.bed"
BOUND_MUT_CNTR_PRO="./data/ssm.open.${RUN_ID}_${MUT_DATASET}_pro_bound_centered.bed"
BOUND_MUT_CNTR_ENH="./data/ssm.open.${RUN_ID}_${MUT_DATASET}_enh_bound_centered.bed"
UNBOUND_MUT_CNTR="./data/ssm.open.${RUN_ID}_${MUT_DATASET}_unbound_centered.bed"
UNBOUND_MUT_CNTR_PRO="./data/ssm.open.${RUN_ID}_${MUT_DATASET}_pro_unbound_centered.bed"
UNBOUND_MUT_CNTR_ENH="./data/ssm.open.${RUN_ID}_${MUT_DATASET}_enh_unbound_centered.bed"
BENCHMARK_FILE="./benchmark/${RUN_ID}.txt"

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

## MERGED_TFBS_FILE:
#  Transcription factor-binding sites on genome
#  1. chromosome
#  2. chromosome_start
#  3. chromosome_end
#  4. transcription_factor

## DHS_FILE:
#  Cancer-specific DHS locations on genome
#  1. chromosome
#  2. chromosome_start
#  3. chromosome_end
#  4. ?
#  5. ?

## TSS_FILE:
#  Transcription start sites
#  1. chromosome
#  2. location
#  3. location

## ENH_FILE:
#  Enhancer regions
#  1. chromosome
#  2. chromosome_start
#  3. chromosome_end
#  4. chromosome:chromosome_start-chromosome_end
#  5. ?
#  6. ?
#  7. ?
#  8. ?



# Sort DHS_FILE lexicographically before intersecting.
TOTAL_DHS="./data/supplementary/${DHS_ID}-DNase.hotspot.fdr0.01.peaks.v2_sorted.bed"
sort -V "$DHS_FILE" |  # sort
  uniq > "$TOTAL_DHS"  # remove duplicates

# Build active TFBSs from merged TFBSs and cancer-specific DHSs.
BOUND_DHS="./data/supplementary/boundDHS_${TFBS_DATASET}.bed"
sort -V "$MERGED_TFBS_FILE" |  # sort
  uniq |  # remove duplicates
  if [[ "$PACKAGE" == "bedtools" ]]; then
    bedtools intersect -a - -b "$TOTAL_DHS" -wa -sorted -g "$GEN_FILE"  # intersect with cacner-specific DHSs
  elif [[ "$PACKAGE" == "bedops" ]]; then
    exit 1 # not yet implemented
  fi |
  sort -V |
  uniq > "$BOUND_DHS"
# BOUND_DHS_CNTR="$BOUND_DHS"
BOUND_DHS_CNTR="./data/supplementary/boundDHS_${TFBS_DATASET}_center1000.bed"
awk '{center=int(($2+$3)/2); print $1"\t"(center>=1000 ? center-1000 : 0)"\t"(center+1000)"\t"$4}' "$BOUND_DHS" |  # transform into centers ±1000 bp
  sort -V |
  uniq > "$BOUND_DHS_CNTR"

# Find parts of DHSs not found in bound DHSs. (background)
UNBOUND_DHS="./data/supplementary/unboundDHS_${TFBS_DATASET}.bed"
bedtools subtract -a "$TOTAL_DHS" -b "$BOUND_DHS" |
  cut -f1-4 |
  sort -V |
  uniq > "$UNBOUND_DHS"
# UNBOUND_DHS_CNTR="$UNBOUND_DHS"
UNBOUND_DHS_CNTR="./data/supplementary/unboundDHS_${TFBS_DATASET}_center1000.bed"
awk '{center=int(($2+$3)/2); print $1"\t"(center>=1000 ? center-1000 : 0)"\t"(center+1000)"\t"$4}' "$UNBOUND_DHS" |  # transform into centers ±1000 bp
  sort -V |
  uniq > "$UNBOUND_DHS_CNTR"
UNBOUND_FA="./data/supplementary/unboundDHS_${TFBS_DATASET}_center1000.fa"
bedtools getfasta -fi "$GEN_FA" -bed "$UNBOUND_DHS" -fo "$UNBOUND_FA"

## BOUND_DHS, UNBOUND_DHS:
#  CNTR: Region of ±1000 bp around center of each region
#  1. chromosome
#  2. region_start_neg1000
#  3. region_end_pos1000
#  4. transcription_factor



# Transform TSSs into their upstream regions.
TSS_REG="./data/supplementary/refseq_TSS_up${UPSTREAM}-down${DOWNSTREAM}.bed"
cut -f1-2 "$TSS_FILE" |  # select cols
 awk -v up=$UPSTREAM -v down=$DOWNSTREAM '{print $1"\t"($2>=up ? $2-up : 0)"\t"($2+down)}' |
 sort -V |  # sort
 uniq > "$TSS_REG"

## TSS_REG:
#  Assumed promoter regions: up/downstream region of each TSS
#  1. chromosome
#  2. location
#  3. location

# Process enhancers data.
ENH_PROC="./data/supplementary/permissive_enhancers_proc.bed"
cut -f1-3 "$ENH_FILE" |
  tail -n +2 |  # remove header
  sort -V > "$ENH_PROC"

## ENH_PROC:
#  Enhancer data without headers
#  1. chromosome
#  2. chromosome_start
#  3. chromosome_end



# Benchmark start, in ms.
if [[ $_BENCHMARK -eq 0 ]]; then
  start_time=`python -c "from time import time; print(int(time()*1000))"`
fi

# Prepare mutations data
MUT_PROC="./data/supplementary/simple_somatic_mutation.open.${MUT_DATASET}_proc.bed"
cut -f9-11,16,17 "$MUT_FILE" |  # select cols
  sed -e 1d |  # remove header
  sed -e $'s/\t/>/4' |  # preprocess to BED format
  sed -e 's/^/chr/' |
  sort -V |
  uniq > "$MUT_PROC"  # remove duplicates

# Intersect further with promoters or enhancers.
intersect_further() {
  _intr_file="$1"
  _reg_file="$2"
  _reg_cntr_file="$3"

  if [[ "$PACKAGE" == "bedtools" ]]; then
    bedtools intersect -a "$_intr_file" -b "$_reg_file" -wa -sorted -g "$GEN_FILE"
  elif [[ "$PACKAGE" == "bedops" ]]; then
    exit 1  # not yet implemented
  fi |
    cut -f1-2,4,6-8 |
    awk '{dist=$2-int(($4+$5)/2); print $1"\t"dist"\t"dist"\t"$3"\t"$6}' |
    sort -V |
    uniq > "$_reg_cntr_file"
}

# Intersect.
intersect() {
  in_file="$1"
  intr_file="$2"
  cntr_file="$3"
  pro_cntr_file="$4"
  enh_cntr_file="$5"

  # Intersect with mutation data.
  if [[ "$PACKAGE" == "bedtools" ]]; then
    bedtools intersect -a "$MUT_PROC" -b "$in_file" -wa -wb -sorted -g "$GEN_FILE"  # intersect with TFBS/non-TFBS regions
  elif [[ "$PACKAGE" == "bedops" ]]; then
    exit 1  # not yet implemented
  fi > "$intr_file"

  # Reexpress mut locations as distances from centers of TFBS/non-TFBS regions.
  cut -f1-2,4,6-8 "$intr_file" |
    awk '{dist=$2-int(($4+$5)/2); print $1"\t"dist"\t"dist"\t"$3"\t"$6}' |
    sort -V |
    uniq > "$cntr_file"
  
  intersect_further "$intr_file" "$TSS_REG" "$pro_cntr_file"
  intersect_further "$intr_file" "$ENH_PROC" "$enh_cntr_file"
}

BOUND_MUT_INTR="./data/supplementary/ssm.open.${RUN_ID}_${MUT_DATASET}_bound_intersected.bed"
intersect "$BOUND_DHS_CNTR" "$BOUND_MUT_INTR" "$BOUND_MUT_CNTR" "$BOUND_MUT_CNTR_PRO" "$BOUND_MUT_CNTR_ENH"

UNBOUND_MUT_INTR="./data/supplementary/ssm.open.${RUN_ID}_${MUT_DATASET}_unbound_intersected.bed"
intersect "$UNBOUND_DHS_CNTR" "$UNBOUND_MUT_INTR" "$UNBOUND_MUT_CNTR" "$UNBOUND_MUT_CNTR_PRO" "$UNBOUND_MUT_CNTR_ENH"

## *_MUT_INTR:
#  Mutations in TFBS/non-TFBS regions
#  1. mutation_chromosome       1
#  2. mutation_location         2
#  3. mutation_location
#  4. mutation_allele           3
#  5. TFBS_chromosome
#  6. TFBS_region_start_neg1000 4
#  7. TFBS_region_end_pos1000   5
#  8. transcription_factor      6

## *_MUT_CNTR:
#  Mut locations as distances from centers of TFBS/non-TFBS regions
#  1. chromosome
#  2. mutation_distance_from_center
#  3. mutation_distance_from_center
#  4. mutation_allele
#  5. transcription_factor

# Benchmark end, in ms.
if [[ $_BENCHMARK -eq 0 ]]; then
  end_time=`python -c "from time import time; print(int(time()*1000))"`
  echo "${RUN_ID}_${MUT_DATASET}" >> "$BENCHMARK_FILE"
  echo "$((end_time-start_time)) ms" >> "$BENCHMARK_FILE"  # duration
  echo >> "$BENCHMARK_FILE"  # newline
fi
