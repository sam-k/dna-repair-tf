#!/usr/bin/env bash
#SBATCH --job-name mut-prof_merged
#SBATCH --mail-user sdk18@duke.edu
#SBATCH --mail-type END,FAIL
#SBATCH --time 12:00:00
#SBATCH -c 4
#SBATCH --output logs/mut-profile_merged.out.txt
#SBATCH --error logs/mut-profile_merged.err.txt

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

MUT_FILE="../datasets/simple_somatic_mutation.open.${MUT_DATASET}.tsv"
MERGED_TFBS_FILE="../datasets/merged_ENCODE.tf.bound.union.bed"
DHS_FILE="../datasets/${DHS_ID}-DNase.hotspot.fdr0.01.peaks.v2.bed"
TSS_FILE="../datasets/refseq_TSS_hg19_170929.bed"
ENH_FILE="../datasets/permissive_enhancers.bed"
GEN_FILE="../datasets/bedtools_hg19_sorted.txt"

MUT_CNTR="./data/ssm.open.${RUN_ID}_${MUT_DATASET}_centered.bed"
MUT_CNTR_PRO="./data/ssm.open.${RUN_ID}_${MUT_DATASET}_pro_centered.bed"
MUT_CNTR_ENH="./data/ssm.open.${RUN_ID}_${MUT_DATASET}_enh_centered.bed"

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
DHS_SORTED="./data/supplementary/${DHS_ID}-DNase.hotspot.fdr0.01.peaks.v2_sorted.bed"
sort -V "$DHS_FILE" |  # sort
  uniq > "$DHS_SORTED"  # remove duplicates

# Build active TFBSs from merged TFBSs and cancer-specific DHSs.
TFBS_FILE="./data/supplementary/activeTFBS_${TFBS_DATASET}.bed"
sort -V "$MERGED_TFBS_FILE" |  # sort
  uniq |  # remove duplicates
  if [[ "$PACKAGE" == "bedtools" ]]; then
    bedtools intersect -a - -b "$DHS_SORTED" -wa -sorted -g "$GEN_FILE"  # intersect with cacner-specific DHSs
  elif [[ "$PACKAGE" == "bedops" ]]; then
    exit 1 # not yet implemented
  fi |
  sort -V > "$TFBS_FILE"

# Transform TFBSs into TFBS centers ±1000 bp.
TFBS_CNTR="./data/supplementary/activeTFBS_${TFBS_DATASET}_center1000.bed"
awk '{center=int(($2+$3)/2); print $1"\t"(center-1000)"\t"(center+1000)"\t"$4}' "$TFBS_FILE" |
  sort -V > "$TFBS_CNTR"

## TFBS_CNTR:
#  Region of ±1000 bp around center of each TFBS
#  1. chromosome
#  2. region_start_neg1000
#  3. region_end_pos1000
#  4. transcription_factor

# Transform TSSs into their upstream regions.
TSS_REG="./data/supplementary/refseq_TSS_up${UPSTREAM}-down${DOWNSTREAM}.bed"
cut -f1-2 "$TSS_FILE" |  # select cols
 sort -V |  # sort
 awk '{print $1"\t"($2-"'$UPSTREAM'")"\t"($2+"'$DOWNSTREAM'")}' > "$TSS_REG"

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

# Intersect with mutation data.
MUT_INTR="./data/supplementary/ssm.open.${RUN_ID}_${MUT_DATASET}_intersected.bed"
cut -f9-11,16,17 "$MUT_FILE" |  # select cols
  sed -e 1d |  # remove header
  sed -e $'s/\t/>/4' |  # preprocess to BED format
  sed -e 's/^/chr/' |
  sort -V |
  uniq | # remove duplicates
  if [[ "$PACKAGE" == "bedtools" ]]; then
    bedtools intersect -a - -b "$TFBS_CNTR" -wa -wb -sorted -g "$GEN_FILE"  # intersect with TFBS ±1000bp regions
  elif [[ "$PACKAGE" == "bedops" ]]; then
    exit 1  # not yet implemented
  fi > "$MUT_INTR"

## MUT_INTR:
#  Mutations in ±1000bp active TFBS regions
#  1. mutation_chromosome
#  2. mutation_location
#  3. mutation_location
#  4. mutation_allele
#  5. TFBS_chromosome
#  6. TFBS_region_start_neg1000
#  7. TFBS_region_end_pos1000
#  8. transcription_factor

# Reexpress mut locations as distances from centers of ±1000bp TFBS regions
cut -f1-2,4,6,8 "$MUT_INTR" |
  awk '{dist=$2-$4-1000; print $1"\t"dist"\t"dist"\t"$3"\t"$5}' |
  sort -V > "$MUT_CNTR"

## MUT_CNTR:
#  Mut locations as distances from centers of ±1000bp TFBS regions
#  1. chromosome
#  2. mutation_distance_from_center
#  3. mutation_distance_from_center
#  4. mutation_allele
#  5. transcription_factor

# Intersect further with promoters or enhancers.
intersect_further() {
  in_file="$1"
  out_file="$2"

  if [[ "$PACKAGE" == "bedtools" ]]; then
    bedtools intersect -a "$MUT_INTR" -b "$in_file" -wa -sorted -g "$GEN_FILE"
  elif [[ "$PACKAGE" == "bedops" ]]; then
    exit 1  # not yet implemented
  fi |
    cut -f1-2,4,6,8 |
    awk '{dist=$2-$4-1000; print $1"\t"dist"\t"dist"\t"$3"\t"$5}' |
    sort -V |
    uniq > "$out_file"
}
intersect_further "$TSS_REG" "$MUT_CNTR_PRO"
intersect_further "$ENH_PROC" "$MUT_CNTR_ENH"

# Benchmark end, in ms.
if [[ $_BENCHMARK -eq 0 ]]; then
  end_time=`python -c "from time import time; print(int(time()*1000))"`
  echo "${RUN_ID}_${MUT_DATASET}" >> "$BENCHMARK_FILE"
  echo "$((end_time-start_time)) ms" >> "$BENCHMARK_FILE"  # duration
  echo >> "$BENCHMARK_FILE"  # newline
fi
