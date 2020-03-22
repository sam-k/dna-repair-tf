#!/usr/bin/env bash
#SBATCH --job-name all-mut-profs
#SBATCH --mail-user sdk18@duke.edu
#SBATCH --mail-type END,FAIL
#SBATCH --time 12:00:00
#SBATCH -c 3
#SBATCH --output logs/all-mut-profiles.cl.out
#SBATCH --error logs/all-mut-profiles.cl.err

### Calls mut-profile_TYPE.cl.sh on all datasets, as specified.

module load python
module load Anaconda

RUN_TYPE="enhancers"  # run type

TFBS_TYPE="distal"  # proximal, distal
TFBS_DHS="DHS"  # DHS, noDHS
CDS_FILE_ID=""  # coding regions file ID

WHICH_DATA="skcm"  # data group name
PACKAGE="bedtools"  # package to use

# Debug flags: 0 for true, 1 for false
_GENERATE_PROFILES=0
_GENERATE_FIGURES=1
_BENCHMARK=1

FILENAME="./mut-profile_${RUN_TYPE}.cl.sh"

## RUN_TYPE:
#  Bash script filename to be called.
#  all
#  enhancers
#  noncoding
#  tss
#  wgs

## TFBS_DHS:
#  Whether to use active (DHS) or inactive (noDHS) TFBSs.
#  DHS
#  noDHS

## TFBS_TYPE:
#  Whether to use proximal or distal (melanoma only) TFBSs.
#  proximal
#  distal

## WHICH_DATA:
#  Which group of somatic mutation data to use.
#  all
#  small
#  skcm

## CDS_FILE_ID:
#  Which coding regions file to use.
#  1: coding_exons.bed
#  2: cds.regions

## PACKAGE:
#  Which bioinformatics package to use.
#  bedops
#  bedtools



## Check arguments before proceeding any further.

invalid_arg_flag=1

check_args() {
  local type="$1"
  local arg="$2"
  local arr="$3[@]"
  for x in "${!arr}"; do
    if [[ $x == "$arg" ]]; then
      return
    fi
  done
  echo "Invalid ${type} argument: ${arg}"
  invalid_arg_flag=0
}

declare -a args=("" "enhancers" "noncoding" "tss" "wgs")
check_args "RUN_TYPE" "$RUN_TYPE" args

declare -a args=("DHS" "noDHS")
check_args "TFBS_DHS" "$TFBS_DHS" args

declare -a args=("distal" "proximal")
check_args "TFBS_TYPE" "$TFBS_TYPE" args

declare -a args=("" "1" "2")
check_args "CDS_FILE_ID" "$CDS_FILE_ID" args

declare -a args=("bedops" "bedtools")
check_args "PACKAGE" "$PACKAGE" args

# Select datasets to use
case "$WHICH_DATA" in
  # All mutation/TFBS datasets
  all )
    declare -a mut=(
      "BLCA"  "BRCA"  "COAD"  "COCA"  "HNSC"  "LUAD"  "LUSC"
      "MELA"  "READ"  "SKCA"  "SKCM"
    )
    declare -a tfbs=(
      "blca"  "brca"  "crc"   "crc"   "hnsc"  "luad_lusc" "luad_lusc"
      "skcm"  "crc"   "skcm"  "skcm"
    )
    ;;
  # Only small (<1 GB) mutation datasets, and their TFBSs
  small )
    declare -a mut=(
      "BLCA"  "COAD"  "HNSC"  "LUAD" "READ"
    )
    declare -a tfbs=(
      "blca"  "crc"   "hnsc"  "luad_lusc" "crc"
    )
    ;;
  # Only skin cancers
  skcm )
    declare -a mut=(
      "MELA"  "SKCA"  "SKCM"
    )
    declare -a tfbs=(
      "skcm"  "skcm"  "skcm"
    )
    ;;
  # Anything else
  *)
    echo "Invalid WHICH_DATA argument: ${WHICH_DATA}"
    invalid_arg_flag=0
    ;;
esac



## Queue scripts on cluster.

# If any argument was invalid, then quit
if [[ $invalid_arg_flag -eq 0 ]]; then
  exit 1
fi

# Build identifier: e.g., proximal-DHS_all
declare -A run_codes=(
  ["all"]="all"
  ["enhancers"]="enh"
  ["noncoding"]="NC"
  ["tss"]="TSS"
  ["wgs"]="WGS"
)
run_id="${TFBS_TYPE}-${TFBS_DHS}_${run_codes[${RUN_TYPE}]}"
if [[ "$RUN_TYPE" == "noncoding" ]]; then
  run_id="${run_id}${CDS_FILE_ID}"
fi

# Queue bash scripts
if [[ $_GENERATE_PROFILES -eq 0 ]]; then
  for ((i=0; i<${#mut[@]}; i++)); do
    sbatch "$FILENAME" "${mut[i]}" "${tfbs[i]}" "$run_id" "$PACKAGE" "$_BENCHMARK"
  done
fi



## Generate figures.

# Queue python script
if [[ $_GENERATE_FIGURES -eq 0 ]]; then
  python "./mut-profile.py" "$run_id" "$WHICH_DATA"
fi
