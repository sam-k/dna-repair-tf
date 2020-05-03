#!/usr/bin/env bash
#SBATCH --job-name all-mut-profs
#SBATCH --mail-user sdk18@duke.edu
#SBATCH --mail-type FAIL
#SBATCH --time 12:00:00
#SBATCH -c 1
#SBATCH --output logs/OUT_all-mut-profiles.txt
#SBATCH --error logs/ERR_all-mut-profiles.txt

### Calls mut-profile_TYPE.cl.sh on all datasets, as specified.

RUN_TYPE="all"  # run type

TFBS_TYPE="proximal"  # proximal, distal
TFBS_DHS="DHS"  # DHS, noDHS
CDS_FILE_ID=""  # coding regions file ID

WHICH_DATA="MELA"  # data group name
PACKAGE="bedops"

# Debug flags: 0 for true, 1 for false
_GENERATE_PROFILES=0
_GENERATE_FIGURES=1
_BENCHMARK=0

FILENAME="./mut-profile_${RUN_TYPE}.cl.sh"

## RUN_TYPE:
#  Bash script filename to be called.
#  all
#  enhancers
#  merged (dhs, small_dhs)
#  merged-bg (dhs, small_dhs)
#  noncoding
#  tss
#  wgs

## TFBS_DHS:
#  Whether to use active or inactive TFBSs.
#  DHS
#  noDHS

## TFBS_TYPE:
#  Whether to use proximal or distal TFBSs.
#  proximal
#  distal (skcm)

## WHICH_DATA:
#  Which group of somatic mutation data to use.
#  all
#  dhs
#  small
#  small_dhs
#  skcm

## CDS_FILE_ID:
#  Which coding regions file to use.
#  (none)
#  1: coding_exons.bed
#  2: cds.regions

## PACKAGE:
#  Which bioinformatics package to use.
#  bedops
#  bedtools



## Check arguments before proceeding any further.

# Flag for whether any argument is invalid
invalid_arg_flag=1
check_args() {
  local type="$1"
  local arg="$2"
  local arr="$3[@]"
  local _flag="$4"
  for x in "${!arr}"; do
    if [[ $x == "$arg" ]]; then
      return
    fi
  done
  if [[ "$_flag" -eq 0 ]]; then
    echo "Invalid ${type} argument: ${arg}"
    invalid_arg_flag=0
  fi
}

declare -a args=("all" "enhancers" "merged" "merged-bg" "noncoding" "tss" "wgs")
check_args "RUN_TYPE" "$RUN_TYPE" args 0

declare -a args=("DHS" "noDHS")
check_args "TFBS_DHS" "$TFBS_DHS" args 0

declare -a args=("distal" "proximal")
check_args "TFBS_TYPE" "$TFBS_TYPE" args 0

declare -a args=("" "1" "2")
check_args "CDS_FILE_ID" "$CDS_FILE_ID" args 0

declare -a args=("bedops" "bedtools" "bedtools-sorted" "bedtools-unsorted")
check_args "PACKAGE" "$PACKAGE" args 0

# Select datasets to use
case "$WHICH_DATA" in
  # All cancers
  all )
    declare -a mut=(
      "BLCA" "BRCA" "COAD" "COCA" "HNSC" "LUAD" "LUSC"
      "MELA" "READ" "SKCA" "SKCM"
    )
    declare -a tfbs=(
      "blca" "brca" "crc"  "crc"  "hnsc" "luad_lusc" "luad_lusc"
      "skcm" "crc"  "skcm" "skcm"
    );;
  # Only cancers with DHS data
  dhs )
    declare -a mut=(
      "BRCA" "COAD" "COCA" "LUAD" "LUSC"
      "MELA" "READ" "SKCA" "SKCM"
    )
    declare -a tfbs=(
      "brca" "crc"  "crc"  "luad_lusc" "luad_lusc"
      "skcm" "crc"  "skcm" "skcm"
    );;
  # Only small (<1 GB) cancers
  small )
    declare -a mut=( "BLCA" "COAD" "HNSC" "LUAD" "READ")
    declare -a tfbs=("blca" "crc"  "hnsc" "luad_lusc" "crc");;
  # Only small (<1 GB) cancers with DHS data 
  small_dhs )
    declare -a mut=( "COAD" "LUAD" "READ")
    declare -a tfbs=("crc"  "luad_lusc" "crc");;
  # Only skin cancers
  skcm )
    declare -a mut=( "MELA" "SKCA" "SKCM")
    declare -a tfbs=("skcm" "skcm" "skcm");;
  # Individual cancer types
  BLCA ) declare -a mut=("BLCA"); declare -a tfbs=("blca");;
  BRCA ) declare -a mut=("BRCA"); declare -a tfbs=("brca");;
  COAD ) declare -a mut=("COAD"); declare -a tfbs=("crc");;
  COCA ) declare -a mut=("COCA"); declare -a tfbs=("crc");;
  READ ) declare -a mut=("READ"); declare -a tfbs=("crc");;
  HNSC ) declare -a mut=("HNSC"); declare -a tfbs=("hnsc");;
  LUAD ) declare -a mut=("LUAD"); declare -a tfbs=("luad_lusc");;
  LUSC ) declare -a mut=("LUSC"); declare -a tfbs=("luad_lusc");;
  MELA ) declare -a mut=("MELA"); declare -a tfbs=("skcm");;
  SKCA ) declare -a mut=("SKCA"); declare -a tfbs=("skcm");;
  SKCM ) declare -a mut=("SKCM"); declare -a tfbs=("skcm");;
  # Anything else
  *)
    echo "Invalid WHICH_DATA argument: ${WHICH_DATA}"
    invalid_arg_flag=0;;
esac

# If any argument was invalid, then quit
if [[ $invalid_arg_flag -eq 0 ]]; then
  exit 1
fi



## Queue scripts on cluster.

# Build identifier: e.g., proximal-DHS_WGS
declare -A run_codes=(
  ["all"]="all"
  ["enhancers"]="enh"
  ["merged"]="merged"
  ["merged-bg"]="mergedbg"
  ["noncoding"]="NC"
  ["tss"]="TSS"
  ["wgs"]="WGS"
)
run_id="${TFBS_TYPE}-${TFBS_DHS}_${run_codes[${RUN_TYPE}]}${CDS_FILE_ID}"

# Initialize benchmark file
if [[ $_BENCHMARK -eq 0 ]]; then
  BENCHMARK_FILE="./benchmark/${run_id}.txt"
  > "$BENCHMARK_FILE"  # clear file
fi

# Queue scripts
for ((i=0; i<${#mut[@]}; i++)); do
  sbatch "./all-mut-profiles_helper.cl.sh" "$FILENAME" "${mut[i]}" "${tfbs[i]}" "$run_id" "$WHICH_DATA" "$PACKAGE" $_GENERATE_PROFILES $_GENERATE_FIGURES $_BENCHMARK
done
