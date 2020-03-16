#!/bin/bash

module load python

### Calls mut-profile_TYPE.cl.sh on all datasets, as specified.

TYPE="enhancers" # run type
TFBS_DHS="DHS" # DHS, noDHS
TFBS_TYPE="distal" # proximal, distal
WHICH_DATA="skcm" # data group name (see below)

FILENAME="./mut-profile_${TYPE}.cl.sh"

## FILENAME:
#  Bash script filename to be called.
#  ./mut-profile.cl.sh
#  ./mut-profile_wgs.cl.sh
#  ./mut-profile_noncoding.cl.sh

## DHS:
#  Whether to use DHS or noDHS TFBS file.
#  DHS
#  noDHS

case "${WHICH_DATA}" in
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
    echo "Invalid argument"
    exit 1
    ;;
esac

# Queue scripts on cluster
for ((i=0; i<${#mut[@]}; i++)); do
    sbatch "${FILENAME}" "${mut[i]}" "${tfbs[i]}" "${TFBS_DHS}" "${TFBS_TYPE}"
done

# Generate figures
python "./mut-profile.cl.py" "${TYPE}" "${TFBS_DHS}" "${TFBS_TYPE}" "${WHICH_DATA}"