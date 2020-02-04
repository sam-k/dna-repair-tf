#!/bin/bash

### Calls mut-profile_TYPE.cl.sh on all datasets, as specified.

FILENAME="./mut-profile_noncoding.cl.sh"
DHS="DHS"

## FILENAME:
#  Bash script filename to be called.
#  ./mut-profile.cl.sh
#  ./mut-profile_wgs.cl.sh
#  ./mut-profile_noncoding.cl.sh

## DHS:
#  Whether to use DHS or noDHS TFBS file.
#  DHS
#  noDHS

# declare -a mut=(
#   "BLCA"  "BRCA"  "COAD"  "COCA"  "HNSC"  "LUAD"  "LUSC"
#   "MELA"  "READ"  "SKCA"  "SKCM"
# )
declare -a mut=(
  "MELA"  "SKCA"  "SKCM"
)
# declare -a tfbs=(
#   "blca"  "brca"  "crc"   "crc"   "hnsc"  "luad_lusc" "luad_lusc"
#   "skcm"  "crc"   "skcm"  "skcm"
# )
declare -a tfbs=(
  "skcm"
)

for ((i=0;i<${#mut[@]};++i)); do
    sbatch "${FILENAME}" "${mut[i]}" "${tfbs[i]}" "${DHS}"
done