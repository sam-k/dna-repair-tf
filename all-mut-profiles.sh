#!/bin/bash

### Calls mut-profile.sh on all datasets.
### Run on all datasets using cluster.

declare -a mut=(
  "BLCA"  "BRCA"  "COAD"  "COCA"  "HNSC"  "LUAD"  "LUSC"
  "MELA"  "SKCA"  "SKCM"
)
declare -a tfbs=(
  "blca"  "brca"  "crc"   "crc"   "hnsc"  "luad_lusc" "luad_lusc"
  "skcm"  "skcm"  "skcm"
)

for ((i=0;i<${#mut[@]};++i)); do
    sbatch ./mut-profile.cl.sh "${mut[i]}" "${tfbs[i]}" "noDHS"
done