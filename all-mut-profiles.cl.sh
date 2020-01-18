#!/bin/bash
#SBATCH --job-name mut-profile-all
#SBATCH --mail-user sdk18@duke.edu
#SBATCH --mail-type END,FAIL
#SBATCH --time 12:00:00
#SBATCH -c 3
#SBATCH --output logs/mut-profile-all-DHS.cl.out
#SBATCH --error logs/mut-profile-all-DHS.cl.err

### Calls mut-profile.cl.sh.
### Run on all datasets using cluster.
# Does not work at the moment.

declare -a mut=(
  "BLCA-CN" "BLCA-US" "BRCA-EU" "BRCA-FR" "BRCA-KR" "BRCA-UK" "COAD-US"
  "COCA-CN" "HNSC-US" "LUAD-US" "LUSC-CN" "LUSC-KR" "LUSC-US" "MELA-AU"
  "SKCA-BR" "SKCM-US")
declare -a tfbs=(
  "blca"    "blca"    "brca"    "brca"    "brca"    "brca"    "crc"
  "crc"     "hnsc"    "luad_lusc" "luad_lusc" "luad_lusc" "luad_lusc" "skcm"
  "skcm"  "skcm")

for ((i=0;i<${#mut[@]};++i)); do
    echo "${mut[i]} ${tfbs[i]}"
done