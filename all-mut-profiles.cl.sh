#!/bin/bash
#SBATCH --job-name mut-profile-all
#SBATCH --mail-user sdk18@duke.edu
#SBATCH --mail-type END,FAIL
#SBATCH --time 12:00:00
#SBATCH -c 2
#SBATCH --output logs/mut-profile-all-DHS.cl.out
#SBATCH --error logs/mut-profile-all-DHS.cl.err

### Run on all datasets using cluster

declare -a mut=("BLCA-US" "COCA-CN" "LUAD-US" "LUSC-US" "MELA-AU" "SKCM-US")
declare -a tfbs=("blca" "crc" "luad_lusc" "luad_lusc" "skcm" "skcm")

for ((i=0;i<${#mut[@]};++i)); do
    sh ./mut-profile.cl.sh "${mut[i]}" "${tfbs[i]}"
done