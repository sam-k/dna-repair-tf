#!/usr/bin/env bash
#SBATCH --job-name all-mut-profs_helper
#SBATCH --mail-user sdk18@duke.edu
#SBATCH --mail-type END,FAIL
#SBATCH --time 12:00:00
#SBATCH -c 1
#SBATCH --output logs/all-mut-profiles_helper.cl.out
#SBATCH --error logs/all-mut-profiles_helper.cl.err

module load python
module load Anaconda

FILENAME="$1"
MUT="$2"
TFBS="$3"
RUN_ID="$4"
WHICH_DATA="$5"
PACKAGE="$6"
_GENERATE_PROFILES="$7"
_GENERATE_FIGURES="$8"
_BENCHMARK="$9"

# Queue bash script
if [[ $_GENERATE_PROFILES -eq 0 ]]; then
  sbatch -W "$FILENAME" "$MUT" "$TFBS" "$RUN_ID" "$PACKAGE" $_BENCHMARK
fi

# Run python script
if [[ $_GENERATE_FIGURES -eq 0 ]]; then
  python "./mut-profile.py" "$RUN_ID" "$WHICH_DATA"
fi
