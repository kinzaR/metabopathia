#!/bin/bash

#SBATCH --job-name=merge_SEs
#SBATCH --mail-type=END
#SBATCH --mail-user=rian.kinza@gmail.com
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --output=logs/outputs_merging_%A_%a.out
#SBATCH --error=logs/errors_merging_%A_%a.err

## wake conda up -_-
eval "$(conda shell.bash hook)"
conda activate r-environment

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "Script directory:"
echo "------------"
echo "$SCRIPTDIR"
echo "------------"
echo ""
echo "command:"
echo "------------"
echo "$0 $@"
echo "------------"
echo ""
echo "------------"
echo "Rscript path"
echo "------------"
#module load R/3.1.2
#module load R
#module load R/3.5.1
which Rscript
echo "------------"

Rscript merge_SEs.R
