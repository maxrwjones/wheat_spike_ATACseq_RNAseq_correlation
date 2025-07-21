#!/bin/bash
#SBATCH --mem 32G
#SBATCH --partition=jic-short,nbi-short
#SBATCH --time=0-02:00
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=../logs/%x_%A_%a.out
#SBATCH --error=../logs/%x_%A_%a.err
#SBATCH --array=1-100

source /jic/software/staging/RCSUPPORT-2912/stagingloader

### Change working directory
cd ../outputs

i=$SLURM_ARRAY_TASK_ID

Rscript simulate_enhancers.r $i