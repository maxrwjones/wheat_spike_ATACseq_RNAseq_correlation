#!/bin/bash
#SBATCH --mem 64G
#SBATCH --partition=jic-medium,nbi-medium
#SBATCH --time=0-18:00
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=../logs/%x_%A.out
#SBATCH --error=../logs/%x_%A.err

### Set wd
cd ../outputs

### Load R
source package 713701eb-ca16-4a79-bf10-7a21b51f274d

### Run Rscript which sorts out ImpulseDE2 objects and runs the programme
Rscript 10_wrangle_and_run_transcripts_ImpulseDE2.r