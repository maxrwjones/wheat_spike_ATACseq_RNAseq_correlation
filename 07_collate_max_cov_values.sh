#!/bin/bash
#SBATCH --mem 32G
#SBATCH --partition=jic-medium,nbi-medium
#SBATCH --time=0-04:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=../logs/%x_%A.out
#SBATCH --error=../logs/%x_%A.err

### Load R
source package 713701eb-ca16-4a79-bf10-7a21b51f274d

### Set wd
cd ../

### Run R script
Rscript 07_collate_max_cov_values.r