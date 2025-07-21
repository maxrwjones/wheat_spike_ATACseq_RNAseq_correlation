#!/bin/bash
#SBATCH --mem 32G
#SBATCH --partition=jic-short,nbi-short
#SBATCH --time=0-02:00
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=../logs/%x_%A.out
#SBATCH --error=../logs/%x_%A.err

### Load R
source package 713701eb-ca16-4a79-bf10-7a21b51f274d

### Change working directory
cd ../outputs

### Run R script
Rscript normalise_peaks_and_genes.r