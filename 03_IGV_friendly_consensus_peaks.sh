#!/bin/bash
#SBATCH --mem 32G
#SBATCH --partition=jic-short,nbi-short
#SBATCH --time=0-01:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=../logs/%x_%A.out
#SBATCH --error=../logs/%x_%A.err

### Set wd
cd ../outputs

### Remove BED fields which would prevent IGV working properly
awk 'BEGIN {FS = "\t"; OFS = "\t"} {print $1, $2, $3, $3 - $2}' consensus_peaks.bed > IGV_friendly_consensus_peaks.bed
