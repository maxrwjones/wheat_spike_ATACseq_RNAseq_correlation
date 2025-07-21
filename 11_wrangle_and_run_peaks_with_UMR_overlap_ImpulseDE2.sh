#!/bin/bash
#SBATCH --mem 248G
#SBATCH --partition=jic-medium,nbi-medium
#SBATCH --time=1-00:00
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=../logs/%x_%A.out
#SBATCH --error=../logs/%x_%A.err

### Set name of UMR file
UMR_file=inputs/CS_UMRs_cov5_unsplit_coords.bed

### Set wd
cd ../

### Create simplified consensus peak / parent gene BED file
awk -F'\t' 'BEGIN {OFS="\t"} {printf "%s\t%s\t%s\t%s", $5, $6, $7, $4; for (i=8; i<=NF; i++) printf "\t%s", $i; print ""}' \
outputs/consensus_peaks_with_parent_genes.bed \
> outputs/consensus_peaks_with_parent_genes_rejigged.bed

### Load bedtools
source package b0ed0698-358b-4c9b-9d21-603ea8d6e478



### Run bedtools to add a column tracking the number of UMRs each peak overlaps with
### (0 for none, or 1, 2, 3, etc)
bedtools intersect -c -a outputs/consensus_peaks_with_parent_genes_rejigged.bed -b $UMR_file \
> outputs/consensus_peaks_with_parent_genes_intersecting_UMRs.bed



### Load R
source package 713701eb-ca16-4a79-bf10-7a21b51f274d

### Set wd
cd outputs/

### Run Rscript which sorts out ImpulseDE2 objects and runs the programme
Rscript 11_wrangle_and_run_peaks_with_UMR_overlap_ImpulseDE2.r