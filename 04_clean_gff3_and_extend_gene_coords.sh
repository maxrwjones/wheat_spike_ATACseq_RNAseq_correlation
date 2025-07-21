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

date

### Change into working dir
cd ../

printf "Removing comment lines and retaining only gene-level entries \n"

grep -Fv "^#" inputs/IWGSC_v1.1_HC_20170706.sorted.gff3 | \
awk '$3 == "gene" { print $0 }' > outputs/HC_gene_coords.gff3

printf "Finished cleaning GFF3 file \n"
printf "Extending gene coordinates using custom R script"

### Run R script
Rscript 04_extend_gene_coords.r

printf "Done \n"


date
