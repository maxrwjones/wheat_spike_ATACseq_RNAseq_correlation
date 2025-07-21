#!/bin/bash
#SBATCH --mem 10G
#SBATCH --partition=jic-short,nbi-short
#SBATCH --time=1-00:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=../logs/%x_%A.out
#SBATCH --error=../logs/%x_%A.err

### Set wd
cd ../

### Load bedtools
source package b0ed0698-358b-4c9b-9d21-603ea8d6e478

### Create an array of file names in target folder
readarray -t file_list < <(ls outputs/single_gene_BEDs/)

### Count number of genes in single

### Loop through all genes
for i in "${file_list[@]}"
do
    ### Echo loop number for tracking purposes
    echo $i

    ### Strip off the ".txt" from current target file and store in a variable
    gene_ID=$( echo $i | sed -e "s/\.txt$//")

    ### Run bedtools intersect
    bedtools intersect -wb -sorted -a outputs/single_gene_BEDs/$i -b outputs/consensus_peaks_with_max_cov.bed >> \
    outputs/consensus_peaks_with_parent_genes.bed

done

echo
date