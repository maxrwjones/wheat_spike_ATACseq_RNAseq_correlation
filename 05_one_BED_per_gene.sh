#!/bin/bash
#SBATCH --mem 4G
#SBATCH --partition=jic-short,nbi-short
#SBATCH --time=0-01:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=../logs/%x_%A.out
#SBATCH --error=../logs/%x_%A.err


### Set wd
cd ../

### Make output directory for single gene BED files
mkdir outputs/single_gene_BEDs

input_file="outputs/HC_gene_coords_extended.bed"

while IFS=$'\t' read -r -a fields; do
    fourth_field="${fields[3]}"

    ### Remove newline characters from the fourth field
    fourth_field="${fourth_field%$'\n'}"
    fourth_field="${fourth_field%$'\r'}"
    output_file="${fourth_field}.txt"

    ### Use printf to output tab-delimited fields
    ### tr is then used to remove trailing carriage returns
    ### Finally, sed is used to remove tabs after the fourth field
    printf "%s\t" "${fields[@]}" | \
    tr -d '\r' | \
    sed 's/\t$//' > \
    "outputs/single_gene_BEDs/$output_file"
    
done < $input_file

#echo -e "${fields[@]}" >> "outputs/single_gene_BEDs/$output_file"