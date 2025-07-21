#!/bin/bash
#SBATCH --mem 10G
#SBATCH --partition=jic-short,nbi-short
#SBATCH --time=0-00:30
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=../logs/%x_%A_%a.out
#SBATCH --error=../logs/%x_%A_%a.err
#SBATCH --array=0-6

### Set wd
cd ../

### Make output directory
mkdir -p outputs/streamlined_MACS2_peaks

### Set 'i' to array task number
i=$SLURM_ARRAY_TASK_ID

### Create an array of file names in the MACS2_peaks input directory
readarray -t file_list < <(ls inputs/MACS2_peaks/*.narrowPeak)

### Strip off the ".narrowPeak" from current target file and store in a variable
file_ID=$( echo ${file_list[i]} | sed -e "s/\.narrowPeak$//" | sed 's#inputs/MACS2_peaks/##g')

### Select fields of interest and store in a new file in the sub-directory "streamlined_peaks"
awk 'BEGIN {FS = "\t"; OFS = "\t"} {print $1, $2, $3, $4, $7, $9}' ${file_list[i]} \
> outputs/streamlined_MACS2_peaks/${file_ID}_streamlined.txt