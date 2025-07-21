#!/bin/bash
#SBATCH --mem 123G
#SBATCH --partition=jic-short,nbi-short
#SBATCH --time=0-01:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=../logs/%x_%A_%a.out
#SBATCH --error=../logs/%x_%A_%a.err
#SBATCH --array=0-15

date
echo

##############################################################################################################################
echo "Setting directories and loading packages"


### Load bedtools
source package b0ed0698-358b-4c9b-9d21-603ea8d6e478

### Set wd
cd ../

### Define directory containing coverage files
cov_dir=inputs/scaled_read_coverage


##############################################################################################################################
echo "Defining which coverage file to use for this array task"


### Create list of coverage files
readarray -t file_list < <(ls $cov_dir/*.bedGraph)

### Create shortcut var. for SLURM task ID
i=$SLURM_ARRAY_TASK_ID

### Save coverage file for this array task to a variable
cov_file=$( echo ${file_list[i]} | awk -F/ '{print $NF}' )

### Save timepoint/rep name for this array task to a variable
time_rep_combo=$( echo ${file_list[i]} | awk -F/ '{print $NF}' | cut -d "_" -f 1)


##############################################################################################################################
echo "Running bedtools and awk operations:"


echo "Running bedtools intersect on consensus peaks and selected coverage file"
bedtools intersect -a inputs/consensus_peaks_min_cols.bed -b $cov_dir/$cov_file -wa -wb > outputs/stagewise_max_cov_files/peak_cov_intersect_${time_rep_combo}.bed

echo "Finding and retaining interval with highest coverage for each peak"
sort -k1,1 -k2,2n -k7,7nr outputs/stagewise_max_cov_files/peak_cov_intersect_${time_rep_combo}.bed | \
    awk '!seen[$1,$2,$3]++' | \
    awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > outputs/stagewise_max_cov_files/max_cov_per_peak_${time_rep_combo}.bed
wait

echo "Deleting intermediate peak_cov_intersect file"
rm outputs/stagewise_max_cov_files/peak_cov_intersect_${time_rep_combo}.bed


echo "Finished all processes"

echo
date