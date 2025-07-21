#!/bin/bash
#SBATCH --mem 20G
#SBATCH --partition=jic-short,nbi-short
#SBATCH --time=0-01:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=../logs/%x_%A.out
#SBATCH --error=../logs/%x_%A.err

### Set wd
cd ../

### Load bedtools
source package b0ed0698-358b-4c9b-9d21-603ea8d6e478

### Count the number of MACS2 peak files that are going to be used
file_num=$( ls -1 streamlined_MACS2_peaks/*.txt | wc -l)

### Run bedtools multiinter 
bedtools multiinter -i streamlined_MACS2_peaks/*.txt \
| awk -v file_num="$file_num" '$4 == file_num { print $0 }' \
> consensus_peaks.bed

### Create a version of consensus peaks without the excess 'counter' columns from bedtools multiinter
cut -f1-3 consensus_peaks.bed > consensus_peaks_min_cols.bed


#########################################################################################
### Final file created below not strictly necessary; these are versions of the consensus
### peak BED with added statistics like fold-change versus background and q-value.
#########################################################################################

### Now intersect the new consensus peaks with all original peak files,
### but this time use "bedtools intersect" to report original B-track features.
### This means we get to keep the peak names, fold-change scores, and peak offsets.
### However, each peak file intersecting an A-track feature gets its own line - not ideal.
### For now, lets just store the numbers we are interested in; fold-change
### (field 17) and q-value (field 18).
bedtools intersect -wa -wb -a consensus_peaks.bed -b streamlined_MACS2_peaks/*.txt \
| awk 'BEGIN {FS = "\t"; OFS = "\t"} {print $18, $19}' \
> consensus_peak_stats.txt


### Next, we iterate over groups of lines equal in size to the number of peak files,
### taking the average of each column. This gives us average stats per consensus peak.
cat consensus_peak_stats.txt \
| awk -v file_num="$file_num" -F'\t' '{
    for (i = 1; i <= NF; i++) {
        sum[i] += $i
    }
    if (NR % file_num == 0) {
        for (i = 1; i <= NF; i++) {
            printf "%s\t", sum[i] / file_num
            sum[i] = 0
        }
        printf "\n"
    }
}' \
> consensus_peak_stats_avg.txt

### We then need to column-bind these stats about the consensus peaks to the
### consensus peak file itself.
paste -d'\t' consensus_peaks.bed consensus_peak_stats_avg.txt \
| sed 's/\t$//' \
> consensus_peaks_plus_stats.txt

### It would be great to also column-bind the stats from each individual peak
### file to the final consensus peaks file... these will be needed for plots and
### correlation of ATAC-seq and RNA-seq data later.
awk -v file_num="$file_num" -F'\t' \
'{ a[NR] = $1 } END { for (i = 1; i <= NR; i++) { printf "%s\t", a[i]; if (i % file_num == 0) printf "\n"; } printf "\n" }' consensus_peak_stats.txt \
| sed 's/\t$//' \
> per_timepoint_foldchanges.txt

awk -v file_num="$file_num" -F'\t' \
'{ a[NR] = $2 } END { for (i = 1; i <= NR; i++) { printf "%s\t", a[i]; if (i % file_num == 0) printf "\n"; } printf "\n" }' consensus_peak_stats.txt \
| sed 's/\t$//' \
> per_timepoint_qvalues.txt

### Create final consensus peaks file with stats and raw fold-changes and qvalues
paste -d'\t' consensus_peaks_plus_stats.txt per_timepoint_foldchanges.txt per_timepoint_qvalues.txt \
> consensus_peaks_plus_stats_plus_raw.txt

wait

rm \
    consensus_peak_stats.txt \
    consensus_peak_stats_avg.txt \
    consensus_peaks_plus_stats.txt \
    per_timepoint_foldchanges.txt \
    per_timepoint_qvalues.txt