file_paths <- fs::dir_ls("outputs/stagewise_max_cov_files")

bed_files <- grep(".bed", file_paths, value=TRUE)
bed_files <- unname(bed_files)
bed_files

sample_names <- sub(".bed", "", bed_files)
sample_names <- sub("max_cov_per_peak_", "", sample_names)
sample_names <- sub(".", "_", sample_names, fixed=TRUE)
sample_names

file_contents <- list()
df_total <- data.frame()

for (i in seq_along(bed_files)) {
    file_contents[[i]] <- read.table(file = bed_files[i], sep='\t', header=FALSE)
    file_contents[[i]] <- file_contents[[i]][,-c(4,5,6)]
    colnames(file_contents[[i]]) <- c("chrom", "start", "end", "max_cov")
    names(file_contents[[i]])[names(file_contents[[i]]) == 'max_cov'] <- sample_names[i]
    if (nrow(df_total) == 0){
        df_total <- file_contents[[i]]
    } else {
        df_total <- merge(df_total, file_contents[[i]], by=c("chrom", "start", "end"))
    }
}

max_coverages <- df_total %>% arrange(chrom, start)
max_coverages

max_cov_long  <- max_coverages %>% pivot_longer(-c(chrom,start,end), names_to = "sample_id", values_to = "max_cov") %>%
                        separate(sample_id, into=c("stage","rep"), sep="_")

dim(max_cov_long)
head(max_cov_long)

### Calculate means based on developmental stage
max_cov_mean <- aggregate(max_cov~chrom+start+end+stage, max_cov_long, mean)

### Change name of count column to "mean_count" for accuracy
names(max_cov_mean)[names(max_cov_mean) == 'max_cov'] <- 'mean_max_cov'

### Pivot wider
max_cov_mean_wide <- max_cov_mean %>% pivot_wider(names_from="stage", values_from="mean_max_cov")

## Reorder columns
max_cov_mean_wide <- max_cov_mean_wide[,c(1,2,3,9,5,11,4,10,8,6,7)]

### View df
dim(max_cov_mean_wide)
head(max_cov_mean_wide)


### Remove rows where average in any mean_max_cov column is zero
mark_na <- function(x) ifelse(x == 0,NA,x)

max_cov_mean_trimmed <- max_cov_mean_wide

max_cov_mean_trimmed[,3:11] <- apply(max_cov_mean_trimmed[,3:11],2,mark_na)

max_cov_mean_trimmed <- na.omit(max_cov_mean_trimmed)

dim(max_cov_mean_trimmed)
head(max_cov_mean_trimmed)

### Resort peaks by chrom and start coordinate
max_cov_mean_trimmed_sorted <- max_cov_mean_trimmed %>% arrange(chrom, start)

dim(max_cov_mean_trimmed_sorted)
head(max_cov_mean_trimmed_sorted)

### Export DF as a BED file
write.table(max_cov_mean_trimmed_sorted, "outputs/consensus_peaks_with_mean_max_cov.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


#######################################################
### Also need a version where the reps are not averaged
#######################################################
### Reorder columns
max_cov_reps <- max_coverages[,c(1,2,3,14,15,6,7,18,19,4,5,16,17,12,13,8,9,10,11)]

### Set "start" and "end" as character variables to aid filtering step below
max_cov_reps$start <- as.character(max_cov_reps$start)
max_cov_reps$end <- as.character(max_cov_reps$end)

### View DF
dim(max_cov_reps)
head(max_cov_reps)

### Remove rows where all the peak values are zero
### i.e. keep rows where number of zeroes is less than number of cols - 3 (to account for "chrom", "start", and "end" cols)
max_cov_reps_no_zeros <- max_cov_reps[rowSums(max_cov_reps==0, na.rm=TRUE) < ncol(max_cov_reps) - 3, ]

### Vief DF
dim(max_cov_reps_no_zeros)
head(max_cov_reps_no_zeros)

### Export DF as a BED file
write.table(max_cov_reps_no_zeros, "outputs/consensus_peaks_with_max_cov.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)