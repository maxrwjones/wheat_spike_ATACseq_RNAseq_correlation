suppressMessages(library(tidyverse))

### Read in DF of genes with padj < 0.05
sig_genes <- read.table("sig_0.05_genes.tsv", header=TRUE)

### Read in DF of peaks ImpulseDE2 results
peaks <- read.csv("peaks_ImpulseDE2_results.csv")

### Remove unneeded column
peaks$X <- NULL

### Filter for peaks with padj < 0.05
sig_peaks <- peaks %>% filter(padj < 0.05)

### Split peak names into an ID and a parent gene
sig_peaks <- separate_wider_delim(sig_peaks, cols = Gene, delim = "_", names = c("parent", "peak"))

### Read in DF of count data per gene/timepoint combo (already calculated means across reps)
gene_tpms <- read.table("mean_tpms_HC_only.tsv", sep="\t", header=TRUE)

### Rename column 'gene_id' to 'gene'
colnames(gene_tpms)[colnames(gene_tpms) == "gene_id"] <- "gene"

### Remove 'SAM' column as we won't use it for correlations
gene_tpms$SAM <- NULL

### Read in DF of max_cov data per peak/timepoint/rep combo
peak_covs <- read.csv("peak_profiles_ImpulseDE2_format.csv")

### Remove SAM columns as we won't use them for correlations
peak_covs$SAM_1 <- NULL
peak_covs$SAM_2 <- NULL


### Calculate means for each timepoint for each peak
peak_covs_long  <- peak_covs %>% pivot_longer(-X, names_to = "sample_id", values_to = "cov") %>%
                        separate(sample_id, into=c("stage","rep"), sep="_")
peak_covs_long_means <- aggregate(cov~X+stage, peak_covs_long, mean)
peak_mean_covs <- peak_covs_long_means %>% pivot_wider(names_from="stage", values_from="cov")

### Rearrange column order
peak_mean_covs <- peak_mean_covs[,c(1,3,8,2,7,6,4,5)]

### Rename to match Waddington stage names used for gene profiles
colnames(peak_mean_covs) <- c("X", "T", "EDR", "LDR", "GP", "LP", "FP", "TS")

### Split peak names into an ID and a parent gene
peak_mean_covs <- separate_wider_delim(peak_mean_covs, cols = X, delim = "_", names = c("parent", "peak"))

### Filter 'gene_counts' DF based on 'sig_genes' to retain only significantly DE genes
sig_gene_tpms <- gene_tpms %>% filter(gene %in% sig_genes$gene)

### Filter 'peak_covs' DF based on 'sig_peaks' to retain only significantly DE genes
sig_peak_covs <- peak_mean_covs %>% filter(peak %in% sig_peaks$peak)

### Define a custom function to normalise by mean centering and dividing by standard deviation
normalise_function <- function(row) {
  row_array <- as.numeric(row)
  mean_centered <- row_array - mean(row_array, na.rm = TRUE)
  std_dev <- sd(row_array, na.rm = TRUE)
  if (std_dev != 0) {  # To avoid division by zero
    normalised <- mean_centered / std_dev
  } else {
    normalised <- mean_centered
  }
  return(normalised)
}

### Create copies of sig_gene_tpms and sig_peak_covs to store normalised results in
sig_gene_tpms_norm <- sig_gene_tpms
sig_peak_covs_norm <- sig_peak_covs

### Apply the custom function row-wise to columns 2 to 8 of sig_gene_tpms to normalise the gene tpm data
sig_gene_tpms_norm[, 2:8] <- t(apply(sig_gene_tpms[, 2:8], 1, normalise_function))

### Apply the custom function row-wise to columns 3 to 9 of sig_peak_covs to normalise the peak coverage data
sig_peak_covs_norm[, 3:9] <- t(apply(sig_peak_covs[, 3:9], 1, normalise_function))

### Create a copy of sig_peak_covs_norm which can then be inverted in the next step
sig_peak_covs_norm_inverted <- sig_peak_covs_norm

### Multiply all peak coverage values by -1 to invert the trajectory pattern (so that we can find silencers)
sig_peak_covs_norm_inverted[, c(-1,-2)] <- sig_peak_covs_norm_inverted[, c(-1,-2)] * -1


### Export normalised datasets
write.table(sig_gene_tpms_norm, "sig_gene_tpms_norm.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(sig_peak_covs_norm, "sig_peak_covs_norm.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(sig_peak_covs_norm_inverted, "sig_peak_covs_norm_inverted.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)