args <- commandArgs()
array_ID <- args[6]

suppressMessages(library(tidyverse))

sig_gene_tpms_norm <- read.table("sig_gene_tpms_norm.tsv", header=TRUE)
sig_peak_covs_norm_inverted <- read.table("sig_peak_covs_norm_inverted.tsv", header=TRUE)

### Randomly assign a sig DE gene to be the parent of each peak
### Use the array task as the seed for random number draws to ensure each array task produces different results
set.seed(array_ID)
sig_peak_covs_shuff <- transform(sig_peak_covs_norm_inverted, parent = sample(sig_peak_covs_norm_inverted$parent, replace=FALSE))

SS_profile_list <- list()
k <- 0

for (i in 1:nrow(sig_gene_tpms_norm)){
    
    ### Put the numeric columns for the current gene (i.e. the tpm trajectory) in a vector
    gene_profile <- as.numeric(sig_gene_tpms_norm[i,] %>% select_if(is.numeric))
    
    ### Store the gene's ID
    gene_id <- sig_gene_tpms_norm[i,1]
    
    ### Find all significant peaks whose parent gene matches the current gene
    assoc_peaks <- sig_peak_covs_shuff %>% filter(parent == gene_id)
    
    ### If there are no significant peaks for current gene, move to next loop iteration (i.e. next gene)
    ### Otherwise, loop through each peak in turn
    if (nrow(assoc_peaks) == 0) {
        next
    } else {
        for (j in 1:nrow(assoc_peaks)){

            k <- k + 1
            
            ### Put the numeric columns for the current peak (i.e. the coverage trajectory) in a vector
            peak_profile <- as.numeric(assoc_peaks[j,] %>% select_if(is.numeric))
            
            ### Store the peak's ID
            peak_id <- assoc_peaks$peak[j]
            
            ### At each timepoint, calculate the mean of the gene and peak profiles
            mean_profile <- (gene_profile + peak_profile) / 2
            
            ### Calculate the sum of squares from the mean for each timepoint
            SS_profile <- (gene_profile - mean_profile)^2 + (peak_profile - mean_profile)^2
            
            ### Add current peak's data to a growing list of vectors
            SS_profile_list[[k]] <- c(gene_id, peak_id, SS_profile)
            
        }
    }
}

### Bind the individual rows from SS_profile_list into a dataframe
all_SS_profiles <- as.data.frame(do.call(rbind, SS_profile_list))

### Rename columns
colnames(all_SS_profiles) <- c("gene", "peak", "T", "EDR", "LDR", "GP", "LP", "FP", "TS")

### Change the data type of the SS columns to numeric
all_SS_profiles[,3:9] <- sapply(all_SS_profiles[,3:9], as.numeric)

### Calculate the sum of SS for each gene-peak combo
all_SS_profiles$sum <- apply(all_SS_profiles[,3:9],1,sum)

### Calculate the mean SS for each gene-peak combo
all_SS_profiles$mean <- apply(all_SS_profiles[,3:9],1,mean)

### Calculate the median SS for each gene-peak combo
all_SS_profiles$median <- apply(all_SS_profiles[,3:9],1,median)


### Export putative silencers
filename <- paste0("silencer_simulations/simulated_silencer_SS_profiles_", array_ID, ".tsv")
write.table(all_SS_profiles, filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)