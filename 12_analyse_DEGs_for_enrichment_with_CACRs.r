suppressMessages(library(tidyverse))

################################################################################
### Load in the gene and peak Impulse results and do some minor formatting #####
################################################################################

gene_results <- read.csv("RNAseq_ImpulseDE2_results.csv.csv", header=TRUE)
peak_results <- read.csv("peaks_ImpulseDE2_results.csv", header=TRUE)

gene_results$Gene <- NULL
peak_results$Gene <- NULL

colnames(gene_results)[colnames(gene_results) == "X"] <- "gene"
colnames(peak_results)[colnames(peak_results) == "X"] <- "peak"

dim(gene_results)

### Remove genes on ChrUn (consensus peaks were not calculated for these genes as "adjacent genes" concept not really applicable)
gene_results <- gene_results %>% filter(!grepl("TraesCSU02G", gene))

dim(gene_results)
head(gene_results)

dim(peak_results)
head(peak_results)

### Split peak names into an ID and a parent gene

peak_parents <- separate_wider_delim(peak_results, cols = peak, delim = "_", names = c("parent", "peak"))

peak_parents


################################################################################
### Start analysing genes for DE across developmental trajectory ###############
################################################################################

### Filter for the subset of significantly DE genes
sig_genes <- gene_results %>% filter(padj < 0.05)

dim(sig_genes)
head(sig_genes)

### Filter for the subset of non-DE genes
nonsig_genes <- gene_results %>% filter(padj >= 0.05)

dim(nonsig_genes)
head(nonsig_genes)

### Filter for the subset of non-DE genes (including unexpressed genes)
nonsig_and_NA_genes <- gene_results %>% filter(padj > 0.05 | is.na(padj))

dim(nonsig_and_NA_genes)
head(nonsig_and_NA_genes)

### Filter for the subset of unexpressed genes
NA_genes <- gene_results %>% filter(is.na(padj))

dim(NA_genes)
head(NA_genes)


### Filter with a more stringent significance threshold
sig_genes_0001 <- gene_results %>% filter(padj < 0.001)
dim(sig_genes_0001)

### Export DF of padj < 0.05 significant genes
write.table(sig_genes, "sig_0.05_genes.tsv", col.names=TRUE, row.names=FALSE, quote = FALSE, sep = "\t")

### Export DF of padj < 0.001 significant genes
write.table(sig_genes_0001, "sig_0.001_genes.tsv", col.names=TRUE, row.names=FALSE, quote = FALSE, sep = "\t")


################################################################################
### Start analysing peaks for significations DE parent genes ###################
################################################################################

#### Pull out peaks with significantly DE parent genes
peaks_with_sig_parents <- subset(peak_parents, parent %in% sig_genes$gene)

cat("Total number of any peaks with sig parents:")
nrow(peaks_with_sig_parents)

cat("Total number of sig peaks with sig parents:")
nrow(peaks_with_sig_parents %>% filter(padj < 0.05))

cat("Average number of any peaks per sig parent")
nrow(peaks_with_sig_parents) / nrow(sig_genes)

cat("Average number of sig peaks per sig parent:")
nrow(peaks_with_sig_parents %>% filter(padj < 0.05)) / nrow(sig_genes)

dim(peaks_with_sig_parents)
head(peaks_with_sig_parents)

################################################################################

### Pull out peaks with non-DE parents
peaks_with_nonsig_parents <- subset(peak_parents, parent %in% nonsig_genes$gene)

cat("Total number of any peaks with non-sig parents:")
nrow(peaks_with_nonsig_parents)

cat("Total number of sig peaks with non-sig parents:")
nrow(peaks_with_nonsig_parents %>% filter(padj < 0.05))

cat("Average number of any peaks per non-sig parent")
nrow(peaks_with_nonsig_parents) / nrow(nonsig_genes)

cat("Average number of sig peaks per non-sig parent:")
nrow(peaks_with_nonsig_parents %>% filter(padj < 0.05)) / nrow(nonsig_genes)

dim(peaks_with_nonsig_parents)
head(peaks_with_nonsig_parents)

################################################################################

### Pull out peaks with non-DE or unexpressed parents
peaks_with_nonsig_or_NA_parents <- subset(peak_parents, parent %in% nonsig_and_NA_genes$gene)

cat("Total number of any peaks with non-sig or unexpressed parents:")
nrow(peaks_with_nonsig_or_NA_parents)

cat("Total number of sig peaks with non-sig or unexpressed parents:")
nrow(peaks_with_nonsig_or_NA_parents %>% filter(padj < 0.05))

cat("Average number of any peaks per non-sig or unexpressed parent")
nrow(peaks_with_nonsig_or_NA_parents) / nrow(nonsig_and_NA_genes)

cat("Average number of sig peaks per non-sig or unexpressed parent:")
nrow(peaks_with_nonsig_or_NA_parents %>% filter(padj < 0.05)) / nrow(nonsig_and_NA_genes)

dim(peaks_with_nonsig_or_NA_parents)
head(peaks_with_nonsig_or_NA_parents)

################################################################################
### Pull out peaks with unexpressed parents
peaks_with_NA_parents <- subset(peak_parents, parent %in% NA_genes$gene)

cat("Total number of any peaks with unexpressed parents:")
nrow(peaks_with_NA_parents)

cat("Total number of sig peaks with unexpressed parents:")
nrow(peaks_with_NA_parents %>% filter(padj < 0.05))

cat("Average number of any peaks per unexpressed parent")
nrow(peaks_with_NA_parents) / nrow(NA_genes)

cat("Average number of sig peaks per unexpressed parent:")
nrow(peaks_with_NA_parents %>% filter(padj < 0.05)) / nrow(NA_genes)

dim(peaks_with_NA_parents)
head(peaks_with_NA_parents)


################################################################################
### Shuffle peaks' parent genes randomly to generate a null distribution #######
################################################################################

### Create a set of empty vectors which will be populated with summary statistics from the peak sets with random parents
AP_SG <- c()    # Vector for "average number of any peaks per sig parent"
SP_SG <- c()    # Vector for "average number of sig peaks per sig parent"

AP_NSG <- c()    # Vector for "average number of any peaks per non-sig parent"
SP_NSG <- c()    # Vector for "average number of sig peaks per non-sig parent"

AP_NSUG <- c()    # Vector for "average number of any peaks per non-sig or unexpressed parent"
SP_NSUG <- c()    # Vector for "average number of sig peaks per non-sig or unexpressed parent"

AP_UG <- c()    # Vector for "average number of any peaks per unexpressed parent"
SP_UG <- c()    # Vector for "average number of sig peaks per unexpressed parent"

first_parents <- c() # Vector to store the newly assigned 'parent' gene of the first peak in each loop for checking purposes


for (i in 1:500){
    
    ### Shuffle which parent genes are associated with which peak (sampling genes with replacement)
    shuffled_parents <- transform(peak_parents, parent = sample(gene_results$gene, replace=TRUE, size=nrow(peak_parents) ) )
    
    ## Store the new 'parent' of the first peak in a vector after each shuffle for checking purposes
    first_parents <- c(first_parents, head(shuffled_parents$parent)[1])
    
    ###############################################################################################                              
                                  
    ### Filter "shuffled_parents" for peaks with significantly DE 'parent' genes
    shuff_peaks_with_sig_parents <- subset(shuffled_parents, parent %in% sig_genes$gene)
    
    ### Filter "shuffled_parents" for peaks with non-DE 'parent' genes
    shuff_peaks_with_nonsig_parents <- subset(shuffled_parents, parent %in% nonsig_genes$gene)
    
    ### Filter "shuffled_parents" for peaks with non-DE or unexpressed 'parent' genes                              
    shuff_peaks_with_nonsig_or_NA_parents <- subset(shuffled_parents, parent %in% nonsig_and_NA_genes$gene)
    
    ### Filter "shuffled_parents" for peaks with unexpressed 'parent' genes                                   
    shuff_peaks_with_NA_parents <- subset(shuffled_parents, parent %in% NA_genes$gene)
                                  
    ###############################################################################################
    
    ### Calculate average number of any peaks per sig parent
    AP_SG <- c(AP_SG, nrow(shuff_peaks_with_sig_parents) / nrow(sig_genes))
    ### Calculate average number of sig peaks per sig parent
    SP_SG <- c(SP_SG, nrow(shuff_peaks_with_sig_parents %>% filter(padj < 0.05)) / nrow(sig_genes))

    
    ### Calculate average number of any peaks per non-sig parent
    AP_NSG <- c(AP_NSG, nrow(shuff_peaks_with_nonsig_parents) / nrow(nonsig_genes))
    ### Calculate average number of sig peaks per non-sig parent
    SP_NSG <- c(SP_NSG, nrow(shuff_peaks_with_nonsig_parents %>% filter(padj < 0.05)) / nrow(nonsig_genes))
                          
                                  
    ### Calculate average number of any peaks per non-sig or unexpressed parent
    AP_NSUG <- c(AP_NSUG, nrow(shuff_peaks_with_nonsig_or_NA_parents) / nrow(nonsig_and_NA_genes))
    ### Calculate average number of sig peaks per non-sig or unexpressed parent
    SP_NSUG <- c(SP_NSUG, nrow(shuff_peaks_with_nonsig_or_NA_parents %>% filter(padj < 0.05)) / nrow(nonsig_and_NA_genes))
                                  
    
    ### Calculate average number of any peaks per unexpressed parent
    AP_UG  <- c(AP_UG, nrow(shuff_peaks_with_NA_parents) / nrow(NA_genes))
    ### Calculate average number of sig peaks per unexpressed parent
    SP_UG  <- c(SP_UG, nrow(shuff_peaks_with_NA_parents %>% filter(padj < 0.05)) / nrow(NA_genes))

}

####################################################################################################################

cat("Check which genes got assigned to be parent of first peak in the first 20 loops - should be randomly selected from all HC, non-ChrUn genes")
head(first_parents, n=20)

cat("Check the first 20 values for average number of any peaks per sig parent")
head(AP_SG, n=20)
cat("Check the first 20 values for average number of sig peaks per sig parent")
head(SP_SG, n=20)

cat("The average number of any peaks per sig parent across 500 shuffled datasets is: ")
cat(mean(AP_SG))
cat("\nAnd the standard deviation of this is: ")
cat(sd(AP_SG))

cat("\n\n")

cat("The average number of sig peaks per sig parent across 500 shuffled datasets is: ")
cat(mean(SP_SG))
cat("\nAnd the standard deviation of this is: ")
cat(sd(SP_SG))


################################################################################
### Plot histogram of numbers of any peaks per sig DE gene #####################
################################################################################

### Define some variables for lines I want to overlay
real_data_mean <- nrow(peaks_with_sig_parents) / nrow(sig_genes)
mean_peaks <- mean(AP_SG)
mean_plus_sd <- mean(AP_SG) + sd(AP_SG)
mean_minus_sd <- mean(AP_SG) - sd(AP_SG)
mean_plus_5sd <- mean(AP_SG) + 5*sd(AP_SG)
mean_minus_5sd <- mean(AP_SG) - 5*sd(AP_SG)

### Plot the data
ggplot() + aes(AP_SG) +
    geom_histogram(binwidth=0.001, colour="#000000", fill="#0099F8") +
    geom_segment(aes(x=real_data_mean, y=0, xend=real_data_mean, yend=50), color="red", linetype="solid", size=1) +
    geom_segment(aes(x=mean_peaks, y=0, xend=mean_peaks, yend=50), color="deeppink", linetype="solid", size=1) +
    geom_segment(aes(x=mean_plus_sd, y=0, xend=mean_plus_sd, yend=50), color="deeppink", linetype="dashed", size=1) +
    geom_segment(aes(x=mean_minus_sd, y=0, xend=mean_minus_sd, yend=50), color="deeppink", linetype="dashed", size=1) +
    geom_segment(aes(x=mean_plus_5sd, y=0, xend=mean_plus_5sd, yend=50), color="purple", linetype="dashed", size=1) +
    geom_segment(aes(x=mean_minus_5sd, y=0, xend=mean_minus_5sd, yend=50), color="purple", linetype="dashed", size=1) +
    labs(title = "Numbers of any CACRs per significantly DE gene in 500 shuffled datasets",
         x = "Number of CACRs",
         y = "Count") +
    theme_classic() +
    theme(axis.title = element_text(size = rel(2.5)),
          axis.text = element_text(size = rel(2)),
          axis.ticks = element_line(linewidth = rel(2)),
          plot.title = element_text(size=rel(3))
         ) +
    scale_y_continuous(expand = c(0,0))

    ggsave("CACRs_per_DEG.svg", height=12, width=50, unit="cm", limitsize = FALSE)

################################################################################
### Plot histogram of numbers of sig peaks per sig DE gene #####################
################################################################################

### Define some variables for lines I want to overlay
real_data_mean <- nrow(peaks_with_sig_parents %>% filter(padj < 0.05)) / nrow(sig_genes)
mean_peaks <- mean(SP_SG)
mean_plus_sd <- mean(SP_SG) + sd(SP_SG)
mean_minus_sd <- mean(SP_SG) - sd(SP_SG)
mean_plus_5sd <- mean(SP_SG) + 5*sd(SP_SG)
mean_minus_5sd <- mean(SP_SG) - 5*sd(SP_SG)

### Plot the data
ggplot() + aes(SP_SG) +
    geom_histogram(binwidth=0.001, colour="#000000", fill="#00bfc4") +
    geom_segment(aes(x=real_data_mean, y=0, xend=real_data_mean, yend=50), color="#f8766d", linetype="solid", size=1) +
    geom_segment(aes(x=mean_peaks, y=0, xend=mean_peaks, yend=50), color="deeppink", linetype="solid", size=1) +
    geom_segment(aes(x=mean_plus_sd, y=0, xend=mean_plus_sd, yend=50), color="deeppink", linetype="dashed", size=1) +
    geom_segment(aes(x=mean_minus_sd, y=0, xend=mean_minus_sd, yend=50), color="deeppink", linetype="dashed", size=1) +
    geom_segment(aes(x=mean_plus_5sd, y=0, xend=mean_plus_5sd, yend=50), color="purple", linetype="dashed", size=1) +
    geom_segment(aes(x=mean_minus_5sd, y=0, xend=mean_minus_5sd, yend=50), color="purple", linetype="dashed", size=1) +
    labs(title = "Numbers of significantly changing peaks per significantly DE gene in 500 shuffled datasets",
         x = "Average number of dCACRs per DEG",
         y = "Number of simulations") +
    theme_classic() +
    theme(axis.title = element_text(size = rel(2.5)),
          axis.text = element_text(size = rel(2), face="bold"),
          axis.ticks = element_line(linewidth = rel(2)),
          plot.title = element_text(size=rel(3))
         ) +
    scale_y_continuous(expand = c(0,0))

ggsave("dCACRs_per_DEG.svg", height=12, width=50, unit="cm", limitsize = FALSE)