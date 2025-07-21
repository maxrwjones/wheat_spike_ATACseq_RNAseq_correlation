peaks <- read.table("consensus_peaks_with_parent_genes.bed", sep='\t', header=FALSE)

dim(peaks)
head(peaks)

peaks_red <- peaks[,-1:-3]
colnames(peaks_red) <- c("parent_gene", "chrom", "start", "end", "SAM_1", "SAM_2", "EL_1", "EL_2", "SR_1", "SR_2", "DR_1", "DR_2", "SMI_1", "SMI_2", "GPD_1", "GPD_2", "FMI_1", "FMI_2", "FOP_1", "FOP_2")
peaks_sorted <- peaks_red %>% arrange(chrom, start)

dim(peaks_sorted)
head(peaks_sorted)

### Create a copy of "peaks_sorted"
peaks_final <- peaks_sorted

### Create a unique identifier for each peak and store as the row names
rownames(peaks_final) <- paste0(peaks_final$parent_gene, "_P", rownames(peaks_final))

### Remove "parent_gene", "chrom", "start", and "end" columns
peaks_final <- peaks_final[,5:20]

### Round to nearest integer
peaks_final <- peaks_final %>% 
    mutate_if(is.numeric, ~round(.,0))

### Export this DF
write.csv(peaks_final, "peak_profiles_ImpulseDE2_format.csv")

### View df
dim(peaks_final)
head(peaks_final)


### Create an annotation DF to run ImpulseDE2 in case-only mode on these peak profiles
### Omit SAM in the annotation

Sample <- colnames(peaks_final[3:16])
Condition <- rep("case", each=14)
Batch <- rep("BNULL", each=14)
Time <- rep(c(1,2,3,4,5,6,7), each=2)

annotation <- data.frame(Sample, Condition, Batch, Time)

rownames(annotation) <- annotation$Sample

annotation

write.csv(annotation, "peak_profiles_ImpulseDE2_annotation.csv")


#####################################################################################
### Run ImpulseDE2 on peak data #####################################################
#####################################################################################


library(ImpulseDE2)

annotation <- read.csv("peak_profiles_ImpulseDE2_annotation.csv", header=TRUE, row.names=1)
counts <- read.csv("peak_profiles_ImpulseDE2_format.csv", header=TRUE, row.names=1)
counts  <- as.matrix(counts)

objectImpulseDE2 <- runImpulseDE2(
  matCountData    = counts,
  dfAnnotation    = annotation,
  boolCaseCtrl    = FALSE, # Is this a case-control analysis?
  vecConfounders  = NULL, # Any factors to account for on batch correction?
  scaNProc        = 4 ) # Number of processes for parallelisation

# Save vecDEgenes (genes identified as differentially expressed by ImpusleDE at threshold scaQThres)
write.csv(objectImpulseDE2$VecDEgenes, "peaks_vecDEGenes.csv")

# Save dfDEAnalysis (dataframe of samples x reported characteristics)
write.csv(objectImpulseDE2$dfDEAnalysis, "peaks_DEAnalysis.csv")

# Save dfImpulseDE2Results (dataframe of results per gene such as p, padj, loglik)
write.csv(objectImpulseDE2$dfImpulseDE2Results, "peaks_ImpulseDE2_results.csv")