setwd("WORKING_DIR_HERE")


### A vector is made that lists all the file paths from the cluster folder where I stored the abundance files output by Kallisto.
### Next, the vector column names are removed using "unname" because they get in the way. The file paths are then edited to remove most of the path and the .tsv suffix, leaving just the sample names.
### An empty dataframe "df_total" is generated.
#### Then each .tsv file is read into an element of the list "file_contents". All columns except "target_id" and "est counts" are removed. The "est_counts" columns are then renamed using the "sample_names" vector defined in the previous cell.
### The if statement causes the first "target_id" column to be retained while all others are removed. (Using "merge".)

file_paths <- fs::dir_ls("inputs/RNAseq_abundance_files")
file_paths2 <- unname(file_paths)
file_paths2

sample_names <- sub(".tsv", "", file_paths2)
sample_names <- sub("inputs/abundance_files/", "", sample_names)
sample_names

file_contents <- list()
df_total <- data.frame()

for (i in seq_along(file_paths)) {
    file_contents[[i]] <- read.table(file = file_paths[i], sep='\t', header=TRUE)
    file_contents[[i]] <- select(file_contents[[i]], 'target_id', 'est_counts')
    names(file_contents[[i]])[names(file_contents[[i]]) == 'est_counts'] <- sample_names[i]
    if (nrow(df_total) == 0){
        df_total <- file_contents[[i]]
    } else {
        df_total <- merge(df_total, file_contents[[i]], by='target_id')
    }
}

df_total


### A copy of df_total is made and manipulated such that all transcripts of a given gene are collapsed, i.e. counts summed.
counts <- df_total
counts <- counts %>% separate(target_id, c('gene_id', 'b'))
counts$b <- NULL
counts <- aggregate(counts[,2:25], by=counts['gene_id'], sum)

colnames(counts) <- c("gene_id", "SAM_1", "SAM_2", "SAM_3", "T_1", "T_2", "T_3", "EDR_1", "EDR_2", "EDR_3", "LDR_1", "LDR_2", "LDR_3", "GP_1", "GP_2", "GP_3", "LP_1", "LP_2", "LP_3", "FP_1", "FP_2", "FP_3", "TS_1", "TS_2", "TS_3")

dim(counts)
head(counts)


### Round counts to nearest integer (ImpulseDE2 only accepts integers)
roundedcounts <- counts %>% 
    mutate_if(is.numeric, ~round(.,0))

### Create a version of the rounded counts dataset with just the HC genes
roundedcounts_HC <- roundedcounts %>% filter(!grepl("LC", gene_id))

rownames(roundedcounts_HC) <- roundedcounts_HC$gene_id
roundedcounts_HC$gene_id <- NULL

### Export df of collapsed, rounded gene counts for HC genes only
write.csv(roundedcounts_HC, "RNAseq_counts_ImpulseDE2_format.csv")

### View df
dim(roundedcounts_HC)
head(roundedcounts_HC)



### Create an annotation DF to run ImpulseDE2 in case-only mode on this RNA-seq count data
Sample <- colnames(roundedcounts_HC[5:25])
Condition <- rep("case", each=21)
Batch <- rep("BNULL", each=21)
Time <- rep(c(1,2,3,4,5,6,7),each=3)

annotation <- data.frame(Sample, Condition, Batch, Time)

rownames(annotation) <- annotation$Sample

annotation

write.csv(annotation, "RNAseq_counts_ImpulseDE2_annotation.csv")



#####################################################################################
### Run ImpulseDE2 on transcript data ###############################################
#####################################################################################

library(ImpulseDE2)

annotation <- read.csv("RNAseq_counts_ImpulseDE2_annotation.csv", header=TRUE, row.names=1)
counts <- read.csv("RNAseq_counts_ImpulseDE2_format.csv", header=TRUE, row.names=1)
counts  <- as.matrix(counts)

objectImpulseDE2 <- runImpulseDE2(
  matCountData    = counts,
  dfAnnotation    = annotation,
  boolCaseCtrl    = FALSE, # Is this a case-control analysis?
  vecConfounders  = NULL, # Any factors to account for on batch correction?
  scaNProc        = 4 ) # Number of processes for parallelisation

# Save vecDEgenes (genes identified as differentially expressed by ImpusleDE at threshold scaQThres)
write.csv(objectImpulseDE2$VecDEgenes, "RNAseq_vecDEGenes.csv")

# Save dfDEAnalysis (dataframe of samples x reported characteristics)
write.csv(objectImpulseDE2$dfDEAnalysis, "RNAseq_DEAnalysis.csv")

# Save dfImpulseDE2Results (dataframe of results per gene such as p, padj, loglik)
write.csv(objectImpulseDE2$dfImpulseDE2Results, "RNAseq_ImpulseDE2_results.csv")




#################################################################################
# Repeat wrangling of RNA-seq data but for TPMs, not counts. Needed for step 14 #
#################################################################################

### A vector is made that lists all the file paths from the cluster folder where I stored the abundance files output by Kallisto.
### Next, the vector column names are removed using "unname" because they get in the way. The file paths are then edited to remove most of the path and the .tsv suffix, leaving just the sample names.
### An empty dataframe "df_total" is generated.
### Then each .tsv file is read into an element of the list "file_contents". All columns except "target_id" and "tpm" are removed. The "est_counts" columns are then renamed using the "sample_names" vector defined in the previous cell.
### The if statement causes the first "target_id" column to be retained while all others are removed. (Using "merge".)

file_paths <- fs::dir_ls("abundance_files")
file_paths2 <- unname(file_paths)
file_paths2

sample_names <- sub(".tsv", "", file_paths2)
sample_names <- sub("abundance_files/", "", sample_names)
sample_names

file_contents <- list()
df_total <- data.frame()

for (i in seq_along(file_paths)) {
    file_contents[[i]] <- read.table(file = file_paths[i], sep='\t', header=TRUE)
    file_contents[[i]] <- select(file_contents[[i]], 'target_id', 'tpm')
    names(file_contents[[i]])[names(file_contents[[i]]) == 'tpm'] <- sample_names[i]
    if (nrow(df_total) == 0){
        df_total <- file_contents[[i]]
    } else {
        df_total <- merge(df_total, file_contents[[i]], by='target_id')
    }
}

df_total


### A copy of df_total is made and manipulated such that all transcripts of a given gene are collapsed, i.e. counts summed.
### The DF is saved as a .csv
tpms <- df_total
tpms <- tpms %>% separate(target_id, c('gene_id', 'b'))
tpms$b <- NULL
tpms <- aggregate(tpms[,2:25], by=tpms['gene_id'], sum)

colnames(tpms) <- c("gene_id", "SAM_1", "SAM_2", "SAM_3", "T_1", "T_2", "T_3", "EDR_1", "EDR_2", "EDR_3", "LDR_1", "LDR_2", "LDR_3", "GP_1", "GP_2", "GP_3", "LP_1", "LP_2", "LP_3", "FP_1", "FP_2", "FP_3", "TS_1", "TS_2", "TS_3")

dim(tpms)
head(tpms)
write.table(tpms, "collapsed_genes_tpms.tsv", sep="\t", quote=FALSE, row.names=FALSE)


### Pivot longer then split sample name into stage and rep
tpms_long  <- tpms %>% pivot_longer(-gene_id, names_to = "sample_id", values_to = "tpm") %>%
                        separate(sample_id, into=c("stage","rep"), sep="_")

dim(tpms_long)
head(tpms_long)

### Calculate means based on developmental stage
mean_tpms <- aggregate(tpm~gene_id+stage, tpms_long, mean)

### Change name of tpm column to "mean_tpm" for accuracy
names(mean_tpms)[names(mean_tpms) == 'tpm'] <- 'mean_tpm'

### View df
dim(mean_tpms)
head(mean_tpms)

### Pivot wider and reorder columns
mean_tpms_wide <- mean_tpms %>% pivot_wider(names_from="stage", values_from="mean_tpm")

mean_tpms_wide <- mean_tpms_wide[,c(1,7,8,2,5,4,6,3,9)]

### Export df
write.table(mean_tpms_wide, "mean_tpms.tsv", sep="\t", quote=FALSE, row.names=FALSE)

### View df
dim(mean_tpms_wide)
head(mean_tpms_wide)

mean_tpms_wide_HC <- mean_tpms_wide %>% filter(!grepl("LC", gene_id))

### Export df
write.table(mean_tpms_wide_HC, "mean_tpms_HC_only.tsv", sep="\t", quote=FALSE, row.names=FALSE)

### View df
dim(mean_tpms_wide_HC)
head(mean_tpms_wide_HC)