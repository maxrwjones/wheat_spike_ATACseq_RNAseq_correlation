### Import gene list and name columns (standard GFF3 column names)
GFF3 <- read.table("outputs/HC_gene_coords.gff3", header=FALSE, sep = "\t")
colnames(GFF3) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

### Print dimensions and a few lines for checking purposes
dim(GFF3)
head(GFF3)
tail(GFF3)

### Extract TraesID from the "attributes" field then delete the latter
GFF3_mod <- GFF3
GFF3_mod$field1 <- str_split_i(GFF3_mod$attributes, ";", 1)
GFF3_mod$gene_ID <- str_split_i(GFF3_mod$field1, "=", 2)

GFF3_mod$attributes <- NULL
GFF3_mod$field1 <- NULL

head(GFF3_mod)



### Create a new DF to put results into
final_df <- as.data.frame(GFF3_mod$gene_ID)
colnames(final_df) <- "gene_ID"

### Add blank columns for adjacent genes' coordinates
final_df$chrom <- ""
final_df$extnd_start <- ""
final_df$extnd_end <- ""

### Print dimensions and a few lines for checking purposes
dim(final_df)
head(final_df)

#######################################################################################################

### For each gene except very first and last in annotation (these edge cases are set separately below)
for (i in 2:(nrow(final_df)-1)) {
    
    ### Store the TraesID of the target gene of interest in a variable
    target_gene <- final_df$gene_ID[i]
    
    ### Find the row indices of the target gene and its adjacent genes in the "GFF3" dataframe
    target_index <- which(GFF3_mod$gene_ID == target_gene)
    prev_index <- target_index-1
    next_index <- target_index+1 
    
    ### Extract the target gene's chromosome name from the "GFF3" DF and store it in the "targets" DF
    final_df$chrom[i] <- GFF3_mod$seqid[target_index]
    
########################################################################################################
    
    ### First need to check if target gene is the first one on its chromosome
    ### If it is, then "prev_index" will yield a different "seqid" to "target_index"
    ### This would lead to the wrong extended coord being applied.
    ### If the seqids are the same, then proceed to check if prev gene end is before current gene's start.
            ### If it is, take that gene's end coord as the extended start coord.
    ### If the seqid are not the same, then it must be the first gene on its chrom, so set the extended start coord to zero

    if (GFF3_mod$seqid[target_index] == GFF3_mod$seqid[prev_index]) {
        
        if (GFF3_mod$end[prev_index] < GFF3_mod$start[target_index]) {
            final_df$extnd_start[i] <- GFF3_mod$end[prev_index]
        } else {
            final_df$extnd_start[i] <- GFF3_mod$start[target_index]
        }
        
    } else {
        final_df$extnd_start[i] <- 0
    }

    
########################################################################################################
    
    ### First need to check if target gene is the last one on its chromosome
    ### If it is, then "next_index" will yield a different "seqid" to "target_index"
    ### This would lead to the wrong extended coord being applied.
    ### If the seqids are the same, then proceed to check if next gene start is after current gene's end.
            ### If it is, take that gene's start coord as the extended end coord.
    ### If the seqid are not the same, then it must be the last gene on its chrom, so set the extended start
    ### coord to 1 trillion (i.e. arbitrarily large)

    if (GFF3_mod$seqid[target_index] == GFF3_mod$seqid[next_index]) {
        
        if (GFF3_mod$start[next_index] > GFF3_mod$end[target_index]) {
            final_df$extnd_end[i] <- GFF3_mod$start[next_index]
        } else {
            final_df$extnd_end[i] <- GFF3_mod$end[target_index]
        }
        
    } else {
        final_df$extnd_end[i] <- 1000000000000
    }

}

########################################################################################################


### Input the chrom, extnd_start, and extnd_end coordinates of the very first gene in the annotation.
### Couldn't do this in above loop (without lots of mostly redundant checks) because, there isn't a valid "prev_index"
### for the very first gene.
target_gene <- final_df$gene_ID[1]
target_index <- which(GFF3_mod$gene_ID == target_gene)
next_index <- target_index+1
final_df$chrom[1] <- GFF3_mod$seqid[target_index]
final_df$extnd_start[1] <- 0
final_df$extnd_end[1] <- GFF3_mod$start[next_index]


### Input the chrom, extnd_start, and extnd_end coordinates of the very last gene in the annotation.
### Couldn't do this in above loop (without lots of mostly redundant checks) because, there isn't a valid "next_index"
### for the very last gene.
target_gene <- final_df$gene_ID[nrow(final_df)]
target_index <- which(GFF3_mod$gene_ID == target_gene)
prev_index <- target_index-1
final_df$chrom[nrow(final_df)] <- GFF3_mod$seqid[target_index]
final_df$extnd_start[nrow(final_df)] <- GFF3_mod$end[prev_index]
final_df$extnd_end[nrow(final_df)] <- 1000000000000


### Put the gene_ID column last in the DF to make it a valid BED file
final_df <- final_df %>% relocate(gene_ID, .after = last_col())

# View DF
dim(final_df)
head(final_df)

write.table(final_df, "outputs/HC_gene_coords_extended.bed", sep = "\t", quote = FALSE, col.names=FALSE, row.names=FALSE)