setwd("C:/Users/sandr/OneDrive/Desktop/RNAseq")
library(tximport)
library(readr)
library(dplyr)

# Read in mapping
tx2gene <- read.table("tx2gene.txt", 
                      col.names = c("transcript_id", "gene_id"),
                      stringsAsFactors = FALSE)

# read in counts
counts <- read.table("counts_matrix.sf", header = TRUE, sep = "\t", row.names = 1)

# quickcheck of matrix
head(counts)

# Sort the sample names
colnames(counts) <- sub("_.*", "", colnames(counts))


# add transcript ids as a col pre-merge
counts$transcript_id <- rownames(counts)

# Merge with tx2gene mapping
counts_with_genes <- merge(counts, tx2gene, by = "transcript_id", all.x = TRUE)

### 
#> matched <- sum(!is.na(counts_with_genes$gene_id))
#> matched
#[1] 29931

#> total <- nrow(counts_with_genes)
#> total
#[1] 34479
###

# Remove rows without gene mapping
counts_with_genes_filtered <- counts_with_genes[!is.na(counts_with_genes$gene_id), ]

colnames(counts_with_genes_filtered)

#rm transcript cols
counts_with_genes_filtered$transcript_id <- NULL
counts_with_genes_filtered$transcript <- NULL

gene_counts <- counts_with_genes_filtered %>%
  group_by(gene_id) %>%
  summarise(across(everything(), sum)) #dplyr

# make it a df and nullify gene_id col
gene_counts <- as.data.frame(gene_counts)
rownames(gene_counts) <- gene_counts$gene_id
gene_counts$gene_id <- NULL

dim(gene_counts)
head(gene_counts[, 1:5])

### Save matrix ###
#write.table(gene_counts, "gene_level_counts.txt", sep = "\t")
write.table(gene_counts, "gene_level_counts.txt", sep = "\t", quote = F)
