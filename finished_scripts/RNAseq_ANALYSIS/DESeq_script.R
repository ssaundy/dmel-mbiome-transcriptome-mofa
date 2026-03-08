##### THIS IS MY RNAseq ANALYSIS FOR OUR 117 DMEL GUTS #####
setwd("C:/Users/sandr/OneDrive - University of Glasgow/PROJECT/DATA/TRANSCRIPTOME")
source("OMICS_rnaseq.R")
### Please note that the counts matrix was created on the HPC using the below command ###
### salmon quantmerge --quants DS*_L001 --column NumReads --output counts_matrix.sf ###

# Bring in the gene level matrix
em = read.table("gene_level_counts.txt",
                    header=TRUE,
                    sep="\t",
                    row.names = 1)
em

# Sort the sample names
colnames(em) <- sub("_.*", "", colnames(em))

# bring in sample metadata (rows = samples, cols = population, mito, nuclear, replicate #)
samp = read.table("C:/Users/sandr/OneDrive/Desktop/mbio/samples.csv",
                  header=TRUE,
                  sep = ",",
                  row.names = 1)
samp

##### Sort the gene symbol issue (lack of) #####
# Load libraries
library(org.Dm.eg.db)
library(dplyr)
library(AnnotationDbi)

# Get the gene IDs from rownames
gene_ids <- rownames(em)

# remove gene-Dmel_ from rownames - retaining the flybase "CG" ID
clean_gene_ids <- gsub("gene-Dmel_", "", gene_ids)

# Map CG numbers to gene symbols
gene_mapping <- AnnotationDbi::select(org.Dm.eg.db, 
                       keys = clean_gene_ids,
                       columns = c("SYMBOL", "FLYBASE"),
                       keytype = "FLYBASECG")

# Create a named vector for mapping
symbol_map <- setNames(gene_mapping$SYMBOL, gene_mapping$FLYBASECG)

# Replace rownames with gene symbols (keep original if no symbol found)
new_rownames <- ifelse(is.na(symbol_map[clean_gene_ids]), 
                       gene_ids, 
                       symbol_map[clean_gene_ids])

rownames(em) <- new_rownames

# save em with symbols
write.table(em, 
            file = "em_symbols_update.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = TRUE)

######

library(DESeq2)
library(ggplot2)
library (pheatmap)

# match samples between em and samp
samp <- samp[colnames(em), ]

### NA error - removed "rna-NM_001260299.1" ###
em <- em[!rownames(em) %in% "rna-NM_001260299.1", ]
sum(is.na(em))

### Error in DESeqDataSet(se, design = design, ignoreRank) : some values in assay are not integers ###
em <- round(em)

#############################
# Quickstart - exploratory
dds <- DESeqDataSetFromMatrix(countData = em,
                              colData = samp,
                              design = ~population)
dds

# Run DESeq
dds <- DESeq(dds)

### PCA ###
# variance_stabilised_data <- var_stab_transform(dds, ignore exp design during transformation)
# FOR INTEGRATION
# Create the VST-transformed object
vsd <- vst(dds, blind = TRUE)
# extract matrix for mofa
vst_counts <- assay(vsd)
# save it
write.table(vst_counts, 
            file = "em_symbols_vst.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

pca_data <- plotPCA(vsd, intgroup = c("population","nuclear"), returnData = TRUE) # returnData T = return coordinates also
pca2 <- prcomp(t(assay(vsd)))
summary(pca_data)

# Plot
p1 <- ggplot(pca_data, aes(PC1, PC2, color = population)) +
  geom_point(size = 3, alpha = 0.7) +
  nature_omics_theme +
  ggtitle("PCA by Population")
print(p1)


p2 <- ggplot(pca_data, aes(PC1, PC2, color = nuclear)) +
  geom_point(size = 3, alpha = 0.7) +
  nature_omics_theme +
  ggtitle("PCA by nDNA")
print(p2)
##############################

#
##
###  run diff exp analysis
##
#
# set up with likelihood ratio test, compare ~population (full) model to null model
dds <- DESeq(dds, test="LRT", reduced=~1)
res <- results(dds)
res

# Check - x1 NA
head(res)

table(res$padj <= 0.05)

# collect all padj under .05 and not NA
sig_diff <- (res$padj < 0.05 & !is.na(res$padj)) 
head(sig_diff)

summary(sig_diff) # 17918 NAs

de_genes <- rownames(res)[sig_diff]
length(de_genes)


#Get normalised counts for sig de
# extract from object first
normalized_counts <- counts(dds, normalized=TRUE)

de_norm_matrix <- normalized_counts[de_genes, ]

dim(de_norm_matrix)

### HEATMAP ###
# Assign annotations
col_annotations <- data.frame(
  nDNA = samp$nuclear,
  mtDNA = samp$mito,
  row.names = rownames(samp)
)
######################################################################
####################### NOT USED #####################################
# log2 transform
de_norm_matrix_log <- log2(de_norm_matrix + 1)

# Plot HM
hm = pheatmap(de_norm_matrix_log, 
         annotation_col = col_annotations,
         scale = "row",
         show_rownames = F,
         show_colnames = F, # too many
         cluster_rows = FALSE, # remove LHS dendrogram
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Differentially Expressed Genes Heatmap",
         fontsize = 10)
hm
#####################################################################

#get top 100 
top_100 <- order(res[sig_diff, ]$padj)[1:100] # order() ranks them most sig -least
top_de_genes <- rownames(res[sig_diff, ])[top_100]
top_de_matrix_log <- log2(normalized_counts[top_de_genes, ] + 1)

pheatmap(top_de_matrix_log,
         annotation_col = col_annotations,
         scale = "row",
         show_rownames = FALSE, # top 100
         show_colnames = FALSE,
         clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation",
         border_color = NA, # rm gridlines
         main = "Top 100 Differentially Expressed Genes",
         fontsize = 8,
         fontsize_row = 4,
         height = 17, width = 10)

##### EXTRACT MOFA2 INPUT #####
  ##########################
    #####################
      ################
        ###########
          #######
            ###
             #
#write.csv(de_norm_matrix, "de_genes_normalized_counts.csv")
#write.csv(de_norm_matrix_log, "de_genes_log2_normalized_counts.csv")



### UMAP ###
library(umap)
library(gridExtra)
library(grid) 
vsd_matrix <- assay(vsd)

# transpose for umap expectations (samples as rws)
vsd_t <- t(vsd_matrix)

set.seed(42)
umap_result <- umap(vsd_t)

# Extract umap coordinates
umap_df <- data.frame(
  UMAP1 = umap_result$layout[,1],
  UMAP2 = umap_result$layout[,2],
  sample = rownames(vsd_t)
)

# merge with sample sheet
umap_df <- merge(umap_df, samp, by.x = "sample", by.y = "row.names")

### UMAP plots ###
# UMAP by population
u1 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = population)) +
  geom_point(size = 3, alpha = 0.7) +
  nature_omics_theme +
  ggtitle("UMAP by Population") +
  xlab("UMAP1") +
  ylab("UMAP2") +
  labs(color = "Population")

print(u1)

#UMAP by nuclear
u2 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = nuclear)) +
  geom_point(size = 3, alpha = 0.7) +
  nature_omics_theme +
  ggtitle("UMAP by nDNA") +
  xlab("UMAP1") +
  ylab("UMAP2") +
  labs(color = "nDNA")

print(u2)

#UMAP by mito
u3 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = mito)) +
  geom_point(size = 3, alpha = 0.7) +
  nature_omics_theme +
  ggtitle("UMAP by mtDNA") +
  xlab("UMAP1") +
  ylab("UMAP2") +
  labs(color = "mtDNA")

print(u3)



# UMAP by Mitonuclear (using population)
u4 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = population)) +
  geom_point(size = 3, alpha = 0.7) +
  nature_omics_theme +
  ggtitle("UMAP by Mitonuclear") +
  xlab("UMAP1") +
  ylab("UMAP2") +
  labs(color = "Population")

print(u4)

# Combine mito, nuc, and pop by mitonuc plts
combined_plot <- grid.arrange(
  u4, u3, u2,
  ncol = 3, nrow = 1,
  top = textGrob("Dimensionality Reduction Visualization of Mitonuclear Genotype Effects on Expression Profiles", gp = gpar(fontsize = 16, fontface = "bold"))
)

#############
#PERMANOVA
#############
library(vegan)
# We have already transposed the vsd matrix so that samples are rows 

# Create a distance matrix - bray-curt diss for this 
dist_matrix <- vegdist(vsd_t, method = "bray")

# test for population effect
perm_pop <- adonis2(dist_matrix ~ population, data = samp, permutations = 999) #999 balances statstical precision with comp efficiency - smallest p = 0.001 (1/1000)
print(perm_pop)

# perm addative
perm_additive <- adonis2(dist_matrix ~ nuclear + mito, data = samp, permutations = 999)
print(perm_additive)

# test for interaction effects
perm_mitonuclear <- adonis2(dist_matrix ~ nuclear * mito, data = samp, permutations = 999)
print(perm_mitonuclear)

#### Lets create tables for this data - using our function from the microbiome script ####
### functions have been placed in the OMICS.R script (should be in the same dir and sourced at the top) ###

# Load required library
library(kableExtra)

# Generate tables using the generalized functions
population_results <- create_population_table_general(perm_pop)
mitonuclear_results <- create_mitonuclear_table_general(perm_additive, perm_mitonuclear)

# Create formatted tables
create_formatted_kable(population_results, 
                       "on Gene Expression Profiles (Bray-Curtis dissimilarity)")

create_formatted_kable(mitonuclear_results, 
                       "on Gene Expression Profiles (Bray-Curtis dissimilarity)")



