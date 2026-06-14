rm(list=ls())
#set wd
setwd("C:/Users/sandr/OneDrive/Desktop/mbio")
# link theme file
source("OMICS_mbio.R")

#load packages
#if (!require("pacman")) install.packages("pacman")
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
pacman::p_load(vegan, pheatmap, ggplot2, ape, SummarizedExperiment, gridExtra, grid, kableExtra)


theme_msb <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      plot.title       = element_text(size = base_size, face = "bold", hjust = 0),
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size - 1, color = "black"),
      legend.title     = element_text(size = base_size - 1, face = "bold"),
      legend.text      = element_text(size = base_size - 1),
      legend.key.size  = unit(0.4, "cm"),
      panel.border     = element_blank(),
      axis.line        = element_line(color = "black", linewidth = 0.4),
      axis.ticks       = element_line(color = "black", linewidth = 0.4),
      plot.margin      = margin(8, 8, 8, 8),
      strip.background = element_blank(),
      strip.text       = element_text(size = base_size - 1, face = "bold")
    )
}

geno_colors_msb <- c("AA" = "#E69F00", "AB" = "#CC4444",
                     "BA" = "#F0B8B8", "BB" = "#7B4F9E")

# Set out dir
OUT_DIR <- "C:/Users/sandr/OneDrive/Desktop/mbio"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

#
##
###
# set the metadata making any strings factors 
ind <- read.csv("samples.csv", stringsAsFactors=T)
# make a new col for mitonuclear by concatenating the mito and nuc genotype identifiers
# factor() converts the resulting concat strings to factor variables
ind$mitonuclear <- with(ind, factor(paste(mito, nuclear, sep="")))

# sort dodgy naming of sample col.....(i..samp)
colnames(ind)[1] <- "samp"


# load in the data at your chosen filtering threshold - we went with 1% abund in 10% samps
# at least 1percent abundance in at least 10percent of samples - it offered the clearest resolution.
d <- read.table("noWD_featuretable_L5_1_10.txt", 
                header = TRUE,
                sep = "\t",
                skip = 1,
                fill = TRUE,  
                check.names = FALSE)

# add "DS" prefix to column names (skip first column which is "OTU ID")
colnames(d)[2:ncol(d)] <- paste0("DS", colnames(d)[2:ncol(d)])

head(d)
dim(d)
###
##
#

#check correspondence
which(ind$samp %in% colnames(d))

#remove lines of index for which we have no data
nrow(ind)
ind <- droplevels(subset(ind, samp %in% colnames(d)))
nrow(ind)

#remove dodgy samples (run heatmap without this line to see)
ind <- droplevels(subset(ind, !samp %in% c("DS37", "DS2")))

#turn the data into just a matrix
rowInd <- d[,c(1,109)]
head(rowInd)
d <- d[,2:108]
dim(d)

# Assign rownames
#rownames(d) <- rowInd$OTU
rownames(d) <- rowInd$`OTU ID`

rownames(ind) <- ind$samp

#line up the columns of the matrix to the metadata
d <- d[,match(ind$samp, colnames(d))]
colnames(d)
#####
####
###
#

# log2 transform
relAb <- t(apply(d, 1, y=colSums(d), function(x, y){x/y}))

##################################################
## ASV TAXON STACKED BARPLOT for Fig 1
## relAb is already calculated 
##################################################
library(ggplot2)
library(dplyr)
library(tidyr)

# Make df
relAb_df <- as.data.frame(relAb)
relAb_df$taxon <- rowInd$`OTU ID`  # ` ` needed due to space in col name

relAb_long <- relAb_df %>%
  pivot_longer(-taxon, names_to = "samp", values_to = "rel_ab") %>%
  left_join(ind[, c("samp", "mitonuclear")], by = "samp")

#Order samples by mitonuclear genotype 
sample_order <- ind %>%
  arrange(mitonuclear, samp) %>%
  pull(samp)

relAb_long$samp <- factor(relAb_long$samp, levels = sample_order)

# Colour palette for 12 taxa 
taxa_list <- unique(relAb_long$taxon)

taxon_colours <- setNames(
  c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7", "#999999", "#44AA99", "#882255",
    "#DDCC77", "#117733"),
  taxa_list
)

relAb_long$taxon <- factor(relAb_long$taxon, levels = taxa_list)

#Genotype group divider positions 
group_positions <- relAb_long %>%
  distinct(samp, mitonuclear) %>%
  arrange(samp) %>%
  group_by(mitonuclear) %>%
  summarise(
    xmin = min(as.numeric(samp)),
    xmax = max(as.numeric(samp)),
    xmid = mean(as.numeric(samp))
  )

# Plotit
taxon_barplot <- ggplot(relAb_long,
                        aes(x = samp, y = rel_ab, fill = taxon)) +
  geom_bar(stat = "identity", width = 1.0, colour = NA) +
  
  geom_vline(data = group_positions[-1, ],
             aes(xintercept = xmin - 0.5),
             colour = "black", linewidth = 0.5) +
  
  annotate("text",
           x     = group_positions$xmid,
           y     = 1.04,
           label = group_positions$mitonuclear,
           size  = 4, fontface = "bold", vjust = 0) +
  
  scale_fill_manual(values = taxon_colours, name = "Family") +
  
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0, 0),
                     limits = c(0, 1.10)) +
  
  labs(x     = "Sample (by mitonuclear genotype)",
       y     = "Relative Abundance",
       title = "A  Gut Microbiota Composition by Mitonuclear Genotype") +
  
  theme_msb() +
  theme(
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    plot.title      = element_text(hjust = 0, size = 11, face = "bold"),
    legend.text     = element_text(size = 8, face = "italic"),
    legend.title    = element_text(size = 9, face = "bold"),
    legend.key.size = unit(0.4, "cm")
  )

print(taxon_barplot)

ggsave(file.path(OUT_DIR, "Fig1A_taxon_barplot_mitonuclear.pdf"), taxon_barplot,
       width = 12, height = 5, units = "in")
ggsave(file.path(OUT_DIR, "Fig1A_taxon_barplot_mitonuclear.png"), taxon_barplot,
       width = 12, height = 5, units = "in", dpi = 300)

#####
####
###
##
#
##################################################
## SHANNON DIVERSITY PLOT for Fig 1
##################################################
library(vegan)

# calculate shan on raw counts (before log transform)
shannon <- diversity(t(d), index = "shannon")

# Build df
shannon_df <- data.frame(
  samp        = ind$samp,
  shannon     = shannon,
  mitonuclear = ind$mitonuclear,
  mito        = ind$mito
)

# Wilcoxon test mito A vs B for annotation
shannon_wilcox <- wilcox.test(shannon ~ ind$mito)
shannon_p <- round(shannon_wilcox$p.value, 3)
p_label <- ifelse(shannon_p < 0.001, "p < 0.001", paste0("p = ", shannon_p))# mito A vs B (Wilcoxon)

library(ggplot2)
library(ggpubr)


shannon_plot <- ggplot(shannon_df, aes(x = mitonuclear, y = shannon,
                                       fill = mitonuclear, colour = mitonuclear)) +
  ##############################################################
  #### Shaded background bands to visually group mito A vs B ###
  ##############################################################
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = geno_colors_msb[["AA"]], alpha = 0.06) +
  annotate("rect", xmin = 2.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
           fill = geno_colors_msb[["BA"]], alpha = 0.06) +
  
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.45,
               linewidth = 0.4, colour = "grey30") +
  geom_jitter(shape = 19, width = 0.15, size = 1.8, alpha = 0.85,
              stroke = 0.3, colour = "grey20") +
  
  # Mito group brackets
  geom_bracket(xmin = "AA", xmax = "AB", y.position = 2.35,
               label = "Mitochondrial A", inherit.aes = FALSE,
               label.size = 2.8, tip.length = 0.02, linewidth = 0.4) +
  geom_bracket(xmin = "BA", xmax = "BB", y.position = 2.35,
               label = "Mitochondrial B", inherit.aes = FALSE,
               label.size = 2.8, tip.length = 0.02, linewidth = 0.4) +
  
  # p-value from Wilcoxon mito A vs B
  annotate("text", x = 2.5, y = 2.62, 
           label = p_label, size = 2.8, hjust = 0.5) +
  
  scale_fill_manual(values = geno_colors_msb, name = "Genotype") +
  scale_colour_manual(values = geno_colors_msb, name = "Genotype") +
  scale_y_continuous(limits = c(NA, 2.7), expand = c(0.02, 0)) +
  
  labs(x = "Mitonuclear genotype",
       y = "Shannon diversity index",
       title = "B") +
  
  theme_msb() +
  theme(
    legend.position  = "none",
    axis.ticks.length = unit(2, "mm"),
    plot.margin      = margin(8, 8, 8, 8, "mm")
  )

print(shannon_plot)


ggsave(file.path(OUT_DIR, "Fig1B_shannon_mitonuclear.pdf"), shannon_plot,
       width = 90, height = 80, units = "mm")
#ggsave(file.path(OUT_DIR, "Fig1B_shannon_mitonuclear.png"), shannon_plot,
       #width = 90, height = 80, units = "mm", dpi = 300)


##################################################
## PCoA PLOT Fig 1
##################################################

# BrayCurtis on raw counts... before log transform
dDist <- vegdist(t(d), method = "bray")
PcoA  <- cmdscale(dDist, eig = TRUE)  # eig=TRUE to extract % variance explained

# Extract % variance explained for axis labels
eig_vals  <- PcoA$eig
pct_var   <- round(eig_vals / sum(eig_vals[eig_vals > 0]) * 100, 1)

PcoDF <- data.frame(
  ind,
  Dim1 = PcoA$points[, 1],
  Dim2 = PcoA$points[, 2]
)


pcoa_plot <- ggplot(PcoDF, aes(x = Dim1, y = Dim2, 
                               colour = mitonuclear, fill = mitonuclear)) +
  stat_ellipse(geom = "polygon", alpha = 0.08, linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_colour_manual(values = geno_colors_msb, name = "Genotype") +
  scale_fill_manual(values = geno_colors_msb, name = "Genotype") +
  labs(x     = paste0("PCoA Axis 1 (", pct_var[1], "% variance)"),
       y     = paste0("PCoA Axis 2 (", pct_var[2], "% variance)"),
       title = "C  Gut Microbiota Community Structure by Mitonuclear Genotype") +
  theme_msb() +
  theme(
    plot.title      = element_text(hjust = 0, size = 11, face = "bold")
  )

print(pcoa_plot)
ggsave(file.path(OUT_DIR, "Fig1C_pcoa_mitonuclear.pdf"), pcoa_plot, width = 7, height = 5, units = "in")
ggsave(file.path(OUT_DIR, "Fig1C_pcoa_mitonuclear.png"), pcoa_plot, width = 7, height = 5, units = "in", dpi = 300)
#
##
###
####
#####




### PERM ###
#test for mito-nuclear effects in microbiota
# make this into a table on it's own, then 
set.seed(1984)
perm1 <- adonis2(t(d) ~ population, data=ind, by="terms", method = "bray") # I used bray-curtis matrix here to match perm updates as this is optimal
perm1	#the microbiota do differ among these populations (this cmmt for 1_10)

perm2 <- adonis2(t(d) ~ mito * nuclear, data=ind, by="terms")
perm2	#there is no mito-nuclear interaction (1_10)

#simplify model to see if there are independent effects of mito and nuc
perm3 <- adonis2(t(d) ~ mito + nuclear, data=ind, by="terms")
perm3	#no. there is only an effect of mitochondria. (1_10)

###
######
########

# PERM UPDATE
# PERMANOVA Analysis for Microbiota Community Structure
# Using 1% in 10% samples filtering level
#
library(gridExtra)
library(kableExtra)

set.seed(1984)
# pop is perm1 above
# Mitonuclear interaction
perm_mitonuclear <- adonis2(t(d) ~ mito * nuclear, data=ind, by="terms", method="bray")
print(perm_mitonuclear)

# Mito + nuclear
perm_additive <- adonis2(t(d) ~ mito + nuclear, data=ind, by="terms", method="bray")
print(perm_additive)

### Create summary tables ###
# Population effect

#
create_population_table <- function() {
  # pull our relevant stats from the analysis
  pop_r2 <- round(perm1$R2[1], 3)  # Using perm1 from initial
  pop_f <- round(perm1$F[1], 2)
  pop_p <- perm1$`Pr(>F)`[1]
  # reusable func to assign sig groups for <0.001, <0.01, <0.05
  format_p <- function(p) {
    if (p < 0.001) return("< 0.001")
    else if (p < 0.01) return("< 0.01")
    else if (p < 0.05) return("< 0.05")
    else return(sprintf("%.3f", p)) # round not sig values to 3 sig figs
  }
  # Create a df for table generation
  population_table <- data.frame(
    Factor = "Population",
    R2 = pop_r2,
    F_statistic = pop_f,
    P_value = format_p(pop_p),
    # set sig constraint visual as in adonis2 output as I will set a formatted legend for final table
    Significance = ifelse(pop_p < 0.001, "***", ifelse(pop_p < 0.01, "**", ifelse(pop_p < 0.05, "*", "ns")))
  )
  
  return(population_table)
}

# same again for mitonuclear effects (without pop) using the new analysis
create_mitonuclear_table <- function() {
  mito_r2 <- round(perm_additive$R2[1], 3)
  mito_f <- round(perm_additive$F[1], 2)
  mito_p <- perm_additive$`Pr(>F)`[1]
  
  nuclear_r2 <- round(perm_additive$R2[2], 3)
  nuclear_f <- round(perm_additive$F[2], 2)
  nuclear_p <- perm_additive$`Pr(>F)`[2]
  
  interaction_r2 <- round(perm_mitonuclear$R2[3], 3)
  interaction_f <- round(perm_mitonuclear$F[3], 2)
  interaction_p <- perm_mitonuclear$`Pr(>F)`[3]
  
  format_p <- function(p) {
    if (p < 0.001) return("< 0.001")
    else if (p < 0.01) return("< 0.01")
    else if (p < 0.05) return("< 0.05")
    else return(sprintf("%.3f", p))
  }
  
  mitonuclear_table <- data.frame(
   # Written as "Mito x Nuclear" directly here!!! #
    Factor = c("Mitochondrial", "Nuclear", "Mito x Nuclear"),
    R2 = c(mito_r2, nuclear_r2, interaction_r2),
    F_statistic = c(mito_f, nuclear_f, interaction_f),
    P_value = c(format_p(mito_p), format_p(nuclear_p), format_p(interaction_p)),
    # sig constraints
    Significance = c(
      ifelse(mito_p < 0.001, "***", ifelse(mito_p < 0.01, "**", ifelse(mito_p < 0.05, "*", "ns"))),
      ifelse(nuclear_p < 0.001, "***", ifelse(nuclear_p < 0.01, "**", ifelse(nuclear_p < 0.05, "*", "ns"))),
      ifelse(interaction_p < 0.001, "***", ifelse(interaction_p < 0.01, "**", ifelse(interaction_p < 0.05, "*", "ns")))
    )
  )
  
  return(mitonuclear_table)
}

# Generate both tables
population_results <- create_population_table()
mitonuclear_results <- create_mitonuclear_table()

options(knitr.table.format = "simple")

# Formatted with kable
# formatted as "R2" here due to ?.
kable(population_results,
      format = "html",
      caption = "PERMANOVA Results Testing Population Effects on Microbiota Community Structure (Bray-Curtis dissimilarity)",
      col.names = c("Factor", "R2", "F", "P-value", "Sig."),
      align = c('l', 'c', 'c', 'c', 'c')) %>% # set the factor name to LHS and centre data cols for separation
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "bordered"), 
                full_width = FALSE,
                position = "center") %>%
  footnote(general = "Significance codes: *** P < 0.001, ** P < 0.01, * P < 0.05, ns = not significant")


kable(mitonuclear_results,
      format = "html",
      caption = "PERMANOVA Results Testing Mitonuclear Effects on Microbiota Community Structure (Bray-Curtis dissimilarity)", 
      col.names = c("Factor", "R2", "F", "P-value", "Sig."),
      align = c('l', 'c', 'c', 'c', 'c')) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "bordered"), 
                full_width = FALSE,
                position = "center") %>%
  footnote(general = "Significance codes: *** P < 0.001, ** P < 0.01, * P < 0.05, ns = not significant")


###############################################
## PAPER FIGURES - PERMANOVA TABLES (Fig1D/E)
###############################################

pop_grob <- tableGrob(
  population_results,
  rows  = NULL,
  theme = ttheme_minimal(
    base_size   = 9,
    base_family = "sans",
    core    = list(fg_params = list(hjust = 0.5, x = 0.5)),
    colhead = list(fg_params = list(hjust = 0.5, x = 0.5, fontface = "bold"))
  )
)

pop_titled <- arrangeGrob(
  textGrob("D  PERMANOVA: Microbiota community structure ~ population (Bray-Curtis)",
           x = 0, hjust = 0,
           gp = gpar(fontsize = 9, fontface = "bold", fontfamily = "sans")),
  pop_grob,
  nrow = 2, heights = c(0.15, 1)
)

ggsave(file.path(OUT_DIR, "Fig1D_permanova_population.pdf"),
       plot = pop_titled, width = 16, height = 4,
       units = "cm", device = "pdf", useDingbats = FALSE)


mito_grob <- tableGrob(
  mitonuclear_results,
  rows  = NULL,
  theme = ttheme_minimal(
    base_size   = 9,
    base_family = "sans",
    core    = list(fg_params = list(hjust = 0.5, x = 0.5)),
    colhead = list(fg_params = list(hjust = 0.5, x = 0.5, fontface = "bold"))
  )
)

mito_titled <- arrangeGrob(
  textGrob("E  PERMANOVA: Microbiota community structure ~ mitonuclear genotype (Bray-Curtis)",
           x = 0, hjust = 0,
           gp = gpar(fontsize = 9, fontface = "bold", fontfamily = "sans")),
  mito_grob,
  nrow = 2, heights = c(0.15, 1)
)

ggsave(file.path(OUT_DIR, "Fig1E_permanova_mitonuclear.pdf"),
       plot = mito_titled, width = 16, height = 5,
       units = "cm", device = "pdf", useDingbats = FALSE)

message("All Figure 1 outputs saved to: ", OUT_DIR)

########
#####
###


###################################################################
##################### NOT USED ####################################
#pheatmap(log(d+1),
 #        annotation_col=ind[,c("mito", "nuclear", "population")],
   #      cluster_cols=F,
     #    labels_row=rowInd$OTU.ID#,
         #	scale="row"
#)

###################################################################
############ PCA NOT USED ########################################
#PCA
pca1 <- prcomp(t(log2(d+1)), scale=T, center=T)
summary(pca1)
PCADF <- data.frame(ind, pca1$x)
p1 <- ggplot(data=PCADF, aes(x=PC1, y=PC2, col=mitonuclear))
p1 + 
  geom_point() +
  facet_grid(nuclear ~ mito) +
  theme_bw()
###################################################################

#PCoA on LOG microbiome reads
dDist <- vegdist(as.matrix(t(log2(d+1))), method = "bray") # Bray-Curtis dissimilarity matrix used in accordance with perm
PcoA <- cmdscale(dDist) # cmdscale = classical multidimensional scaling
PcoDF <- data.frame(ind, Dim1=PcoA[,1], Dim2=PcoA[,2]) # top 2 axes that explain the most variation
# Plot it...
p2 <- ggplot(data=PcoDF, aes(x=Dim1, y=Dim2, col=mitonuclear))
p2 + 
  geom_point() +
  facet_grid(nuclear ~ mito) +
  theme_bw() +
  scale_y_continuous(expand = c(0,0))


# Corresponding box plots
p3 <- ggplot(data=PcoDF, aes(y=Dim1, x=mitonuclear, col=mitonuclear))
p3 + 
  geom_boxplot(outliers=F) +
  geom_jitter() +
  nature_omics_theme

p4 <- ggplot(data=PcoDF, aes(y=Dim2, x=mitonuclear, col=mitonuclear))
p4 + 
  geom_boxplot(outliers=F) +
  geom_jitter() +
  nature_omics_theme

##
### PCoA fig update
##
# Calculate diversity indices
shannon <- diversity(t(d), index="shannon")
richness <- specnumber(t(d))

### set margins
par(mar=c(5,4,4,2))

layout(matrix(c(1,2,3,4), nrow=2, byrow=TRUE), heights = c(1,1,0.5,0.5))

# Define mitonuclear grouping colours
mitonuclear_colors <- c("AA" = "red", "AB" = "darkorchid1", "BA" = "orange", "BB" = "blue")

# Calculate sig
shannon_t_test <- t.test(shannon ~ ind$mito)
p_value <- round(shannon_t_test$p.value, 3)
print(shannon_t_test)

# Shannon diversity by mito
colors_mito <- c("red", "blue")
plot(shannon, col=colors_mito[as.factor(ind$mito)], pch=19, cex=1.2,
     xlab="", ylab="Shannon Diversity", 
     main="Shannon Diversity by Mito",
     ylim=c(0, ceiling(max(shannon))),
     xaxt="n")
mtext(paste0("p = ", p_value), side=3, line=0.5, cex=0.9)
par(xpd=TRUE)
legend("right", inset=c(-0.10, 0),
       legend=levels(as.factor(ind$mito)), 
       col=colors_mito, pch=19, cex=0.9, bty="n")
par(xpd=FALSE)

# PCoA re-plot
plot(PcoDF$Dim1, PcoDF$Dim2, 
     col=mitonuclear_colors[as.factor(PcoDF$mitonuclear)], 
     pch=19,
     xlab="PCoA Dimension 1", 
     ylab="PCoA Dimension 2",
     ylim=range(PcoDF$Dim2),
     main="PCoA")
par(xpd=TRUE)
legend("right", inset=c(-0.10, 0),
       legend=names(mitonuclear_colors), 
       col=mitonuclear_colors, pch=19, cex=0.8, bty="n")
par(xpd=FALSE)

# PCoA Dim1 boxplot
boxplot(Dim1 ~ mitonuclear, data=PcoDF, 
        col=mitonuclear_colors[levels(PcoDF$mitonuclear)],
        main="PCoA Dim1",
        ylab="PCoA Dimension 1",
        las=2)
par(xpd=TRUE)
legend("right", inset=c(-0.10, 0),
       legend=names(mitonuclear_colors),
       col=mitonuclear_colors, pch=15, cex=0.7, bty="n")
par(xpd=FALSE)

# PCoA Dim2 boxplot
boxplot(Dim2 ~ mitonuclear, data=PcoDF, 
        col=mitonuclear_colors[levels(PcoDF$mitonuclear)],
        main="PCoA Dim2", 
        ylab="PCoA Dimension 2",
        las=2)
par(xpd=TRUE)
legend("right", inset=c(-0.10, 0),
       legend=names(mitonuclear_colors),
       col=mitonuclear_colors, pch=15, cex=0.7, bty="n")
par(xpd=FALSE)

# Reset layout
layout(1)

##
### Calculate relative abundances
##
### updated
##

### calculate relative abundances
# log2 transform
#relAb <- t(apply(d, 1, y=colSums(d), function(x, y){x/y}))

d <- log2(d + 1)


# Collect p-values to plot each ASV individually
# make empty vectors to store
p_values <- c()
taxa_names <- c()

par(mfrow=c(1,1))
par(bty="n")


# Continue with coefficient of variation analysis...
# Coefficient of variation = variance / mean
coefVar <- apply(d, 1, function(x){var(x)/mean(x)})
par(mar=c(12, 4, 4, 2))
barplot(coefVar, las=2)


# Set up multi-panel plot - adjust based on number of taxa
n_taxa <- nrow(relAb)
n_cols <- ceiling(sqrt(n_taxa))
n_rows <- ceiling(n_taxa / n_cols)
par(mfrow=c(3,4))
#par(mar=c(5,3,2,2))
par(mar=c(4,4.5,4.5,1)) # the right one; had to disable to fix prev plot
par(bty="n")


for(i in 1:nrow(relAb)){
  #ASV <- rowInd$OTU.ID[i]
  ASV <- rowInd$`OTU ID`[i]
  y <- as.numeric(relAb[i,])
  x3 <- ind$mito
  
  # box plots (i) comparing mitoA vs mitoB
  plot(y ~ x3, main=paste(ASV), las=2, 
       xlab="mtDNA", ylab="Relative Abundance", outline=F,
       col=c("red", "blue"), cex.main=0.9, cex.axis=0.8)
  
  y_A <- y[ind$mito=="A"]
  y_B <- y[ind$mito=="B"]
  
  # Wilcoxon test and collect results
  test_result <- wilcox.test(y_A, y_B)
  p_values <- c(p_values, test_result$p.value)
  taxa_names <- c(taxa_names, ASV)
  
  # Add p-value to plot (position under the titles w/ mtext)
  if(test_result$p.value < 0.001) {
    p_text <- "p < 0.001"
  } else {
    p_text <- paste("p =", round(test_result$p.value, 3))
  }
  #mtext(p_text, side=3, line=-2, cex=0.7)
  #mtext(p_text, side=3, line=0.5, cex=0.6)
  mtext(p_text, side=3, line=-0.1, cex=0.7)
  #text(x=1.5, y=max(y)*0.85, labels=p_text, cex=0.7, font=1)
}

# Reset layout
par(mfrow=c(1,1))
par(mar=c(5,4.8,4,2)) 
#par(mar=c(3.5,3.5,3.5,1))

# Apply FDR correction post loop
adjusted_p <- p.adjust(p_values, method="fdr")
results_df <- data.frame(Taxa=taxa_names, P_value=p_values, FDR=adjusted_p)
significant <- results_df[results_df$FDR < 0.05,]
