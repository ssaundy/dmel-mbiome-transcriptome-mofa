setwd("C:/Users/sandr/OneDrive/Desktop/mofa")
source("OMICS.R")
rm(list=ls())
library(ggplot2)
library(MOFA2)
library(pheatmap)

# Assign colours to our genotypes
geno_colors <- c("AA"="orange", "AB"="red", "BA"="pink", "BB"="purple")

# load model
model <- load_model("model.hdf5")
plot_data_overview(model)

Nsamples = sum(model@dimensions$N)
samps <- samples_names(model)[[1]]

# Add metadata to the model
sample_metadata <- read.csv("samples.csv") # columns = sample, pop, mito, nuclear, replicate

# Fix weird sample col name
colnames(sample_metadata)[1] <- "samp"

# Check what columns we actually have
print(colnames(sample_metadata))
print(head(sample_metadata))

# Create mitonuclear factor
sample_metadata$mitonuclear <- with(sample_metadata, factor(paste(mito, nuclear, sep=""))) # factor() ensures discrete factor treatment, rather than continuous, therefore the output can bu used in our statistical models

# Subset
sample_metadata <- droplevels(subset(sample_metadata, samp %in% samps))
sample_metadata <- sample_metadata[match(samps, sample_metadata$samp),]

# rename the column to 'sample' for MOFA compatibility
colnames(sample_metadata)[1] <- "sample"

# Add to model
samples_metadata(model) <- sample_metadata

# Check genotype distribution
table(sample_metadata$mitonuclear)



###############################################
## 1. VARIANCE EXPLAINED - Answer "How much?"
###############################################
############# NOT USED ########################
# Total variance explained (exploratory)
#plot_variance_explained(model, x="view", y="factor", factors=1:15) #
## For paper, 15 too much
plot_variance_explained(model, x="view", y="factor", factors=1:10,
                        add_text = TRUE) +
  theme_classic() +
  labs(title = "Variance Explained") +
  theme(strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5))
###############################################

### UPDATED VARIANCE EXPLAINED PLOT ###

### Get the variance data in a matrix ###
# get_variance_explained(),mofa func- extracts variance statistics from the model
# r2_per_factor being a metric returned by get_variance_explained() 
# group1 being the default first group name, samples are all in one group
# nested inside group1 are vals for each view (fly,micro) for the number of factor specified
var_data <- get_variance_explained(model)$r2_per_factor$group1[1:10, ] 
## $r2_total() also an option... total variance across all factors

# Create the plot and add text
p <- plot_variance_explained(model, x="view", y="factor", factors=1:10) +
  theme_classic() +
  labs(title = "Variance Explained") +
  theme(strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5))

# Add the text values
p + geom_text(data = expand.grid(factor = paste0("Factor", 1:10), 
                                 view = c("fly", "micro")),
              aes(x = view, y = factor, 
                  label = round(as.vector(var_data), 1)), # whole matrix as a single vector
              inherit.aes = FALSE, color = "firebrick", size = 3)


# Simple bar
plot_variance_explained(model, x="group", y="factor", plot_total = TRUE)

# Simple bar with group labels removed
plot_variance_explained(model, x="group", y="factor", plot_total = TRUE) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())



###############################################
## 2. GENOTYPE EFFECTS ON FACTORS 
## Answer Aim 1:	How does host genetic variation (mitochondrial and nuclear) affect gut microbiome composition?
###############################################

# factor values by genotype (violin plot)
plot_factor(model, 
            color_by="mitonuclear",
            factors = 1:6,  # Focus on top factors
            add_violin = TRUE,
            dodge = TRUE
) + scale_fill_manual(values=geno_colors) +
  scale_color_manual(values=geno_colors)

# Statistical test: Do genotypes differ significantly in factor values?
# first get our factor num and corresponding values as columns in the sample metadata sheet
factor_values <- get_factors(model, as.data.frame=TRUE) # using get_factors() to extract latent factors
factor_values_with_meta <- merge(factor_values, sample_metadata, by="sample")
factor_values_with_meta
# calc p values for each factor
unique_factors <- unique(factor_values_with_meta$factor)
results <- data.frame(factor = unique_factors, p_value = NA)

for(i in 1:length(unique_factors)) {
  factor_data <- subset(factor_values_with_meta, factor == unique_factors[i])
  # 1 way anova testing testing if factor vals are sig different between genotype groups (AA, AB, BA, BB)
  aov_result <- aov(value ~ mitonuclear, data = factor_data)
  results$p_value[i] <- summary(aov_result)[[1]][["Pr(>F)"]][1] # [[1]][["Pr(>F)"]][1] to extract p: https://stackoverflow.com/questions/3366506/extract-p-value-from-aov
}

print(results)

################################
### MAKE  REASONIG TABLE ###
################################
# create the table
factor_table <- data.frame(factor = unique_factors, p_value = NA, f_statistic = NA)

for(i in 1:length(unique_factors)) {
  factor_data <- subset(factor_values_with_meta, factor == unique_factors[i])
  aov_result <- aov(value ~ mitonuclear, data = factor_data)
  aov_summary <- summary(aov_result)[[1]]
  # We need to  extract p and F-statistics
  factor_table$p_value[i] <- aov_summary[["Pr(>F)"]][1]
  factor_table$f_statistic[i] <- aov_summary[["F value"]][1]
}

# Sort by p-value (ascending)
factor_table <- factor_table[order(factor_table$p_value), ]

# Add significance codes
factor_table$significance <- ifelse(factor_table$p_value < 0.001, "***",
                                    ifelse(factor_table$p_value < 0.01, "**",
                                           ifelse(factor_table$p_value < 0.05, "*", "ns")))

# Clean up factor names (remove "Factor" prefix for cleaner look)
factor_table$factor <- gsub("Factor", "", factor_table$factor)

# round values for display
factor_table$f_statistic <- round(factor_table$f_statistic, 2)
# If highly sig (<0.001) format it, all else will follow above
factor_table$p_value_display <- ifelse(factor_table$p_value < 0.001, 
                                       format(factor_table$p_value, scientific = TRUE, digits = 2),
                                       round(factor_table$p_value, 3))

# Create final table for display
final_table <- data.frame(
  Factor = factor_table$factor,
  `F-statistic` = factor_table$f_statistic,
  `P-value` = factor_table$p_value_display,
  Significance = factor_table$significance
)


# Print the table
print(final_table)

### CREATE TABLE
### And cut it to 10F
create_formatted_kable_mofa(mofa_results_simple[1:10,], "- ANOVA testing genotype effects on factor values")

###############################################
## 3. FEATURE WEIGHTS - Answer "What drives these factors?" Aim 1
###############################################

# Extract top weighted features for each significant factor (rank by magnitude regardless of positive/negative direction)
get_top_features <- function(model, factor_num, view_name, n_features = 20) {
  weights <- get_weights(model, 
                         views = view_name, 
                         factors = factor_num, 
                         as.data.frame = TRUE,
                         scale = TRUE)
  weights$abs_value <- abs(weights$value)
  top_features <- weights[order(weights$abs_value, decreasing = TRUE), ][1:n_features, ]
  return(top_features)
}

# Analyze factors 1-3
for(factor_num in 1:3) {
  # Plot weights
  p1 <- plot_weights(model,
                     view = "micro",  
                     factor = factor_num,
                     nfeatures = 15,
                     scale = TRUE,
                     abs = FALSE
  ) + ggtitle(paste("Factor", factor_num, "- Microbiome Features")) +
    theme(plot.margin = margin(t = 10, r = 20, b = 10, l = 20),
          axis.text.y = element_text(size = 9, margin = margin(r = 5)),
          panel.grid.major.y = element_line(size = 0.3, color = "grey90"))
  
  p2 <- plot_weights(model,
                     view = "fly", 
                     factor = factor_num,
                     nfeatures = 10,
                     scale = TRUE,
                     abs = FALSE
  ) + ggtitle(paste("Factor", factor_num, "- Transcriptome Features")) +
    theme(plot.margin = margin(t = 10, r = 20, b = 10, l = 20),
          axis.text.y = element_text(size = 9, margin = margin(r = 5)),
          panel.grid.major.y = element_line(size = 0.3, color = "grey90"))
  
  print(p1)
  print(p2)
  ############### These plots are guff, just view the below - but still... ALL MUST BE RUN ###########################

  
  # Get top features
  top_micro <- get_top_features(model, factor_num, "micro", 10) # all factor 3!
  top_fly <- get_top_features(model, factor_num, "fly", 10) # "
}



###############################################
## MICROBIOME-TRANSCRIPTOME CORRELATIONS - Aim 2
###############################################

# Function to create correlation matrix between views
create_cross_correlations <- function(model, factor_num, n_features = 15) {
  # Get top features for each omics type for a factor
  micro_weights <- get_top_features(model, factor_num, "micro", n_features)
  fly_weights <- get_top_features(model, factor_num, "fly", n_features)
  
  # Get the actual data (abundance/expression)
  micro_data <- get_data(model, views = "micro", as.data.frame = TRUE)
  fly_data <- get_data(model, views = "fly", as.data.frame = TRUE)
  
  # Convert feature columns for consistent string matching
  micro_data$feature <- as.character(micro_data$feature)
  fly_data$feature <- as.character(fly_data$feature)
  micro_weights$feature <- as.character(micro_weights$feature)
  fly_weights$feature <- as.character(fly_weights$feature)
  
  # Create empty correlation matrix
  cor_matrix <- matrix(nrow = nrow(micro_weights), 
                       ncol = nrow(fly_weights))
  
  # Fill the corr matrix
  for(i in 1:nrow(micro_weights)) {
    micro_feature <- micro_weights$feature[i]
    micro_vals <- micro_data[micro_data$feature == micro_feature, "value"]
    
    for(j in 1:nrow(fly_weights)) {
      fly_feature <- fly_weights$feature[j]
      fly_vals <- fly_data[fly_data$feature == fly_feature, "value"]
      
      # calculate corr if we have matching data, use = "complete.obs" handles missing vals by only using complete pairs
      if(length(micro_vals) == length(fly_vals) && length(micro_vals) > 0) {
        cor_matrix[i,j] <- cor(micro_vals, fly_vals, use = "complete.obs")
      }
    }
  }
  
  # Add row and col names
  rownames(cor_matrix) <- micro_weights$feature
  colnames(cor_matrix) <- fly_weights$feature
  
  return(cor_matrix)
}

### Run for our first * factors
for(factor_num in 1:10) {
  
  # Create the correlation matrix for this factor
  cor_matrix <- create_cross_correlations(model, factor_num, 5) # top 5 features per omics type - therefore 5x5 matrix
  
  # Make heatmap
  pheatmap(cor_matrix,
           main = paste("Factor", factor_num, ": Microbiome-Transcriptome Correlations"),
           display_numbers = TRUE,
           number_format = "%.2f",
           fontsize_number = 8)
}

#############################
## SUMMARY FOR FACTOR 3 ##
#############################
# Think of a multi plot to show all the views that we have for a given factor (different ways of looking at F to answer our questions)
# show our violin plot, weights, and HM all in one
library(gridExtra)  # for combining panels
library(grid) # for text annotations

# Panel A: Factor activity by genotype
panel_a <- plot_factor(model, 
                       factors = 3,
                       color_by = "mitonuclear",
                       add_violin = TRUE,
                       dodge = TRUE) + 
  scale_fill_manual(values = geno_colors) +
  scale_color_manual(values = geno_colors) +
  ggtitle("Factor 3 Activity by Genotype") +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5), size=12)

# Panel B: Microbiome weights
panel_b <- ggplot(top_micro, aes(x = reorder(feature, abs_value), y = value)) +
  geom_col(aes(fill = value > 0), width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "blue")) +
  coord_flip() +
  labs(title = "Top Microbiome Features",
       size=12,
       x = "Features", 
       y = "Weight") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 8))

# Panel C: Transcriptome weights
panel_c <- ggplot(top_fly, aes(x = reorder(feature, abs_value), y = value)) +
  geom_col(aes(fill = value > 0), width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "blue")) +
  coord_flip() +
  labs(title = "Top Transcriptome Features",
       size = 12,
       x = "Features", 
       y = "Weight") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 8))

# Panel D: Correlations
cor_matrix_f3 <- create_cross_correlations(model, 3, 5) # F3 5x5

# Convert correlation matrix to ggplot-friendly format
cor_df <- expand.grid(Microbiome = rownames(cor_matrix_f3),
                      Transcriptome = colnames(cor_matrix_f3))
cor_df$Correlation <- as.vector(cor_matrix_f3)

panel_d <- ggplot(cor_df, aes(x = Microbiome, y = Transcriptome, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", 
                       midpoint = 0, limits = c(-1, 1)) +
  geom_text(aes(label = round(Correlation, 2)), size = 3) +
  labs(title = "Microbiome-Transcriptome Correlations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        legend.key.size = unit(0.2, "cm"),      # Make legend keys smaller
        legend.key.width = unit(0.1, "cm"),     # Make legend narrower
        legend.title = element_text(size = 8),   # Smaller legend title
        legend.text = element_text(size = 7),    # Smaller legend text
        legend.position = c(1, 0.21))

# Combine all panels
combined_plot <- grid.arrange(
  panel_a, panel_b,
  panel_c, panel_d,
  ncol = 2, nrow = 2,
  top = textGrob("Factor 3 Analysis", gp = gpar(fontsize = 16, fontface = "bold"))
)

###########################
#Show linear relationships
############################

# Clean feature names, trim omics suffixes
current_features <- features_names(model)
current_features$micro <- gsub("_micro", "", current_features$micro)
current_features$fly <- gsub("_fly", "", current_features$fly)
features_names(model) <- current_features

library(car)

# linearise microbime
panel_e <- plot_data_scatter(model,
                             view = "micro",         
                             factor = 3,             
                             features = 5,           
                             add_lm = TRUE,
                             color_by = "mitonuclear"
) + scale_color_manual(values = geno_colors) +
  ggtitle("Microbiome Features vs Factor 3") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    strip.text = element_text(size = 5.6), 
    axis.text = element_text(size = 8)
  )
# linearise transcriptome
panel_f <- plot_data_scatter(model,
                             view = "fly",         
                             factor = 3,             
                             features = 5,           
                             add_lm = TRUE,
                             color_by = "mitonuclear"
) + scale_color_manual(values = geno_colors) +
  ggtitle("Transcriptome Features vs Factor 3") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    strip.text = element_text(size = 8),  
    axis.text = element_text(size = 8)
  )

# Combine all panels for Factor 3
combined_plot <- grid.arrange(
  panel_a, panel_b,
  panel_c, panel_d,
  panel_e, panel_f,
  ncol = 2, nrow = 3,
  heights = c(1, 1, 1.2),
  top = textGrob("Factor 3 Analysis", gp = gpar(fontsize = 16, fontface = "bold"))
)

###################################################
## SUMMARY FOR FACTOR 1 ##
# Copy + paste job from F3 code 
###################################################

# Get top features for Factor 1
top_micro_f1 <- get_top_features(model, 1, "micro", 10)
top_fly_f1 <- get_top_features(model, 1, "fly", 10)

panel_a_f1 <- plot_factor(model, 
                          factors = 1,
                          color_by = "mitonuclear",
                          add_violin = TRUE,
                          dodge = TRUE) + 
  scale_fill_manual(values = geno_colors) +
  scale_color_manual(values = geno_colors) +
  ggtitle("Factor 1 Activity by Genotype") +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))

panel_b_f1 <- ggplot(top_micro_f1, aes(x = reorder(feature, abs_value), y = value)) +
  geom_col(aes(fill = value > 0), width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "blue")) +
  coord_flip() +
  labs(title = "Top Microbiome Features",
       x = "Features", 
       y = "Weight") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 8))

panel_c_f1 <- ggplot(top_fly_f1, aes(x = reorder(feature, abs_value), y = value)) +
  geom_col(aes(fill = value > 0), width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "blue")) +
  coord_flip() +
  labs(title = "Top Transcriptome Features",
       x = "Features", 
       y = "Weight") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 8))

# Correlations
cor_matrix_f1 <- create_cross_correlations(model, 1, 5)
cor_df_f1 <- expand.grid(Microbiome = rownames(cor_matrix_f1),
                         Transcriptome = colnames(cor_matrix_f1))
cor_df_f1$Correlation <- as.vector(cor_matrix_f1)

panel_d_f1 <- ggplot(cor_df_f1, aes(x = Microbiome, y = Transcriptome, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", 
                       midpoint = 0, limits = c(-1, 1)) +
  geom_text(aes(label = round(Correlation, 2)), size = 2) +
  labs(title = "Microbiome-Transcriptome Correlations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.2, "cm"),      # Make legend keys smaller
        legend.key.width = unit(0.1, "cm"),     # Make legend narrower
        legend.title = element_text(size = 8),   # Smaller legend title
        legend.text = element_text(size = 7),    # Smaller legend text
        legend.position = c(1, 0.21))

# FIXED: Remove invalid parameters
panel_e_f1 <- plot_data_scatter(model,
                                view = "micro",         
                                factor = 1,             
                                features = 5,
                                add_lm = TRUE,
                                color_by = "mitonuclear"
) + scale_color_manual(values = geno_colors) +
  ggtitle("Microbiome Features vs Factor 1") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    strip.text = element_text(size = 5.6), 
    axis.text = element_text(size = 8)
  )

panel_f_f1 <- plot_data_scatter(model,
                                view = "fly",         
                                factor = 1,             
                                features = 5,
                                add_lm = TRUE,
                                color_by = "mitonuclear"
) + scale_color_manual(values = geno_colors) +
  ggtitle("Transcriptome Features vs Factor 1") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    strip.text = element_text(size = 8),  
    axis.text = element_text(size = 8)
  )

# Combine all panels for Factor 1
combined_plot_f1 <- grid.arrange(
  panel_a_f1, panel_b_f1,
  panel_c_f1, panel_d_f1,
  panel_e_f1, panel_f_f1,
  ncol = 2, nrow = 3,
  heights = c(1, 1, 1.2),
  top = textGrob("Factor 1 Analysis", gp = gpar(fontsize = 16, fontface = "bold"))
)

###########################
## SUMMARY FOR FACTOR 2 ##
###########################

# Get top features for Factor 2
top_micro_f2 <- get_top_features(model, 2, "micro", 10)
top_fly_f2 <- get_top_features(model, 2, "fly", 10)

panel_a_f2 <- plot_factor(model, 
                          factors = 2,
                          color_by = "mitonuclear",
                          add_violin = TRUE,
                          dodge = TRUE) + 
  scale_fill_manual(values = geno_colors) +
  scale_color_manual(values = geno_colors) +
  ggtitle("Factor 2 Activity by Genotype") +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))

panel_b_f2 <- ggplot(top_micro_f2, aes(x = reorder(feature, abs_value), y = value)) +
  geom_col(aes(fill = value > 0), width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "blue")) +
  coord_flip() +
  labs(title = "Top Microbiome Features",
       x = "Features", 
       y = "Weight") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 8))

panel_c_f2 <- ggplot(top_fly_f2, aes(x = reorder(feature, abs_value), y = value)) +
  geom_col(aes(fill = value > 0), width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "blue")) +
  coord_flip() +
  labs(title = "Top Transcriptome Features",
       x = "Features", 
       y = "Weight") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 8))

# Correlations
cor_matrix_f2 <- create_cross_correlations(model, 2, 5)
cor_df_f2 <- expand.grid(Microbiome = rownames(cor_matrix_f2),
                         Transcriptome = colnames(cor_matrix_f2))
cor_df_f2$Correlation <- as.vector(cor_matrix_f2)

panel_d_f2 <- ggplot(cor_df_f2, aes(x = Microbiome, y = Transcriptome, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", 
                       midpoint = 0, limits = c(-1, 1)) +
  geom_text(aes(label = round(Correlation, 2)), size = 3) +
  labs(title = "Microbiome-Transcriptome Correlations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.2, "cm"),      # Make legend keys smaller
        legend.key.width = unit(0.1, "cm"),     # Make legend narrower
        legend.title = element_text(size = 8),   # Smaller legend title
        legend.text = element_text(size = 7),    # Smaller legend text
        legend.position = c(1, 0.21))

# Remove invalid parameters 
panel_e_f2 <- plot_data_scatter(model,
                                view = "micro",         
                                factor = 2,             
                                features = 5,
                                add_lm = TRUE,
                                color_by = "mitonuclear"
) + scale_color_manual(values = geno_colors) +
  ggtitle("Microbiome Features vs Factor 2") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    strip.text = element_text(size = 5.6), 
    axis.text = element_text(size = 8)
  )

panel_f_f2 <- plot_data_scatter(model,
                                view = "fly",         
                                factor = 2,             
                                features = 5,
                                add_lm = TRUE,
                                color_by = "mitonuclear"
) + scale_color_manual(values = geno_colors) +
  ggtitle("Transcriptome Features vs Factor 2") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    strip.text = element_text(size = 8),  
    axis.text = element_text(size = 8)
  )

# Combine all panels for Factor 2
combined_plot_f2 <- grid.arrange(
  panel_a_f2, panel_b_f2,
  panel_c_f2, panel_d_f2,
  panel_e_f2, panel_f_f2,
  ncol = 2, nrow = 3,
  heights = c(1, 1, 1.2),
  top = textGrob("Factor 2 Analysis", gp = gpar(fontsize = 16, fontface = "bold"))
)



#################################################################
# Isolate percentage of variance explained (linearisation R2) values
#################################################################
# function to extract R-squared values from linear relationships
extract_r_squared <- function(model, factor_num, view_name, n_features = 5) {
  
  # get the data
  data <- get_data(model, views = view_name, as.data.frame = TRUE)
  factor_values <- get_factors(model, factors = factor_num, as.data.frame = TRUE)
  
  # Convert feature columns to character to avoid missmatch issues
  data$feature <- as.character(data$feature)
  
  # Get top weighted features for this factor
  weights <- get_weights(model, 
                        views = view_name, 
                        factors = factor_num, 
                        as.data.frame = TRUE,
                        scale = TRUE)
  weights$feature <- as.character(weights$feature)  # convert to character
  weights$abs_value <- abs(weights$value)
  top_features <- weights[order(weights$abs_value, decreasing = TRUE), ][1:n_features, ]
  
  # Init results dataframe
  r_squared_results <- data.frame(
    factor = paste0("Factor", factor_num),
    view = view_name,
    feature = character(n_features),
    r_squared = numeric(n_features),
    p_value = numeric(n_features),
    stringsAsFactors = FALSE
  )
  
  # Calculate R squared for each top feature
  for(i in 1:n_features) {
    feature_name <- top_features$feature[i]
    
    # Get feature data
    feature_data <- data[data$feature == feature_name, ]
    
    # merge with factor values
    merged_data <- merge(feature_data, factor_values, by = "sample")
    # if we have merged data...
    if(nrow(merged_data) > 0) {
      # fit linear model
      lm_fit <- lm(value.x ~ value.y, data = merged_data)
      
      # Extract r squared and p-value
      r_squared_results$feature[i] <- feature_name
      r_squared_results$r_squared[i] <- summary(lm_fit)$r.squared
      r_squared_results$p_value[i] <- summary(lm_fit)$coefficients[2, 4]  # p-value for slope
    } else {
      # Handle case where no data is found
      r_squared_results$feature[i] <- feature_name
      r_squared_results$r_squared[i] <- NA
      r_squared_results$p_value[i] <- NA
    }
  }
  
  return(r_squared_results)
}

# Extract R squared values for factors of interest
r_squared_table <- data.frame()

# For Factors 1, 2, and 3, both microbiome and transcriptome
for(factor_num in 1:3) {
  # Microbiome
  micro_r2 <- extract_r_squared(model, factor_num, "micro", 5)
  r_squared_table <- rbind(r_squared_table, micro_r2)
  
  # transcriptome (fly gut)
  fly_r2 <- extract_r_squared(model, factor_num, "fly", 5)
  r_squared_table <- rbind(r_squared_table, fly_r2)
}

# Clean up the table
r_squared_table$feature <- gsub("_micro|_fly", "", r_squared_table$feature)
r_squared_table$r_squared <- round(r_squared_table$r_squared, 3)
r_squared_table$p_value <- round(r_squared_table$p_value, 4)

# Add our significance codes
r_squared_table$significance <- ifelse(r_squared_table$p_value < 0.001, "***",
                                      ifelse(r_squared_table$p_value < 0.01, "**",
                                             ifelse(r_squared_table$p_value < 0.05, "*", "ns")))

# Sort by factor and R-squared value
r_squared_table <- r_squared_table[order(r_squared_table$factor, -r_squared_table$r_squared), ]

# Remove any rows with na values
r_squared_table <- r_squared_table[!is.na(r_squared_table$r_squared), ]

# quick check...
print(r_squared_table)

# we want a more formatted version to match the other nice tables
library(knitr)
kable(r_squared_table, 
      caption = "R-squared values for linear relationships between factors and features",
      row.names = FALSE)

library(kableExtra)
mofa_r_squared_results <- data.frame(
  Factor = r_squared_table$factor,
  View = r_squared_table$view,
  Feature = r_squared_table$feature,
  `R-squared` = r_squared_table$r_squared,
  `P-value` = ifelse(r_squared_table$p_value < 0.001, 
                     format(r_squared_table$p_value, scientific = TRUE, digits = 2),
                     r_squared_table$p_value), 
  Significance = r_squared_table$significance, # tap into the pre-made sig codes
  check.names = FALSE
)

# Create the data frame with column names first
mofa_r_squared_simple <- data.frame(
  Factor = r_squared_table$factor,
  View = r_squared_table$view,
  Feature = r_squared_table$feature,
  R_squared = r_squared_table$r_squared,
  P_value = ifelse(r_squared_table$p_value < 0.001, 
                   format(r_squared_table$p_value, scientific = TRUE, digits = 2),
                   r_squared_table$p_value),
  Significance = r_squared_table$significance
)

# Create a simple kable table
kable(mofa_r_squared_simple, 
      caption = "Summary for linear relationships between factors and top weighted features",
      row.names = FALSE,
      col.names = c("Factor", "View", "Feature", "R²", "P-value", "Significance")) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), 
                full_width = FALSE) %>%
  footnote(general = "Significance codes: *** P < 0.001, ** P < 0.01, * P < 0.05, ns = not significant")



# Export as CSV - exporting as a png is a NIGHTMARE
write.csv(mofa_r_squared_simple, 
          "C:/Users/sandr/OneDrive/Desktop/mofa/mofa_table.csv", 
          row.names = FALSE)

### DONE ###