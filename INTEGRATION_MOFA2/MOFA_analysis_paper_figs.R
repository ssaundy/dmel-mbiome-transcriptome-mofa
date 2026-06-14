rm(list=ls())
###
##
#
setwd("C:/Users/sandr/OneDrive - University of Glasgow/PROJECT/DATA/INTEGRATION")
#
##
###
library(ggplot2)
library(MOFA2)
library(pheatmap)
library(gridExtra)
library(grid)
library(car)
library(knitr)
library(kableExtra)
library(dplyr)
library(tidyr)
library(patchwork)

source("OMICS_mofa.R")

# Assign colours to genotypes
geno_colors <- c("AA"="orange", "AB"="red", "BA"="pink", "BB"="purple")

###############################################
## Load model and attach metadata
###############################################

model <- load_model("model.hdf5")
plot_data_overview(model)

samps <- samples_names(model)[[1]]

sample_metadata <- read.csv("samples.csv")
colnames(sample_metadata)[1] <- "samp"

# Create mitonuclear factor
sample_metadata$mitonuclear <- with(sample_metadata,
                                    factor(paste(mito, nuclear, sep="")))

# Subset and align to model sample order
sample_metadata <- droplevels(subset(sample_metadata, samp %in% samps))
sample_metadata <- sample_metadata[match(samps, sample_metadata$samp), ]
colnames(sample_metadata)[1] <- "sample"

samples_metadata(model) <- sample_metadata

# Check genotype distribution
table(sample_metadata$mitonuclear)

# Clean feature name suffixes (used downstream in scatter plots)
current_features <- features_names(model)
current_features$micro <- gsub("_micro", "", current_features$micro)
current_features$fly   <- gsub("_fly",   "", current_features$fly)
features_names(model) <- current_features


###############################################
## FUNCTIONS
###############################################

get_top_features <- function(model, factor_num, view_name, n_features = 20) {
  weights <- get_weights(model,
                         views   = view_name,
                         factors = factor_num,
                         as.data.frame = TRUE,
                         scale   = TRUE)
  weights$abs_value <- abs(weights$value)
  weights[order(weights$abs_value, decreasing = TRUE), ][1:n_features, ]
}

create_cross_correlations <- function(model, factor_num, n_features = 15) {
  micro_weights <- get_top_features(model, factor_num, "micro", n_features)
  fly_weights   <- get_top_features(model, factor_num, "fly",   n_features)
  
  micro_data <- get_data(model, views = "micro", as.data.frame = TRUE)
  fly_data   <- get_data(model, views = "fly",   as.data.frame = TRUE)
  
  micro_data$feature    <- as.character(micro_data$feature)
  fly_data$feature      <- as.character(fly_data$feature)
  micro_weights$feature <- as.character(micro_weights$feature)
  fly_weights$feature   <- as.character(fly_weights$feature)
  
  cor_matrix <- matrix(nrow = nrow(micro_weights), ncol = nrow(fly_weights))
  
  for (i in 1:nrow(micro_weights)) {
    micro_vals <- micro_data[micro_data$feature == micro_weights$feature[i], "value"]
    for (j in 1:nrow(fly_weights)) {
      fly_vals <- fly_data[fly_data$feature == fly_weights$feature[j], "value"]
      if (length(micro_vals) == length(fly_vals) && length(micro_vals) > 0) {
        cor_matrix[i, j] <- cor(micro_vals, fly_vals, use = "complete.obs")
      }
    }
  }
  
  rownames(cor_matrix) <- micro_weights$feature
  colnames(cor_matrix) <- fly_weights$feature
  cor_matrix
}

extract_r_squared <- function(model, factor_num, view_name, n_features = 5) {
  data          <- get_data(model, views = view_name, as.data.frame = TRUE)
  factor_values <- get_factors(model, factors = factor_num, as.data.frame = TRUE)
  data$feature  <- as.character(data$feature)
  
  weights <- get_weights(model,
                         views   = view_name,
                         factors = factor_num,
                         as.data.frame = TRUE,
                         scale   = TRUE)
  weights$feature   <- as.character(weights$feature)
  weights$abs_value <- abs(weights$value)
  top_features <- weights[order(weights$abs_value, decreasing = TRUE), ][1:n_features, ]
  
  results <- data.frame(
    factor      = paste0("Factor", factor_num),
    view        = view_name,
    feature     = character(n_features),
    r_squared   = numeric(n_features),
    p_value     = numeric(n_features),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n_features) {
    feature_data <- data[data$feature == top_features$feature[i], ]
    merged       <- merge(feature_data, factor_values, by = "sample")
    
    if (nrow(merged) > 0) {
      lm_fit <- lm(value.x ~ value.y, data = merged)
      results$feature[i]   <- top_features$feature[i]
      results$r_squared[i] <- summary(lm_fit)$r.squared
      results$p_value[i]   <- summary(lm_fit)$coefficients[2, 4]
    } else {
      results$feature[i]   <- top_features$feature[i]
      results$r_squared[i] <- NA
      results$p_value[i]   <- NA
    }
  }
  results
}

plot_factor_summary <- function(model, factor_num, geno_colors) {
  
  top_micro <- get_top_features(model, factor_num, "micro", 10)
  top_fly   <- get_top_features(model, factor_num, "fly",   10)
  
  panel_a <- plot_factor(model,
                         factors    = factor_num,
                         color_by   = "mitonuclear",
                         add_violin = TRUE,
                         dodge      = TRUE) +
    scale_fill_manual(values  = geno_colors) +
    scale_color_manual(values = geno_colors) +
    ggtitle(paste("Factor", factor_num, "Activity by Genotype")) +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title      = element_text(hjust = 0.5))
  
  panel_b <- ggplot(top_micro, aes(x = reorder(feature, abs_value), y = value)) +
    geom_col(aes(fill = value > 0), width = 0.7) +
    scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "blue")) +
    coord_flip() +
    labs(title = "Top Microbiome Features", x = "Features", y = "Weight") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title      = element_text(hjust = 0.5),
          axis.text.y     = element_text(size = 8))
  
  panel_c <- ggplot(top_fly, aes(x = reorder(feature, abs_value), y = value)) +
    geom_col(aes(fill = value > 0), width = 0.7) +
    scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "blue")) +
    coord_flip() +
    labs(title = "Top Transcriptome Features", x = "Features", y = "Weight") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title      = element_text(hjust = 0.5),
          axis.text.y     = element_text(size = 8))
  
  cor_matrix <- create_cross_correlations(model, factor_num, 5)
  cor_df     <- expand.grid(Microbiome    = rownames(cor_matrix),
                            Transcriptome = colnames(cor_matrix))
  cor_df$Correlation <- as.vector(cor_matrix)
  
  panel_d <- ggplot(cor_df, aes(x = Microbiome, y = Transcriptome, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", midpoint = 0, limits = c(-1, 1)) +
    geom_text(aes(label = round(Correlation, 2)), size = 3) +
    labs(title = "Microbiome-Transcriptome Correlations") +
    theme_minimal() +
    theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y      = element_text(size = 8),
          plot.title       = element_text(hjust = 0.5),
          legend.key.size  = unit(0.2, "cm"),
          legend.key.width = unit(0.1, "cm"),
          legend.title     = element_text(size = 8),
          legend.text      = element_text(size = 7),
          legend.position  = c(1, 0.21))
  
  panel_e <- plot_data_scatter(model,
                               view     = "micro",
                               factor   = factor_num,
                               features = 5,
                               add_lm   = TRUE,
                               color_by = "mitonuclear") +
    scale_color_manual(values = geno_colors) +
    ggtitle(paste("Microbiome Features vs Factor", factor_num)) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          strip.text = element_text(size = 5.6),
          axis.text  = element_text(size = 8))
  
  panel_f <- plot_data_scatter(model,
                               view     = "fly",
                               factor   = factor_num,
                               features = 5,
                               add_lm   = TRUE,
                               color_by = "mitonuclear") +
    scale_color_manual(values = geno_colors) +
    ggtitle(paste("Transcriptome Features vs Factor", factor_num)) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          strip.text = element_text(size = 8),
          axis.text  = element_text(size = 8))
  
  grid.arrange(
    panel_a, panel_b,
    panel_c, panel_d,
    panel_e, panel_f,
    ncol = 2, nrow = 3,
    heights = c(1, 1, 1.2),
    top = textGrob(paste("Factor", factor_num, "Analysis"),
                   gp = gpar(fontsize = 16, fontface = "bold"))
  )
}


###############################################
## FIGURE 3
###############################################

## VARIANCE EXPLAINED

plot_variance_explained(model, x = "view", y = "factor", factors = 1:10) +
  theme_classic() +
  labs(title = "Variance Explained") +
  theme(strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5))

var_data <- get_variance_explained(model)$r2_per_factor$group1[1:10, ]

p_var <- plot_variance_explained(model, x = "view", y = "factor", factors = 1:10) +
  theme_classic() +
  labs(title = "Variance Explained") +
  theme(strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5))

p_var + geom_text(
  data = expand.grid(factor = paste0("Factor", 1:10),
                     view   = c("fly", "micro")),
  aes(x = view, y = factor, label = round(as.vector(var_data), 1)),
  inherit.aes = FALSE, color = "grey", size = 7
)

var_total <- get_variance_explained(model)$r2_total$group1
print(var_total)

plot_variance_explained(model, x = "group", y = "factor", plot_total = TRUE) +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())


## FACTOR VALUES BY GENOTYPE

plot_factor(model,
            color_by   = "mitonuclear",
            factors    = 1:6,
            add_violin = TRUE,
            dodge      = TRUE) +
  scale_fill_manual(values  = geno_colors) +
  scale_color_manual(values = geno_colors)


## ANOVA Factors ~ mitonuclear genotype 

factor_values           <- get_factors(model, as.data.frame = TRUE)
factor_values_with_meta <- merge(factor_values, sample_metadata, by = "sample")
unique_factors          <- unique(factor_values_with_meta$factor)

factor_anova_table <- data.frame(factor      = unique_factors,
                                 f_statistic = NA,
                                 p_value     = NA)

for (i in seq_along(unique_factors)) {
  factor_data <- subset(factor_values_with_meta, factor == unique_factors[i])
  aov_summary <- summary(aov(value ~ mitonuclear, data = factor_data))[[1]]
  factor_anova_table$f_statistic[i] <- aov_summary[["F value"]][1]
  factor_anova_table$p_value[i]     <- aov_summary[["Pr(>F)"]][1]
}



factor_anova_table <- factor_anova_table[order(factor_anova_table$p_value), ]
factor_anova_table$significance <- ifelse(factor_anova_table$p_value < 0.001, "***",
                                          ifelse(factor_anova_table$p_value < 0.01,  "**",
                                                 ifelse(factor_anova_table$p_value < 0.05,  "*", "ns")))
factor_anova_table$factor       <- gsub("Factor", "", factor_anova_table$factor)
factor_anova_table$f_statistic  <- round(factor_anova_table$f_statistic, 2)
factor_anova_table$p_value_display <- ifelse(factor_anova_table$p_value < 0.001,
                                             format(factor_anova_table$p_value,
                                                    scientific = TRUE, digits = 2),
                                             round(factor_anova_table$p_value, 3))

final_anova_table <- data.frame(
  Factor        = factor_anova_table$factor,
  `F-statistic` = factor_anova_table$f_statistic,
  `P-value`     = factor_anova_table$p_value_display,
  Significance  = factor_anova_table$significance
)

print(final_anova_table)

mofa_results_simple <- data.frame(
  Factor       = final_anova_table$Factor,
  F_stat       = final_anova_table$F.statistic,
  P_value      = final_anova_table$P.value,
  Significance = final_anova_table$Significance
)
attr(mofa_results_simple, "title") <- "MOFA Factor Associations with Mitonuclear Genotype"
create_formatted_kable_mofa(mofa_results_simple[1:10, ],
                            "- ANOVA testing genotype effects on factor values")


## ANOVA Factors ~ microbiota 

micro_data_wide <- get_data(model, views = "micro", as.data.frame = TRUE)
micro_data_wide$feature <- as.character(micro_data_wide$feature)

N_ASVS <- 10 #
agg     <- aggregate(value ~ feature, data = micro_data_wide, FUN = mean, na.rm = TRUE)
agg     <- agg[order(agg$value, decreasing = TRUE), ]
top_asvs <- agg$feature[1:N_ASVS]

micro_wide <- micro_data_wide %>%
  dplyr::filter(feature %in% top_asvs) %>%
  pivot_wider(names_from = feature, values_from = value)

micro_anova_table <- data.frame(
  factor      = character(),
  f_statistic = numeric(),
  p_value     = numeric(),
  r_squared   = numeric()
)

for (factor_num in 1:10) {
  factor_data <- subset(factor_values_with_meta,
                        factor == paste0("Factor", factor_num))
  merged      <- merge(factor_data, micro_wide, by = "sample")
  formula_str <- paste("value ~", paste(paste0("`", top_asvs, "`"), collapse = " + "))
  lm_fit      <- lm(as.formula(formula_str), data = merged)
  fstat       <- summary(lm_fit)$fstatistic
  p_value     <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
  
  micro_anova_table <- rbind(micro_anova_table, data.frame(
    factor      = paste0("Factor", factor_num),
    f_statistic = round(fstat[1], 2),
    p_value     = p_value,
    r_squared   = round(summary(lm_fit)$r.squared, 3)
  ))
}

micro_anova_table$significance <- ifelse(micro_anova_table$p_value < 0.001, "***",
                                         ifelse(micro_anova_table$p_value < 0.01,  "**",
                                                ifelse(micro_anova_table$p_value < 0.05,  "*", "ns")))
micro_anova_table$p_value <- ifelse(micro_anova_table$p_value < 0.001,
                                    format(micro_anova_table$p_value, scientific = TRUE, digits = 2),
                                    round(micro_anova_table$p_value, 3))
print(micro_anova_table)


## Microbiome-Transcriptome correlation matrices 

for (factor_num in 1:10) {
  cor_matrix <- create_cross_correlations(model, factor_num, 5)
  pheatmap(cor_matrix,
           main            = paste("Factor", factor_num,
                                   ": Microbiome-Transcriptome Correlations"),
           display_numbers = TRUE,
           number_format   = "%.2f",
           fontsize_number = 8)
}


###############################################
## FIGURE 4
###############################################

for (factor_num in c(1, 2, 3, 4)) {
  plot_factor_summary(model, factor_num, geno_colors)
}

for (factor_num in 1:3) {
  print(
    plot_weights(model, view = "micro", factor = factor_num,
                 nfeatures = 15, scale = TRUE, abs = FALSE) +
      ggtitle(paste("Factor", factor_num, "- Microbiome Features")) +
      theme(plot.margin        = margin(10, 20, 10, 20),
            axis.text.y        = element_text(size = 9, margin = margin(r = 5)),
            panel.grid.major.y = element_line(size = 0.3, color = "grey90"))
  )
  print(
    plot_weights(model, view = "fly", factor = factor_num,
                 nfeatures = 10, scale = TRUE, abs = FALSE) +
      ggtitle(paste("Factor", factor_num, "- Transcriptome Features")) +
      theme(plot.margin        = margin(10, 20, 10, 20),
            axis.text.y        = element_text(size = 9, margin = margin(r = 5)),
            panel.grid.major.y = element_line(size = 0.3, color = "grey90"))
  )
}

r_squared_table <- data.frame()

for (factor_num in c(2, 4)) {
  r_squared_table <- rbind(r_squared_table,
                           extract_r_squared(model, factor_num, "micro", 5),
                           extract_r_squared(model, factor_num, "fly",   5))
}

r_squared_table$feature      <- gsub("_micro|_fly", "", r_squared_table$feature)
r_squared_table$r_squared    <- round(r_squared_table$r_squared, 3)
r_squared_table$p_value      <- round(r_squared_table$p_value,   4)
r_squared_table$significance <- ifelse(r_squared_table$p_value < 0.001, "***",
                                       ifelse(r_squared_table$p_value < 0.01,  "**",
                                              ifelse(r_squared_table$p_value < 0.05,  "*", "ns")))
r_squared_table <- r_squared_table[order(r_squared_table$factor,
                                         -r_squared_table$r_squared), ]
r_squared_table <- r_squared_table[!is.na(r_squared_table$r_squared), ]

print(r_squared_table)

mofa_r_squared_simple <- data.frame(
  Factor       = r_squared_table$factor,
  View         = r_squared_table$view,
  Feature      = r_squared_table$feature,
  R_squared    = r_squared_table$r_squared,
  P_value      = ifelse(r_squared_table$p_value < 0.001, "< 0.001",
                        r_squared_table$p_value),
  Significance = r_squared_table$significance
)

kable(mofa_r_squared_simple,
      caption   = "Linear relationships for genotype-associated factors (2 and 4)",
      row.names = FALSE,
      col.names = c("Factor", "View", "Feature", "R²", "P-value", "Significance")) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) %>%
  footnote(general = "Significance codes: *** P<0.001, ** P<0.01, * P<0.05, ns = not significant")

write.csv(mofa_r_squared_simple,
          "C:/Users/sandr/OneDrive/Desktop/mofa/mofa_r_squared_table.csv",
          row.names = FALSE)

### ANALYSIS DONE ###

###############################################
## PAPER FIGURES
###############################################

library(patchwork)
### Enter OUTPUT DIRECTORY ###
 ##                        ##
  ##                     ##
   # Desktop in my case #
      ################
OUT_DIR <- "C:/Users/sandr/OneDrive/Desktop/mofa/figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

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

sig_factors <- 1:10

save_plot <- function(plot_obj, filename, width = 8, height = 6) {
  ggsave(
    filename    = file.path(OUT_DIR, filename),
    plot        = plot_obj,
    width       = width,
    height      = height,
    units       = "cm",
    device      = "pdf",
    useDingbats = FALSE
  )
  message("Saved: ", filename)
}


###############################################
## FIGURE 3 
###############################################

### SAVE ANOVA

# Save microbiota ANOVA table
micro_anova_display <- data.frame(
  Factor       = micro_anova_table$factor,
  F_statistic  = micro_anova_table$f_statistic,
  R_squared    = micro_anova_table$r_squared,
  P_value      = micro_anova_table$p_value,
  Significance = micro_anova_table$significance
)

micro_anova_grob <- tableGrob(
  micro_anova_display,
  rows  = NULL,
  theme = ttheme_minimal(
    base_size   = 9,
    base_family = "sans",
    core    = list(fg_params = list(hjust = 0.5, x = 0.5)),
    colhead = list(fg_params = list(hjust = 0.5, x = 0.5, fontface = "bold"))
  )
)

micro_anova_titled <- arrangeGrob(
  textGrob("C2  ANOVA: Factor values ~ microbiota composition",
           x = 0, hjust = 0,
           gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "sans")),
  micro_anova_grob,
  nrow = 2, heights = c(0.08, 1)
)

ggsave(file.path(OUT_DIR, "Fig3C2_micro_anova_table.pdf"),
       plot   = micro_anova_titled,
       width  = 12, height = 10,
       units  = "cm", device = "pdf", useDingbats = FALSE)
message("Saved: Fig3C2_micro_anova_table.pdf")

## Variance explained

var_data <- get_variance_explained(model)$r2_per_factor$group1[1:10, ]

var_df <- as.data.frame(var_data) %>%
  mutate(factor = paste0("Factor", 1:10)) %>%
  pivot_longer(-factor, names_to = "view", values_to = "r2") %>%
  mutate(
    factor = factor(factor, levels = paste0("Factor", 10:1)),
    view   = dplyr::recode(view, "fly" = "Transcriptome", "micro" = "Microbiome")
  )

fig3a <- ggplot(var_df, aes(x = view, y = factor, fill = r2)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(r2, 1),
                colour = r2 > 5),   # threshold: adjust to your scale
            size = 2.8) +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
  scale_fill_gradientn(
    colors = c("#F7FBFF", "#6BAED6", "#08306B"),
    name   = "R² (%)",
    limits = c(0, NA)
  ) +
  scale_x_discrete(position = "top") +
  labs(title = "A", x = NULL, y = "Factor") +
  theme_msb() +
  theme(axis.line = element_blank(), axis.ticks = element_blank())

save_plot(fig3a, "Fig3A_variance_explained.pdf", width = 8, height = 9)


##  Factor violins — significant factors only (main figure)

fig3b <- plot_factor(model,
                     factors    = c(2, 4, 6),
                     color_by   = "mitonuclear",
                     add_violin = TRUE,
                     dodge      = TRUE) +
  scale_fill_manual(values  = geno_colors_msb, name = "Genotype") +
  scale_color_manual(values = geno_colors_msb, name = "Genotype") +
  labs(title = "B", x = "Genotype", y = "Factor value") +
  theme_msb() +
  theme(legend.position = "right",
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank())

save_plot(fig3b, "Fig3B_factor_violins_sig.pdf", width = 12, height = 6)


##  Factor violins — all 10 factors (supplemental)

fig3b_all <- plot_factor(model,
                         factors    = 1:10,
                         color_by   = "mitonuclear",
                         add_violin = TRUE,
                         dodge      = TRUE) +
  scale_fill_manual(values  = geno_colors_msb, name = "Genotype") +
  scale_color_manual(values = geno_colors_msb, name = "Genotype") +
  labs(title = "B (all factors)", x = "Genotype", y = "Factor value") +
  theme_msb() +
  theme(legend.position = "right",
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank())

save_plot(fig3b_all, "Fig3B_factor_violins_all10.pdf", width = 30, height = 6)


## ANOVA table 

anova_display <- final_anova_table %>%
  mutate(Factor = paste0("Factor ", Factor)) %>%
  head(10)

anova_grob <- tableGrob(
  anova_display,
  rows  = NULL,
  theme = ttheme_minimal(
    base_size   = 9,
    base_family = "sans",
    core    = list(fg_params = list(hjust = 0.5, x = 0.5)),
    colhead = list(fg_params = list(hjust = 0.5, x = 0.5, fontface = "bold"))
  )
)

anova_titled <- arrangeGrob(
  textGrob("C  ANOVA: Factor values ~ mitonuclear genotype",
           x = 0, hjust = 0,
           gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "sans")),
  anova_grob,
  nrow = 2, heights = c(0.08, 1)
)

ggsave(file.path(OUT_DIR, "Fig3C_anova_table.pdf"),
       plot = anova_titled, width = 10, height = 8,
       units = "cm", device = "pdf", useDingbats = FALSE)
message("Saved: Fig3C_anova_table.pdf")


## Correlation matrices 

for (factor_num in sig_factors) {
  cor_matrix <- create_cross_correlations(model, factor_num, 5)
  rownames(cor_matrix) <- gsub("_micro", "", rownames(cor_matrix))
  colnames(cor_matrix) <- gsub("_fly",   "", colnames(cor_matrix))
  
  p <- pheatmap(
    cor_matrix,
    main            = paste0("D  Factor ", factor_num,
                             ": Microbiome-Transcriptome Correlations"),
    display_numbers = TRUE,
    number_format   = "%.2f",
    fontsize        = 9,
    fontsize_number = 8,
    color           = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks          = seq(-1, 1, length.out = 101),
    border_color    = "white",
    treeheight_row  = 0,
    treeheight_col  = 0,
    silent          = TRUE
  )
  
  pdf(file.path(OUT_DIR, paste0("Fig3D_correlation_F", factor_num, ".pdf")),
      width = 10 / 2.54, height = 7 / 2.54,
      family = "sans", useDingbats = FALSE)
  grid::grid.newpage()
  grid::grid.draw(p$gtable)
  dev.off()
  message("Saved: Fig3D_correlation_F", factor_num, ".pdf")
}


###############################################
## FIGURE 4
###############################################

plot_loadings_clean <- function(model, factor_num, view_name,
                                n_features = 10, panel_label = "") {
  top_feats <- get_top_features(model, factor_num, view_name, n_features)
  top_feats$feature <- gsub("_micro|_fly", "", top_feats$feature)
  view_label <- ifelse(view_name == "micro", "Microbiome", "Transcriptome")
  
  ggplot(top_feats, aes(x = reorder(feature, value), y = value,
                        fill = value > 0)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC")) +
    coord_flip() +
    labs(title = paste0(panel_label, "  Factor ", factor_num,
                        " — ", view_label, " loadings"),
         x = NULL, y = "Weight (scaled)") +
    theme_msb() +
    theme(legend.position = "none")
}

plot_scatter_clean <- function(model, factor_num, view_name,
                               panel_label = "") {
  view_label <- ifelse(view_name == "micro", "Microbiome", "Transcriptome")
  
  # Get top features and factor values
  top_feats         <- get_top_features(model, factor_num, view_name, 5)
  top_feats$feature <- gsub("_micro|_fly", "", top_feats$feature)
  feat_names        <- top_feats$feature
  
  data_long         <- get_data(model, views = view_name, as.data.frame = TRUE)
  data_long$feature <- gsub("_micro|_fly", "", data_long$feature)
  factor_vals       <- get_factors(model, factors = factor_num, as.data.frame = TRUE)
  
  plot_df <- data_long %>%
    dplyr::filter(feature %in% feat_names) %>%
    merge(factor_vals, by = "sample") %>%
    merge(sample_metadata[, c("sample", "mitonuclear")], by = "sample") %>%
    mutate(feature = factor(feature, levels = feat_names))
  
  # Compute perfeature R and p, anchored at data minimum
  stats_df <- plot_df %>%
    group_by(feature) %>%
    summarise(
      R    = cor(value.y, value.x, use = "complete.obs"),
      pval = cor.test(value.y, value.x)$p.value,
      xpos = min(value.y, na.rm = TRUE),
      ypos = min(value.x, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      label = paste0("R = ", round(R, 2), "\np = ",
                     ifelse(pval < 0.001, "< 0.001", round(pval, 3)))
    )
  
  ggplot(plot_df, aes(x = value.y, y = value.x, colour = mitonuclear)) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, colour = "grey40",
                linewidth = 0.5, fill = "grey80") +
    geom_text(data = stats_df,
              aes(x = xpos, y = ypos, label = label),
              inherit.aes = FALSE,
              hjust = 0, vjust = 0, size = 2.5, colour = "black") +
    scale_colour_manual(values = geno_colors_msb, name = "Genotype") +
    facet_wrap(~ feature, nrow = 1) +
    labs(title = paste0(panel_label, "  Factor ", factor_num,
                        " — ", view_label, " abundance vs factor value"),
         x     = paste("Factor", factor_num, "value"),
         y     = paste(view_label, "abundance")) +
    theme_msb() +
    theme(strip.text = element_text(size = 7),
          axis.text  = element_text(size = 7))
}

# Save each panel individually per factor
for (factor_num in sig_factors) {
  
  p_micro_load <- plot_loadings_clean(model, factor_num, "micro", 10, "A")
  p_fly_load   <- plot_loadings_clean(model, factor_num, "fly",   10, "B")
  p_micro_scat <- plot_scatter_clean(model, factor_num, "micro", "C")
  p_fly_scat   <- plot_scatter_clean(model, factor_num, "fly",   "D")
  
  save_plot(p_micro_load,
            paste0("Fig4_Factor", factor_num, "_A_micro_loadings.pdf"),
            width = 10, height = 10)
  
  save_plot(p_fly_load,
            paste0("Fig4_Factor", factor_num, "_B_fly_loadings.pdf"),
            width = 10, height = 10)
  
  save_plot(p_micro_scat,
            paste0("Fig4_Factor", factor_num, "_C_micro_scatter.pdf"),
            width = 30, height = 8)
  
  save_plot(p_fly_scat,
            paste0("Fig4_Factor", factor_num, "_D_fly_scatter.pdf"),
            width = 30, height = 8)
}


###############################################
## GO/KEGG ENRICHMENT
###############################################

library(clusterProfiler)
library(org.Dm.eg.db)

N_GENES <- 50 # enough for strength but not too dilute

run_enrichment <- function(model, factor_num, n_genes = 50) {
  top_genes <- get_top_features(model, factor_num, "fly", n_genes)
  gene_list  <- gsub("_fly", "", as.character(top_genes$feature))
  
  entrez_ids <- bitr(gene_list,
                     fromType = "SYMBOL",
                     toType   = "ENTREZID",
                     OrgDb    = org.Dm.eg.db)

  
  if (nrow(entrez_ids) == 0) {
    return(NULL)
  }
  
  go_bp <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Dm.eg.db,
                    ont = "BP", pAdjustMethod = "BH",
                    pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  
  go_mf <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Dm.eg.db,
                    ont = "MF", pAdjustMethod = "BH",
                    pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  
  go_cc <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Dm.eg.db,
                    ont = "CC", pAdjustMethod = "BH",
                    pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  
  kegg  <- enrichKEGG(gene = entrez_ids$ENTREZID, organism = "dme",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  
  list(go_bp = go_bp, go_mf = go_mf, go_cc = go_cc, kegg = kegg, genes = gene_list)
}

enrichment_results <- list()

for (factor_num in sig_factors) {
  enrichment_results[[paste0("Factor", factor_num)]] <- run_enrichment(model, factor_num, N_GENES)
}

plot_enrichment <- function(enrich_obj, title, n_terms = 15) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
    return(NULL)
  }
  dotplot(enrich_obj, showCategory = n_terms) +
    labs(title = title) +
    theme_msb(base_size = 9) +
    theme(axis.text.y     = element_text(size = 7),
          plot.title      = element_text(size = 9, face = "bold"),
          legend.position = "right")
}

for (factor_num in sig_factors) {
  res <- enrichment_results[[paste0("Factor", factor_num)]]
  if (is.null(res)) next
  
  p_go_bp <- plot_enrichment(res$go_bp, paste0("Factor ", factor_num, " — GO Biological Process"))
  p_go_mf <- plot_enrichment(res$go_mf, paste0("Factor ", factor_num, " — GO Molecular Function"))
  p_go_cc <- plot_enrichment(res$go_cc, paste0("Factor ", factor_num, " — GO Cellular Component"))
  p_kegg  <- plot_enrichment(res$kegg,  paste0("Factor ", factor_num, " — KEGG Pathways"))
  
  if (!is.null(p_go_bp)) save_plot(p_go_bp, paste0("Fig4_Factor", factor_num, "_GO_BP.pdf"),  width = 14, height = 10)
  if (!is.null(p_go_mf)) save_plot(p_go_mf, paste0("Fig4_Factor", factor_num, "_GO_MF.pdf"),  width = 14, height = 10)
  if (!is.null(p_go_cc)) save_plot(p_go_cc, paste0("Fig4_Factor", factor_num, "_GO_CC.pdf"),  width = 14, height = 10)
  if (!is.null(p_kegg))  save_plot(p_kegg,  paste0("Fig4_Factor", factor_num, "_KEGG.pdf"),   width = 14, height = 10)
}

# Export enrichment CSV tables
for (factor_num in sig_factors) {
  res <- enrichment_results[[paste0("Factor", factor_num)]]
  if (is.null(res)) next
  
  if (!is.null(res$go_bp) && nrow(as.data.frame(res$go_bp)) > 0) {
    go_bp_df <- as.data.frame(res$go_bp) %>%
      dplyr::select(Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID) %>%
      arrange(p.adjust)
    write.csv(go_bp_df,
              file.path(OUT_DIR, paste0("Factor", factor_num, "_GO_BP_terms.csv")),
              row.names = FALSE)
  }
  
  if (!is.null(res$go_mf) && nrow(as.data.frame(res$go_mf)) > 0) {
    go_mf_df <- as.data.frame(res$go_mf) %>%
      dplyr::select(Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID) %>%
      arrange(p.adjust)
    write.csv(go_mf_df,
              file.path(OUT_DIR, paste0("Factor", factor_num, "_GO_MF_terms.csv")),
              row.names = FALSE)
  }
  
  if (!is.null(res$go_cc) && nrow(as.data.frame(res$go_cc)) > 0) {
    go_cc_df <- as.data.frame(res$go_cc) %>%
      dplyr::select(Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID) %>%
      arrange(p.adjust)
    write.csv(go_cc_df,
              file.path(OUT_DIR, paste0("Factor", factor_num, "_GO_CC_terms.csv")),
              row.names = FALSE)
  }
  
  if (!is.null(res$kegg) && nrow(as.data.frame(res$kegg)) > 0) {
    kegg_df <- as.data.frame(res$kegg) %>%
      dplyr::select(Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID) %>%
      arrange(p.adjust)
    write.csv(kegg_df,
              file.path(OUT_DIR, paste0("Factor", factor_num, "_KEGG_terms.csv")),
              row.names = FALSE)
  }
}

