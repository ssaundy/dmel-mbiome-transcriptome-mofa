### THIS IS A SECTION FOR THEMES###

# CUSTOM THEME

my_theme = theme(
  
  plot.title = element_text(size=30),
  
  axis.text.x = element_text(size=16),
  
  axis.text.y = element_text(size=16),
  
  axis.title.x = element_text(size=25),
  
  axis.title.y = element_text(size=25)
  
)

# NATURE

nature_omics_theme = theme(
  # Title formatting
  plot.title = element_text(size = 18, face = "bold", margin = margin(b = 10)),
  
  # Axis formatting
  axis.title = element_text(size = 12, face = "plain"),
  axis.text = element_text(size = 10, colour = "black"),
  axis.line = element_line(colour = "black", size = 0.5),
  
  # Panel formatting
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_line(colour = "grey90", size = 0.2),
  panel.grid.minor = element_blank(),
  
  # Legend formatting
  legend.title = element_text(size = 10, face = "plain"),
  legend.text = element_text(size = 9),
  legend.key = element_rect(fill = "white"),
  legend.background = element_rect(fill = "white", colour = NA),
  legend.position = "right",
  
  # Margins and spacing
  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
  
  # Remove unnecessary boxes
  panel.border = element_blank(),
  
  # Facet formatting (if used)
  strip.background = element_rect(fill = "grey95"),
  strip.text = element_text(size = 10)
)


###########################################################################
###### Generalized table formatting functions #####
###########################################################################

# function to format p-values
format_p_value <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  else if (p < 0.01) return("< 0.01")
  else if (p < 0.05) return("< 0.05")
  else return(sprintf("%.3f", p))
}

# function to add significance indicator stars
add_significance <- function(p) {
  if (is.na(p)) return("NA")
  ifelse(p < 0.001, "***", 
         ifelse(p < 0.01, "**", 
                ifelse(p < 0.05, "*", "ns")))
}

# PERMANOVA table creator
create_permanova_table <- function(permanova_result, factor_names = NULL, table_title = "PERMANOVA Results") {
  
  # Extract results - handle both single and multiple factor cases
  if (is.null(factor_names)) {
    # Auto-detect factor names from the permanova result
    factor_names <- rownames(permanova_result)[!rownames(permanova_result) %in% c("Residual", "Total")]
  }
  
  n_factors <- length(factor_names)
  
  # Extract values for each factor
  r2_values <- round(permanova_result$R2[1:n_factors], 3)
  f_values <- round(permanova_result$F[1:n_factors], 2)
  p_values <- permanova_result$`Pr(>F)`[1:n_factors]
  
  # Create the table
  results_table <- data.frame(
    Factor = factor_names,
    R2 = r2_values,
    F_statistic = f_values,
    P_value = sapply(p_values, format_p_value),
    Significance = sapply(p_values, add_significance),
    stringsAsFactors = FALSE
  )
  
  # Add table title as attribute
  attr(results_table, "title") <- table_title
  
  return(results_table)
}

# Wrapper functions for my project specific analyses
create_population_table_general <- function(perm_result) {
  return(create_permanova_table(perm_result, 
                                factor_names = "Population",
                                table_title = "Population Effects"))
}

create_mitonuclear_table_general <- function(perm_additive, perm_interaction) {
  # Combine additive and interaction results
  combined_table <- rbind(
    create_permanova_table(perm_additive, 
                           factor_names = c("Nuclear", "Mitochondrial"),
                           table_title = "Mitonuclear Effects")[1:2, ],
    create_permanova_table(perm_interaction, 
                           factor_names = c("Nuclear", "Mitochondrial", "Nuclear ? Mitochondrial"),
                           table_title = "Mitonuclear Effects")[3, ]
  )
  
  # Reorder to match sensible order
  combined_table <- combined_table[c(2, 1, 3), ]  # Mito, Nuclear, Interaction
  combined_table$Factor <- c("Mitochondrial", "Nuclear", "Nuclear ? Mitochondrial")
  
  return(combined_table)
}

# Generic function to create formatted kable tables
create_formatted_kable <- function(results_table, caption_suffix = "") {
  
  table_title <- attr(results_table, "title")
  if (is.null(table_title)) table_title <- "PERMANOVA Results"
  
  full_caption <- paste0(table_title, " ", caption_suffix)
  
  kable(results_table, 
        caption = full_caption,
        col.names = c("Factor", "R?", "F", "P-value", "Sig."),
        align = c('l', 'c', 'c', 'c', 'c'),
        row.names = FALSE) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "bordered"), 
                  full_width = FALSE,
                  position = "center") %>%
    footnote(general = "Significance codes: *** P < 0.001, ** P < 0.01, * P < 0.05, ns = not significant")
}
##########################################
### reformatted mofa formatter for table
##########################################
create_formatted_kable_mofa <- function(results_table, caption_suffix = "") {
  
  table_title <- attr(results_table, "title")
  if (is.null(table_title)) table_title <- "MOFA Results"
  
  full_caption <- paste0(table_title, " ", caption_suffix)
  
  kable(results_table, 
        caption = full_caption,
        col.names = c("Factor", "F-statistic", "P-value", "Sig."),
        align = c('l', 'c', 'c', 'c'),
        row.names = FALSE) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "bordered"), 
                  full_width = FALSE,
                  position = "center") %>%
    footnote(general = "Significance codes: *** P < 0.001, ** P < 0.01, * P < 0.05, ns = not significant")
}

# Then create a simpler table without R?
#mofa_results_simple <- data.frame(
 # Factor = factor_table$factor,
  #F_stat = factor_table$f_statistic,
  #P_value = factor_table$p_value_display,
  #Significance = factor_table$significance
#)

#attr(mofa_results_simple, "title") <- "MOFA Factor Associations with Mitonuclear Genotype"

# Use the modified function
#create_formatted_kable_mofa(mofa_results_simple, "- ANOVA testing genotype effects on factor values")
###