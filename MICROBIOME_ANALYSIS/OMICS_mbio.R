### THIS IS A SECTION FOR THEMES###
library(ggplot2)

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
