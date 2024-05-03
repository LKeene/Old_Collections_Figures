# define colors to be used in plots
# import (source) this file into other R 
# scripts to keep color use consistent
fresh_frozen_color       = "navyblue"
experimental_dried_color = "firebrick"
old_collection_color     = "orchid4"

theme_this_paper <- function(base_size = 12) {
    theme_minimal(base_size = base_size) +
    theme(panel.border = element_rect(linetype = "solid", fill = NA),
          strip.background = element_rect(colour = "black", fill = "grey95"),
          # strip.text = element_text(face = "bold"),
          # axis.text = element_text(face = "bold"),
          legend.position = "none") 
}

rotate_x_axis_labels <- function() {
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
}
