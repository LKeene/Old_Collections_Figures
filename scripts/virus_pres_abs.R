# This script visualizes identified virus sequences in museum collections as a 
# heatmap.

library(tidyverse)
library(viridis)
library(svglite)

other_viruses <- read_csv("metadata/OtherViruses_known.csv")

other_viruses <- other_viruses %>% 
  mutate(sample_id = factor(sample_id, levels = c("California_2011_2", "Canada_2010_1", 
                                    "NorthCarolina_2006_1", "NorthCarolina_2006_2", 
                                    "Pennsylvania_2003_1", "California_2000_1",
                                    "California_2000_2", "Pennsylvania_1963_1", "Pennsylvania_1963_2", 
                                    "Hawaii_1953_1", "Hawaii_1953_2", 
                                    "NewJersey_1942", "Illinois_1930",
                                    "NewYork_1927", "Minnesota_1919_3", 
                                    "Illinois_1915_1", "Illinois_1915_3", "Illinois_1908", 
                                    "Positive_Control_1", "Positive_Control_2"))) %>% 
  mutate(virus_name = factor(virus_name, levels = c("Puslinch virus", "Drosophila-associated sobemo-like virus", 
                                                    "Tobacco mosaic virus", "Thika virus", "Nora virus", 
                                                    "La Jolla virus", "Drosophila C virus", "Drosophila A virus",
                                                    "Dansoman virus", "Craigies Hill virus", 
                                                    "Bloomfield virus", "Drosophila melanogaster sigmavirus", 
                                                    "Chaq-like virus", "Vera virus", 
                                                    "Chaq virus", "Galbut virus")))

fig <- ggplot(other_viruses, aes(x = sample_id, y = virus_name)) +
  geom_point(aes(color = percent_covered, shape = factor(Drosophila_species)), size = 6) +
  scale_color_viridis(direction = -1, option = "A") +
  scale_shape_manual(labels = c(substitute(paste(italic("D. melanogaster"))), 
                                substitute(paste(italic("D. simulans"))),
                                substitute(paste(italic("Unknown Drosophilidae")))),
                     values = c(15, 17, 19)) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
#        axis.text = element_text(face = "bold"),
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "", y = "", color = "% of Sequence \nCovered", shape = "Species") 

fig
ggsave("plots/Figure3/virus_pres_abs.pdf", units = "in", width = 10, height = 8)
ggsave("plots/Figure3/virus_pres_abs.jpg", units = "in", width = 12, height = 8)
ggsave("plots/Figure3/virus_pres_abs.svg", units = "in", width = 12, height = 8)
