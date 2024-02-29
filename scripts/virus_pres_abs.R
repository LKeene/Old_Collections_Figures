library(tidyverse)
library(viridis)
library(svglite)

other_viruses <- read_csv("metadata/OtherViruses_known.csv")

other_viruses <- other_viruses %>% 
  mutate(sample_id = factor(sample_id, levels = c("California_2011_rep2", "Canada_2010_rep1", 
                                    "NorthCarolina_2006_rep1", "NorthCarolina_2006_rep2", 
                                    "Pennsylvania_2003_rep1", "California_2000_rep1",
                                    "California_2000_rep2", "Pennsylvania_1963_rep1", "Pennsylvania_1963_rep2", 
                                    "Hawaii_1953_rep1", "Hawaii_1953_rep2", 
                                    "NewJersey_1942", "Illinois_1930",
                                    "NewYork_1927", "Minnesota_1919_rep3", 
                                    "Illinois_1915_rep1", "Illinois_1915_rep3", "Illinois_1908"))) %>% 
  mutate(virus_name = factor(virus_name, levels = c("Puslinch virus", "Drosophila-associated sobemo-like virus", 
                                                    "Tobacco mosaic virus", "Nora virus", 
                                                    "La Jolla virus", "Drosophila C virus", 
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
ggsave("plots/virus_pres_abs.pdf", units = "in", width = 10, height = 8)
ggsave("plots/virus_pres_abs.jpg", units = "in", width = 12, height = 8)
ggsave("plots/virus_pres_abs.svg", units = "in", width = 12, height = 8)
