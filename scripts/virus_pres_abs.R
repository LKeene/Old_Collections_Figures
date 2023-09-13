library(tidyverse)
library(viridis)

other_viruses <- read_csv("metadata/OtherViruses_known.csv")

other_viruses <- other_viruses %>% 
  mutate(id = str_c(sample_location," ",year," ",rep)) %>% 
  mutate(virus_name = as.factor(fct_rev(taxon_id)),
         id = factor(id, levels = c("California 2011 1", "Ontario, CAN 2010 1", 
                                    "North Carolina 2006 1", "North Carolina 2006 2", 
                                    "California 2000 1", "California 2000 2",
                                    "Pennsylvania 1963 1", "Pennsylvania 1963 2", "Hawai'i 1954 1", 
                                    "Hawai'i 1953 1", "New Jersey 1942 1", 
                                    "Illinois 1930 1", "New York 1927 1",
                                    "Minnesota 1919 1", "Illinois 1915 1", 
                                    "Illinois 1915 2", "Illinois 1908 1")))

fig <- ggplot(other_viruses, aes(x = id, y = virus_name)) +
  geom_point(aes(color = percent_coverage, shape = factor(Drosophila_species)), size = 5) +
  scale_color_viridis() +
  scale_shape_manual(labels = c(substitute(paste(italic("D. melanogaster"))), 
                                substitute(paste(italic("D. simulans")))),
                     values = c(15, 19)) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Sample ID", y = "Virus Name", color = "% of Sequence Covered", shape = "Species") 

fig
ggsave("plots/virus_pres_abs.pdf", units = "in", width = 10, height = 8)
