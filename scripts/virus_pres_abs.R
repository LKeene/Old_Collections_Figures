library(tidyverse)
library(viridis)

other_viruses <- read_csv("metadata/OtherViruses_known.csv")

other_viruses <- other_viruses %>% 
  mutate(sample_id = factor(sample_id, levels = c("BBDIV1550-12", "PHDIP1072-11", 
                                    "Davidson2006_1", "Davidson2006_2", 
                                    "Frost_4908", "SB_2000_sim_1",
                                    "SB_2000_sim_2", "Frost03211963", "Frost_03071963", 
                                    "Kilauea1953_1_sim", "Kauai1953_1_mel", 
                                    "1004278", "1004283",
                                    "Albany1927", "StPaul1919_3", 
                                    "1004281", "1004279", "1004277")))

fig <- ggplot(other_viruses, aes(x = sample_id, y = virus_name)) +
  geom_point(aes(color = percent_covered, shape = factor(Drosophila_species)), size = 5) +
  scale_color_viridis() +
  scale_shape_manual(labels = c(substitute(paste(italic("D. melanogaster"))), 
                                substitute(paste(italic("D. simulans"))),
                                substitute("Unknown")),
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
  labs(x = "Sample ID", y = "Virus Name", color = "% of Sequence \nCovered", shape = "Species") 

fig
ggsave("plots/virus_pres_abs.pdf", units = "in", width = 10, height = 8)
ggsave("plots/virus_pres_abs.jpg", units = "in", width = 10, height = 8)
