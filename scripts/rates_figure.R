# This script generates figures to visualize the different estimated evolutionary
# rates for the viruses described in this study.

library(tidyverse)
library(RColorBrewer)
library(svglite)

rates <- read_csv("metadata/rates_data.csv")

rates2 <- rates %>% 
  mutate(virus = factor(virus, levels = c("Galbut virus RNA 1", "Galbut virus RNA 2", 
                                          "Galbut virus RNA 3", "Chaq virus", 
                                          "Vera virus RNA 1", "Vera virus RNA 2", 
                                          "Chaq-like virus", 
                                          "Drosophila melanogaster sigmavirus", 
                                          "Craigies Hill virus RNA 1", 
                                          "Craigies Hill virus RNA 2", 
                                          "Dansoman virus RNA 1", "Dansoman virus RNA 2", 
                                          "Drosophila C virus", "La Jolla virus",
                                          "Nora virus", "Tobacco mosaic virus"))) %>% 
  mutate(sample_id = factor(sample_id, levels = c("California_2011_2", 
                                                  "NorthCarolina_2006_1", 
                                                  "California_2000_2", 
                                                  "Pennsylvania_1963_1", 
                                                  "Pennsylvania_1963_2", 
                                                  "NewJersey_1942", "Illinois_1930", 
                                                  "Minnesota_1919_3", "Illinois_1908")))

rate_plot <- ggplot(filter(rates2, percent_nucleotide_similarity >= 90), aes(x = virus, y = rate)) +
#  scale_color_viridis() +
  geom_point(aes(color = sample_id), size = 6, alpha = 0.75) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_log10() +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        #        axis.text = element_text(face = "bold"),
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "", y = "Rate (s/n/y)", color = "Sample ID")

rate_plot
ggsave("plots/rates.pdf", units = "in", width = 10, height = 8)
ggsave("plots/rates.svg", units = "in", width = 10, height = 8)
ggsave("plots/rates.jpg", units = "in", width = 10, height = 8)

# rate as a function of time
rate_v_time <- ggplot(filter(rates2, percent_nucleotide_similarity >= 90), aes(x = year, y = rate)) +
  geom_point(aes(color = virus), size = 6, alpha = 0.75) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_log10() +
  scale_x_log10() +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        #        axis.text = element_text(face = "bold"),
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Old Collection Sequence Date (year)", y = "Rate (s/n/y)", color = "Virus Name")

rate_v_time
ggsave("plots/rate_v_time.pdf", units = "in", width = 10, height = 8)
ggsave("plots/rate_v_time.svg", units = "in", width = 10, height = 8)
ggsave("plots/rate_v_time.jpg", units = "in", width = 10, height = 8)

# rate as time, without single virus sequences
rates3 <- rates2 %>% 
  filter(virus %in% c("Galbut virus RNA 1", "Galbut virus RNA 2", 
                      "Galbut virus RNA 3", "Chaq virus", "Vera virus RNA 1", 
                      "Vera virus RNA 2", "Chaq-like virus")) %>% 
  mutate(parent_virus = str_sub(virus, 1, 6),
         parent_virus = str_replace(parent_virus, "Vera v", "Vera virus"),
         parent_virus = str_replace(parent_virus, "Galbut", "Galbut virus"),
         parent_virus = str_replace(parent_virus, "Chaq v", "Galbut virus"),
         parent_virus = str_replace(parent_virus, "Chaq-l", "Vera virus"))
  
rate_v_time2 <- ggplot(rates3, aes(x = year, y = rate)) +
  geom_point(aes(color = virus), size = 8, alpha = 0.85) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_log10() +
  scale_x_reverse() +
  facet_wrap(~ parent_virus) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        #        axis.text = element_text(face = "bold"),
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Old Collection Sequence Date (year)", y = "Rate (s/n/y)", color = "Virus Name")

rate_v_time2
ggsave("plots/Supplemental5/rate_v_time2.pdf", units = "in", width = 10, height = 8)
ggsave("plots/Supplemental5/rate_v_time2.svg", units = "in", width = 10, height = 8)
ggsave("plots/Supplemental5/rate_v_time2.jpg", units = "in", width = 10, height = 8)
