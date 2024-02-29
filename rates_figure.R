library(tidyverse)
library(viridis)
library(svglite)

rates <- read_csv("metadata/rates_data.csv")

rates2 <- rates %>% 
  mutate(virus = factor(virus, levels = c("Galbut virus RNA 1", "Galbut virus RNA 2", 
                                          "Galbut virus RNA 3", "Chaq virus", 
                                          "Vera virus RNA 1", "Vera virus RNA 2", 
                                          "Chaq-like virus", 
                                          "Drosophila melanogaster sigmavirus", 
                                          "Craigies Hill virus RNA 1", 
                                          "Craigies Hill virus RNA 2", "
                                          Dansoman virus RNA 1", "Dansoman virus RNA 2", 
                                          "Drosophila C virus", "La Jolla virus",
                                          "Nora virus", "Tobacco mosaic virus"))) %>% 
  mutate(sample_id = factor(sample_id, levels = c("California_2011_2", 
                                                  "NorthCarolina_2006_1", 
                                                  "California_2000_2", 
                                                  "Pennsylvania_1963_1", 
                                                  "Pennsylvania_1963_2", 
                                                  "NewJersey_1942", "Illinois_1930", 
                                                  "Minnesota_1919_3", "Illinois_1908")))

rate_plot <- ggplot(rates2, aes(x = virus, y = rate)) +
#  scale_color_viridis() +
  geom_point(aes(color = sample_id), size = 6, alpha = 0.75) +
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
