library(tidyverse)
df <- read.csv("metadata/LocationData.csv")

df <- df %>% 
  select(Date.Collected, Concentration..ng.μl., Extraction.Method, Storage.Type) %>% 
  rename(date_collected = Date.Collected,
         concentration = Concentration..ng.μl.,
         method = Extraction.Method)

ggplot(df, aes(x = date_collected)) + 
  geom_point(aes(y = concentration, color = method, shape = factor(Storage.Type)), 
             size=4.25, stroke=0.1, alpha = 0.75) +
  scale_color_manual(values = c("blueviolet", "darkgreen"), guide = "none") +
  facet_wrap(~method) +
  scale_x_reverse() +
  theme_bw(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  labs(x = "Year Collected", y = "Concentration ng/ul", shape = "Strorage Type")
ggsave("plots/CollectionVsConcentration.pdf", units = "in", width = 10, height = 8)
ggsave("plots/CollectionVsConcentration.jpg", units = "in", width = 10, height = 8)
