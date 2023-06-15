library(tidyverse)
df <- read.csv("metadata/LocationData.csv")

df <- df %>% 
  select(Date.Collected, Concentration..ng.μl., Extraction.Method) %>% 
  rename(date_collected = Date.Collected,
         concentration = Concentration..ng.μl.,
         method = Extraction.Method)

ggplot(df, aes(x = date_collected)) + 
  geom_point(aes(y = concentration, fill = method), shape=21, size=2, stroke=0.1, color="black") +
  scale_fill_manual(values = c("blueviolet", "darkgreen")) +
  #geom_line(aes(y=Concentration..ng.ul.) + 
  theme_bw(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white")) +
  labs(x = "Year Collected", y = "Concentration ng/ul", fill = 'Extraction Method')
ggsave("plots/CollectionVsConcentration.pdf", units = "in", width = 10, height = 8)
