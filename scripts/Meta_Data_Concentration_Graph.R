#This script reads in metadata from an excel sheet with information on the yield of RNA
# from museum specimens and generates two visualizations. It also performs a MLR analysis 
# to determine if there is a statistically significant difference in RNA cocnentration 
# based on time and whether the Hawaii samples are removed.

library(tidyverse)
library(svglite)

df <- read.csv("metadata/LocationData.csv")

df2 <- df %>% 
  select(Date.Collected, Concentration..ng.μl., Extraction.Method, Storage.Type) %>% 
  rename(date_collected = Date.Collected,
         concentration = Concentration..ng.μl.,
         method = Extraction.Method)

OC_conc <- ggplot(df2, aes(x = date_collected)) + 
  geom_point(aes(y = concentration, color = method, shape = factor(Storage.Type)), 
             size=6, stroke=0.1, alpha = 0.75) +
  scale_color_manual(values = c("blueviolet", "darkgreen"), guide = "none") +
  facet_wrap(~method) +
  scale_x_reverse() +
  theme_bw(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        #        axis.title = element_text(face = "bold"),
        #        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  labs(x = "Year Collected", y = "Concentration ng/ul", shape = "Storage Type")

OC_conc
ggsave("plots/Figure2/CollectionVsConcentration.pdf", units = "in", width = 10, height = 8)
ggsave("plots/Figure2/CollectionVsConcentration.jpg", units = "in", width = 10, height = 8)
ggsave("plots/Figure2/CollectionVsConcentration.svg", units = "in", width = 10, height = 8)

df_mlr <- lm(concentration ~ method + Storage.Type, data = df2)
summary(df_mlr)

# Hawaii vs everything else
df3 <- df %>% 
  select(Date.Collected, Concentration..ng.μl., Extraction.Method, Storage.Type,
         Location) %>% 
  rename(date_collected = Date.Collected,
         concentration = Concentration..ng.μl.,
         method = Extraction.Method,
         location = Location) %>% 
  mutate(HI = str_detect(df$Location, "HI"),
         HI = str_replace(HI, "TRUE", "Hawai'i"),
         HI = str_replace(HI, "FALSE", "Other"))

HI_vs_other <- ggplot(df3, aes(x = date_collected, y = concentration)) +
  geom_point(aes(color = method, shape = factor(Storage.Type)), size = 6, 
             stroke = 0.1, alpha = 0.75) +
  scale_color_manual(values = c("blueviolet", "darkgreen"), guide = "none") +
  facet_wrap(~HI) +
  scale_x_reverse() +
  theme_bw(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        #        axis.title = element_text(face = "bold"),
        #        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  labs(x = "Year Collected", y = "Concentration ng/ul", shape = "Storage Type")
 
HI_vs_other  
ggsave("plots/Supplemental4/HIvsOther.pdf", units = "in", width = 10, height = 8)
ggsave("plots/Supplemental4/HIvsOther.svg", units = "in", width = 10, height = 8)
ggsave("plots/Supplemental4/HIvsOther.jpg", units = "in", width = 10, height = 8)

df4 <- df3 %>% 
  filter(HI == "Other")

df_mlr_no_hawaii <- lm(concentration ~ method + Storage.Type, data = df4)
summary(df_mlr_no_hawaii)
