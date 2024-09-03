# This script is used to generate a figure of the negative and positive control 
# from the experimentally dried fly and mosquito RT-qPCRs

library(tidyverse)
library(readxl)
library(ggthemes)
library(ggpubr)
library(viridis)
library(svglite)

# Read in the cleaned control data
controls <- read_xlsx("tidy_formats/qPCR_control.xlsx")

# cleanup
controls <- controls %>% 
  mutate(target_name = str_replace(target_name, "Galbut", "Galbut virus"),
         target_name = str_replace(target_name, "Rpl33", "RpL32"),
         target_name = str_replace(target_name, "Rpl32", "RpL32"),
         target_name = str_replace(target_name, "RpL32", "RpL32 mRNA"),
         target_name = str_replace(target_name, "LaJolla", "La Jolla"),
         target_name = str_replace(target_name, "La Jolla", "La Jolla virus"),
         target_name = str_replace(target_name, "Nora", "Nora virus"),
         target_name = str_replace(target_name, "Thika", "Thika virus"),
         target_name = str_replace(target_name, "ACT", "Actin"),
         target_name = str_replace(target_name, "Actin", "Actin mRNA"),
         target_name = str_replace(target_name, "Verdadero", "Verdadero virus"),
         target_name = str_replace(target_name, "Rennavirus", "Guadeloupe mosquito virus"),
         ct = str_replace(ct, "Undetermined", "0.0")) %>% 
  unite("control_type", type:status, na.rm = TRUE, remove = FALSE) %>% 
  mutate(control_type = str_replace(control_type, "cDNA_Pos", "cDNA Positive"),
         control_type = str_replace(control_type, "cDNA_Neg", "cDNA Negative"),
         control_type = str_replace(control_type, "ExtCtrl", "Extraction Negative"),
         control_type = str_replace(control_type, "qPCR", "RT-qPCR Negative"),
         positive = str_replace(positive, "N", "No"),
         positive = str_replace(positive, "Y", "Yes"))

# fly visualization
fly_ctrl <- ggplot(filter(controls, target_name %in% c("Galbut virus", "Nora virus", 
                                                  "Thika virus", "La Jolla virus",
                                                  "RpL32 mRNA"))) +
  geom_point(aes(x = as.numeric(week), y = as.numeric(ct), 
                 color = factor(control_type, 
                                levels = c("cDNA Positive", "cDNA Negative", 
                                           "RT-qPCR Negative", "Extraction Negative")), 
                 shape = positive), size = 3.5, 
             stroke = 0.1, alpha = 0.6) +
  scale_color_manual(values = c("darkred", "lightskyblue1", "skyblue4", "paleturquoise3")) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        #        axis.title = element_text(face = "bold"),
        #        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  facet_wrap(~factor(target_name, levels = c("RpL32 mRNA", "Galbut virus", 
                                        "Nora virus", "Thika virus", "La Jolla virus")),
             ncol = 1) +
  labs(x = "Weeks of Storage", 
       y = "Cycles to Threshold Detection (ct)", 
       shape = "RT-qPCR \nPositive Status", color = "Type of control")

fly_ctrl
# remove # to save plot
ggsave("plots/Figure1/fly_controls.pdf", units = "in", width = 10, height = 8)
ggsave("plots/Figure1/fly_controls.svg", units = "in", width = 10, height = 8)
ggsave("plots/Figure1/fly_controls.jpg", units = "in", width = 10, height = 8)

# mosquito visualization
mos_ctrl <- ggplot(filter(controls, target_name %in% 
                            c("Actin mRNA", "Verdadero virus", 
                              "Guadeloupe mosquito virus"))) +
  geom_point(aes(x = as.numeric(week), y = as.numeric(ct), 
                 color = factor(control_type, 
                                levels = c("cDNA Positive", "cDNA Negative", 
                                           "RT-qPCR Negative", "Extraction Negative")), 
                 shape = positive), size = 3.5, 
             stroke = 0.1, alpha = 0.6) +
  scale_color_manual(values = c("darkred", "lightskyblue1", "skyblue4", "paleturquoise3")) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        #        axis.title = element_text(face = "bold"),
        #        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  facet_wrap(~factor(target_name, levels = c("Actin mRNA", "Verdadero virus", 
                                             "Guadeloupe mosquito virus")),
             ncol = 1) +
  labs(x = "Weeks of Storage", 
       y = "Cycles to Threshold Detection (ct)", 
       shape = "RT-qPCR \nPositive Status", color = "Type of control")

mos_ctrl
# remove # to save plot
ggsave("plots/Figure1/mos_controls.pdf", units = "in", width = 10, height = 8)
ggsave("plots/Figure1/mos_controls.svg", units = "in", width = 10, height = 8)
ggsave("plots/Figure1/mos_controls.jpg", units = "in", width = 10, height = 8)
