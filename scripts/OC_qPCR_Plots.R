# This script is used to generate a figure of the negative and positive control 
# from the experimentally dried fly and mosquito RT-qPCRs

library(tidyverse)
library(readxl)
library(ggthemes)
library(ggpubr)
library(viridis)
library(svglite)

# Load OC qPCR Data
OC <- read_xlsx("qPCR_data/OC_qPCR_tidy.xlsx")


# cleanup
OC <- OC %>% 
  mutate(TargeT_Name = str_replace(TargeT_Name, "galbut", "Galbut virus"),
         TargeT_Name = str_replace(TargeT_Name, "Rpl", "RpL32 mRNA"),
         overall_positive = str_replace(overall_positive, "N", "No"),
         overall_positive = str_replace(overall_positive, "Y", "Yes")) 

# OC data
OC_fig <- ggplot(filter(OC, status == "sample")) +
  geom_point(aes(x = as.numeric(date_collected), y = CT, color = overall_positive, 
             shape = overall_positive), size = 6, stroke = 0.1, alpha = 0.7) +
  scale_color_manual(values = c("blue3", "olivedrab4")) +
  scale_x_reverse() +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  facet_wrap(~TargeT_Name) +
  labs(x = "Year Collected", 
       y = "Cycles to Threshold Detection (ct)", 
       color = "RT-qPCR \nPositive", shape = "RT-qPCR \nPositive")


OC_fig
ggsave("plots/Figure2/OC_qPCR.pdf", units = "in", width = 10, height = 8)
ggsave("plots/Figure2/OC_qPCR.svg", units = "in", width = 10, height = 8)
ggsave("plots/Figure2/OC_qPCR.jpg", units = "in", width = 10, height = 8)

OC_ctrs <- OC %>% 
  filter(status == "control")

# OC qPCR Controls
OC_ctrl <- ggplot(filter(OC, status == "control")) +
  geom_point(aes(x = step, y = CT, color = control_type, 
                 shape = overall_positive), size = 6, stroke = 0.1, alpha = 0.7) +
  scale_color_manual(values = c("skyblue4", "darkred")) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  facet_wrap(~TargeT_Name) +
  labs(x = "Protocol Step", 
       y = "Cycles to Threshold Detection (ct)", 
       color = "Type of Control", shape = "RT-qPCR \nPositive")

 OC_ctrl
 ggsave("plots/Figure2/OC_qPCR_ctrl.pdf", units = "in", width = 10, height = 8)
 ggsave("plots/Figure2/OC_qPCR_ctrl.svg", units = "in", width = 10, height = 8)
 ggsave("plots/Figure2/OC_qPCR_ctrl.jpg", units = "in", width = 10, height = 8)
 