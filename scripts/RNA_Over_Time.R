# load libraries
library(tidyverse)
library(readxl)
library(ggthemes)
library(writexl)
library(ggrepel)
library(ggpmisc)
library(rstatix)
library(ggpubr)
library(viridis)
library(svglite)
library(broom)
library(knitr)
library(gt)
library(webshot2)
library(performance)

# read in data
fly <- read_xlsx("tidy_formats/fly_data2.xlsx")

mosquito <- read_xlsx("tidy_formats/mosquito_data2.xlsx")

# fly data figures

# fresh for normalizing
fresh <- fly %>% 
  filter(group == "Fresh",
         present == "y",
         week == 2) %>%
  select(week, target, ct) %>% 
  type_convert() %>% 
  group_by(week, target) %>% 
  summarise(mean_ct_fresh = mean(ct))

# pre-processing 
fly <- fly %>%
  filter(present == "y") %>% 
  type_convert() %>% 
  group_by(target, group, week, length) %>%
  mutate(mean_ct = mean(ct, na.rm = TRUE),
         sd_ct = sd(ct, na.rm = TRUE),
         sample_length = paste(target, length, sep = "_")) 

fly_data2 <- left_join(fly, fresh, by = join_by(target)) %>% 
  mutate(delta_ct = mean_ct_fresh - ct)

fly_data3 <- fly_data2 %>% 
  mutate(target = str_replace(target, "Galbut", "Galbut virus"),
         target = str_replace(target, "RpL32", "RpL32 mRNA"),
         target = str_replace(target, "La Jolla", "La Jolla virus"),
         target = str_replace(target, "Nora", "Nora virus"),
         target = str_replace(target, "Thika", "Thika virus"))

# get mean fold change, remove fresh & short samples
fly_fc <- fly_data3 %>% 
  rename(week = week.x) %>% 
  group_by(target, group, week) %>% 
  filter(length != "short",
         group != "Fresh") %>% 
  mutate(fold_change = 2^-(delta_ct),
         mean_fc = mean(fold_change, na.rm = TRUE),
         sd_fc = sd(fold_change, na.rm = TRUE),
         mean_dct = mean(delta_ct, na.rm = TRUE),
         sd_dct = sd(delta_ct, na.rm = TRUE)) 

# delta ct - three targets
rel_fly_dct <- ggplot(filter(fly_fc, target %in% c("Galbut virus", "Nora virus", "RpL32 mRNA")),
       aes(x = week)) +
  geom_point(aes(y = delta_ct, fill = group), shape = 21, size = 1.35, 
             stroke = 0.1, color = "black", alpha = 0) +
  stat_compare_means(aes(y= delta_ct, group = group), label = "p.signif", 
                     hide.ns = TRUE, label.y = 7, size = 3.5, alpha = 0.75, symnum.args = 
                       list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                            symbols = c("D", "C", "B", "A", "ns"))) +
  geom_line(aes(y = mean_dct, group = group, linetype = group, colour = group), 
            linewidth = 0.75, alpha = 1) +
  scale_fill_manual(values = c("firebrick", "navyblue")) +
  scale_colour_manual(values = c("firebrick", "navyblue")) +
  geom_errorbar(aes(ymin = (mean_dct - sd_dct), ymax = (mean_dct + sd_dct), color = group), 
                width = 1, alpha = 0.25) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        #        axis.title = element_text(face = "bold"),
        #        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  facet_wrap(~factor(target, levels = c("RpL32 mRNA", "Galbut virus", 
                                        "Nora virus")), ncol = 1) +
  labs(x = "Weeks After Sample Storage", 
       y = "Delta Ct Relative to Time Point 0 Fresh Fly", 
       fill = "Sample Storage", linetype = "Sample Storage", colour = "Sample Storage")

rel_fly_dct
# remove # to save plot
ggsave("plots/Figure1/Relative_fly_dct_3targets.pdf", units = "in", width = 10, height = 8)
ggsave("plots/Figure1/Relative_fly_dct_3targets.svg", units = "in", width = 10, height = 8)
ggsave("plots/Figure1/Relative_fly_dct_3targets.jpg", units = "in", width = 10, height = 8)

rel_fly_supp <- ggplot(filter(fly_fc, target %in% c("La Jolla virus", "Thika virus")), aes(x = week)) +
  geom_point(aes(y = delta_ct, fill = group), shape = 21, size = 1.35, 
             stroke = 0.1, color = "black", alpha = 0) +
  stat_compare_means(aes(y= delta_ct, group = group), label = "p.signif", 
                     hide.ns = TRUE, label.y = 7, size = 3.5, alpha = 0.75, symnum.args = 
                       list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                            symbols = c("D", "C", "B", "A", "ns"))) +
  geom_line(aes(y = mean_dct, group = group, linetype = group, colour = group), 
            linewidth = 1, alpha = 1) +
  scale_fill_manual(values = c("firebrick", "navyblue")) +
  scale_colour_manual(values = c("firebrick", "navyblue")) +
  geom_errorbar(aes(ymin = (mean_dct - sd_dct), ymax = (mean_dct + sd_dct), color = group), 
                width = 1, alpha = 0.25) +  
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        #        axis.title = element_text(face = "bold"),
        #        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  facet_wrap(~factor(target, levels = c("La Jolla virus", "Thika virus")), 
             ncol = 1, scales = "free_y") +
  labs(x = "Weeks After Sample Storage", 
       y = "Delta Ct Relative to Time Point 0 Fresh Fly", 
       fill = "Sample Storage", linetype = "Sample Storage", colour = "Sample Storage")

rel_fly_supp
ggsave("plots/Supplemental1/Relative_fly_dct_SuppTargets.pdf", units = "in", width = 10, height = 8)
ggsave("plots/Supplemental1/Relative_fly_dct_SuppTargets.svg", units = "in", width = 10, height = 8)
ggsave("plots/Supplemental1/Relative_fly_dct_SuppTargets.jpg", units = "in", width = 10, height = 8)

# short vs long Galbut & Rpl
short_v_long <- fly_data3 %>% 
  rename(week = week.x) %>% 
  filter(target %in% c("Galbut virus", "RpL32 mRNA")) %>% 
  filter(group == "Dry") %>% 
  select(week, group, target, length, mean_ct, sd_ct) %>% 
  group_by(week, target, group, length, mean_ct, sd_ct) %>% 
  summarise() 

short_v_long_wide <- short_v_long %>% 
  pivot_wider(names_from = length, values_from = c(mean_ct, sd_ct)) %>% 
  ungroup()

#Plot min and max
plot_min_x <- 16
plot_max_x <- 30
plot_min_y <- 16
plot_max_y <- 30

long_vs_short <- ggplot(short_v_long_wide) +
  scale_fill_viridis() +
  geom_point(aes(x = mean_ct_long, y = mean_ct_short, fill = week), shape = 21, 
             size = 4, stroke = 0.25, alpha = 0.85) +
  geom_errorbarh(aes(xmin = if_else((mean_ct_long - sd_ct_long) < plot_min_x, 
                     plot_min_x,(mean_ct_long - sd_ct_long)), 
                 xmax = if_else((mean_ct_long + sd_ct_long) > plot_max_x,
                                plot_max_x,(mean_ct_long + sd_ct_long)), 
                 y = mean_ct_short), 
                height = 0.25, color = "slategray3", linewidth = 0.25, alpha = 0.75) +
  geom_errorbar(aes(ymin = if_else((mean_ct_short - sd_ct_short) < plot_min_y, 
                                   plot_min_y,(mean_ct_short - sd_ct_short)), 
                    ymax = if_else((mean_ct_short + sd_ct_short) > plot_max_y, 
                                   plot_max_y,(mean_ct_short + sd_ct_short)), 
                    x = mean_ct_long), 
                width = 0.25, color = "slategray3", linewidth = 0.25, alpha = 0.75) +
  facet_wrap(~factor(target, levels = c("RpL32 mRNA", "Galbut virus"))) +
  stat_poly_line(aes(x = mean_ct_long, mean_ct_short), method = "lm", alpha = 0.5, 
                 se = FALSE, colour = "slategray3", linetype = "dotdash", linewidth = 0.5) +
  stat_poly_eq(aes(x = mean_ct_long, mean_ct_short, label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*"))) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        #        axis.title = element_text(face = "bold"),
        #        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  coord_fixed(xlim = c(plot_min_x, plot_max_x), ylim = c(plot_min_y, plot_max_y)) +
  geom_text_repel(aes(x = mean_ct_long, y = mean_ct_short, label = week), 
                  size = 4, colour = "grey30") +
  geom_abline(intercept = 0, slope = 1, color="grey40", alpha=0.5, linewidth=0.5,
              linetype=2) +
  labs(x = "Mean Ct at Each Time Point (long primer)", 
       y = "Mean Ct at Each Time Point \n(short primer)", fill = "Week")

long_vs_short
# remove # to save plot
ggsave("plots/Figure1/long_vs_short.pdf", units = "in", width = 10, height = 8)
ggsave("plots/Figure1/long_vs_short.svg", units = "in", width = 10, height = 8)
ggsave("plots/Figure1/long_vs_short.jpg", units = "in", width = 10, height = 8)

# dry vs frozen fly
plot_min_x_df <- 16
plot_max_x_df <- 33
plot_min_y_df <- 16
plot_max_y_df <- 33

dry_v_frozen <- fly_data3 %>% 
  rename(week = week.x) %>% 
  filter(group != "Fresh") %>% 
  filter(length == "long") %>% 
  select(week, group, target, mean_ct, sd_ct) %>% 
  group_by(week, target, group, mean_ct, sd_ct) %>% 
  summarise() 

dry_v_frozen_wide <- dry_v_frozen %>% 
  pivot_wider(names_from = group, values_from = c(mean_ct, sd_ct)) %>% 
  ungroup()

dry_vs_frozen_fly <- ggplot(dry_v_frozen_wide) +
  scale_fill_viridis() +
  geom_point(aes(x = mean_ct_Dry, y = mean_ct_Frozen, fill = week), shape = 21, 
             size = 4, stroke = 0.25, alpha = 0.85) + 
  geom_errorbarh(aes(xmin = mean_ct_Dry - sd_ct_Dry, 
                     xmax = mean_ct_Dry + sd_ct_Dry, y = mean_ct_Frozen), 
                 height = 0.25, color = "slategray3", linewidth = 0.25, alpha = 0.75) +
  geom_errorbar(aes(ymin = mean_ct_Frozen - sd_ct_Frozen, 
                    ymax = mean_ct_Frozen + sd_ct_Frozen, x = mean_ct_Dry), 
                width = 0.25, color = "slategray3", linewidth = 0.25, alpha = 75) +
  facet_wrap(~factor(target, levels = c("RpL32 mRNA", "Galbut virus", "La Jolla virus", 
                                        "Nora virus", "Thika virus")), nrow = 2) +
  stat_poly_line(aes(x = mean_ct_Dry, mean_ct_Frozen), method = "lm", alpha = 0.5, 
                 se = FALSE, colour = "slategray3", linetype = "dotdash", 
                 linewidth = 0.5) +
  stat_poly_eq(aes(x = mean_ct_Dry, mean_ct_Frozen, 
                   label = paste(after_stat(eq.label), 
                                 after_stat(rr.label), sep = "*\", \"*"))) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        #        axis.title = element_text(face = "bold"),
        #        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  coord_fixed(xlim = c(plot_min_x_df, plot_max_x_df), 
              ylim = c(plot_min_y_df, plot_max_y_df)) +
  geom_text_repel(aes(x = mean_ct_Dry, y = mean_ct_Frozen, label = week), 
                  size = 4, colour = "grey30") +
  geom_abline(intercept = 0, slope = 1, color="grey40", alpha=0.5, linewidth=0.5,
              linetype=2) +
  labs(x = "Mean Ct at Each Time Point (Dry)", 
       y = "Mean Ct at Each Time Point (Frozen)", fill = "Week")

dry_vs_frozen_fly
# remove # to save plot
ggsave("plots/Supplemental1/dry_vs_frozen.pdf", units = "in", width = 10, height = 8)
ggsave("plots/Supplemental1/dry_vs_frozen.svg", units = "in", width = 10, height = 8)
ggsave("plots/Supplemental1/dry_vs_frozen.jpg", units = "in", width = 10, height = 8)


# Mosquito figures

# Week 4 for normalizing
w4 <- mosquito %>% 
  filter(group == "Frozen",
         present == "y",
         week == 4) %>%
  select(week, target, ct) %>% 
  type_convert() %>% 
  group_by(week, target) %>% 
  summarise(mean_ct_4w = mean(ct))

# pre-processing mosquito plots
mosquito <- mosquito %>% 
  filter(present == "y") %>% 
  type_convert() %>% 
  group_by(target, group, week) %>%
  mutate(mean_ct = mean(ct, na.rm = TRUE),
         sd_ct = sd(ct, na.rm = TRUE),
         sample_group = paste0(target)) %>% 
  ungroup()

mos_data2 <- left_join(mosquito, w4, by = join_by(target)) %>% 
  mutate(delta_ct = mean_ct_4w - ct) %>% 
  rename(week = week.x)

# get mean fold change
mos_data3 <- mos_data2 %>% 
  group_by(target, group, week) %>% 
  mutate(fold_change = 2^-(delta_ct),
         mean_fc = mean(fold_change, na.rm = TRUE),
         sd_fc = sd(fold_change, na.rm = TRUE),
         mean_dct = mean(delta_ct, na.rm = TRUE),
         sd_dct = sd(delta_ct, na.rm = TRUE),
         target = str_replace(target, "Verdadero", "Verdadero virus"),
         target = str_replace(target, "Rennavirus", "Guadeloupe mosquito virus"),
         target = str_replace(target, "Actin", "Actin mRNA"))

# delta ct fold change
rel_mos_dct <- ggplot(mos_data3, aes(x = week)) +
  geom_point(aes(y = delta_ct, fill = group), shape = 21, size = 1.35, 
             stroke = 0.1, color = "black", alpha = 0) +
  stat_compare_means(aes(y= delta_ct, group = group), label = "p.signif",
                     hide.ns = TRUE, label.y = 7, size = 3.5, alpha = 0.85, symnum.args = 
                       list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                            symbols = c("D", "C", "B", "A", "ns"))) +
  geom_line(aes(y = mean_dct, group = group, linetype = group, colour = group), 
            linewidth = 1, alpha = 1) +
  scale_fill_manual(values = c("firebrick", "navyblue")) +
  scale_colour_manual(values = c("firebrick", "navyblue")) +
  geom_errorbar(aes(ymin = (mean_dct - sd_dct), ymax = (mean_dct + sd_dct),
                    color = group), width = 1, alpha = 0.25) +  
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        #        axis.title = element_text(face = "bold"),
        #        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  facet_wrap(~factor(target, levels = c("Actin mRNA", "Verdadero virus", 
                                        "Guadeloupe mosquito virus")), ncol = 1) +
  labs(x = "Weeks After Sample Storage", 
       y = "Delta Ct Relative to Time Point \n4 Week Frozen Mosquito",
       linetype = "Sample Storage",
       fill = "Sample Storage",
       colour = "Sample Storage") 

rel_mos_dct  
# remove # to save plot
ggsave("plots/Supplemental2/Relative_mosquito_dct.pdf", units = "in", width = 10, height = 8)
ggsave("plots/Supplemental2/Relative_mosquito_dct.svg", units = "in", width = 10, height = 8)
ggsave("plots/Supplemental2/Relative_mosquito_dct.jpg", units = "in", width = 10, height = 8)

# dry vs frozen mos
plot_min_x_df_m <- 14
plot_max_x_df_m <- 30
plot_min_y_df_m <- 14
plot_max_y_df_m <- 30

dry_v_frozen_mos <- mos_data3 %>%
  select(week, group, target, mean_ct, sd_ct) %>% 
  group_by(week, target, group, mean_ct, sd_ct) %>% 
  summarise() 

dry_v_frozen_wide_mos <- dry_v_frozen_mos %>% 
  pivot_wider(names_from = group, values_from = c(mean_ct, sd_ct)) %>% 
  ungroup()

dry_vs_frozen_mos <- ggplot(dry_v_frozen_wide_mos) +
  scale_fill_viridis() +
  geom_point(aes(x = mean_ct_Dry, y = mean_ct_Frozen, fill = week), shape = 21, 
             size = 4, stroke = 0.25, alpha = 0.85) +
  geom_errorbarh(aes(xmin = mean_ct_Dry - sd_ct_Dry, 
                     xmax = mean_ct_Dry + sd_ct_Dry, y = mean_ct_Frozen), 
                 height = 0.25, color = "slategray3", linewidth = 0.25, alpha = 0.75) +
  geom_errorbar(aes(ymin = mean_ct_Frozen - sd_ct_Frozen, 
                    ymax = mean_ct_Frozen + sd_ct_Frozen, 
                    x = mean_ct_Dry), 
                width = 0.25, color = "slategray3", linewidth = 0.25, alpha = 0.75) +
  facet_wrap(~factor(target, levels = c("Actin mRNA", "Verdadero virus", 
                                        "Guadeloupe mosquito virus")), ncol = 2) +
  stat_poly_line(aes(x = mean_ct_Dry, mean_ct_Frozen), method = "lm", alpha = 0.5, 
                 se = FALSE, colour = "slategray3", linetype = "dotdash", linewidth = 0.5) +
  stat_poly_eq(aes(x = mean_ct_Dry, mean_ct_Frozen, 
                   label = paste(after_stat(eq.label), 
                                 after_stat(rr.label), sep = "*\", \"*"))) +
  coord_fixed(xlim = c(plot_min_x_df_m, plot_max_x_df_m), 
              ylim = c(plot_min_y_df_m, plot_max_y_df_m)) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        #        axis.title = element_text(face = "bold"),
        #        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  geom_text_repel(aes(x = mean_ct_Dry, y = mean_ct_Frozen, label = week), 
                  size = 4, colour = "grey30") +
  geom_abline(intercept = 0, slope = 1, color="grey40", alpha=0.5, linewidth=0.5,
              linetype=2) +
  labs(x = "Mean Ct at Each Time Point (Dry)", 
       y = "Mean Ct at Each Time Point (Frozen)", fill = "Week")

dry_vs_frozen_mos
# remove # to save plot
ggsave("plots/Supplemental2/dry_vs_frozen_mos.pdf", units = "in", width = 10, height = 8)
ggsave("plots/Supplemental2/dry_vs_frozen_mos.svg", units = "in", width = 10, height = 8)
ggsave("plots/Supplemental2/dry_vs_frozen_mos.jpg", units = "in", width = 10, height = 8)

# Number Positive                                                                                                                                                                                                                                                                             
fly_summary <- fly %>% 
  group_by(target, group, week, length) %>% 
  count(present) %>% 
  filter(group != "Fresh") %>% 
  mutate(target = str_replace(target, "Galbut", "Galbut virus"),
         target = str_replace(target, "RpL32", "RpL32 mRNA"),
         target = str_replace(target, "La Jolla", "La Jolla virus"),
         target = str_replace(target, "Nora", "Nora virus"),
         target = str_replace(target, "Thika", "Thika virus"))

fly_summary$week <- as.numeric(fly_summary$week)


# mosquito targets
mos_summary <- mosquito %>% 
  group_by(target, group, week) %>% 
  count(present) %>% 
  mutate(target = str_replace(target, "Actin", "Actin mRNA"),
         target = str_replace(target, "Rennavirus", "Guadeloupe mosquito virus"),
         target = str_replace(target, "Verdadero", "Verdadero virus"),
         length = "long")

mos_summary$week <- as.numeric(mos_summary$week)

all <- rbind(fly_summary, mos_summary)

all <- all %>% 
  mutate(length = str_replace(length, "long", "Long"),
         length = str_replace(length, "short", "Short"))


n_pos_fly <- ggplot(filter(all, target == "RpL32 mRNA")) +
  geom_point(aes(x = week, y = n, color = group), 
             size = 5, stroke = 0.25, alpha = 0.75) +
  scale_color_manual(values = c("firebrick", "navyblue")) +
  facet_grid(group ~ length, scales = ) +
  ylim(0,6) +
  theme_few(base_size = 11) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        #        axis.title = element_text(face = "bold"),
        #        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  labs(x = "Weeks After Sample Storage", y = "Number Samples Positive by RT-qPCR",
       color = "Sample Storage", shape = "Primer Length")

n_pos_fly
ggsave("plots/Figure1/n_positve_fly_mRNA.pdf", units = "in", width = 10, height = 8)
ggsave("plots/Figure1/n_positve_fly_mRNA.jpg", units = "in", width = 10, height = 8)
ggsave("plots/Figure1/n_positve_fly_mRNA.svg", units = "in", width = 10, height = 8)

n_pos_mos <- ggplot(filter(all, target == "Actin mRNA")) +
  geom_point(aes(x = week, y = n, color = group), 
             size = 5, stroke = 0.25, alpha = 0.75) +
  scale_color_manual(values = c("firebrick", "navyblue")) +
  ylim(0, 6) +
  facet_grid(group ~ length) +
  theme_few(base_size = 11) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        #        axis.title = element_text(face = "bold"),
        #        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  labs(x = "Weeks After Sample Storage", y = "Number Samples Positive by RT-qPCR",
       color = "Sample Storage", shape = "Primer Length")

n_pos_mos
ggsave("plots/Supplemental2/n_positve_mos_mRNA.pdf", units = "in", width = 10, height = 8)
ggsave("plots/Supplemental2/n_positve_mos_mRNA.jpg", units = "in", width = 10, height = 8)
ggsave("plots/Supplemental2/n_positve_mos_mRNA.svg", units = "in", width = 10, height = 8)

# RNA Concentrations Over Time
conc <- read_csv("metadata/OverTimeRNAConcentrations.csv")

ggplot(conc) +
  geom_point(aes(x = week, y = concentration, color = week)) +
  facet_grid(organism~storage)

# Concentration Figure
y_axis <- expression(paste(" Mean Concentration (ng/",mu,"l)"))

conc_fly <- ggplot(filter(conc2, organism == "D. melanogaster")) +
  geom_point(aes(x = week, y = mean_conc, color = storage), size = 6, 
             alpha = 0.75) +
  scale_color_manual(values = c("firebrick", "navyblue")) +
  geom_errorbar(aes(x = week, ymin = (mean_conc - sd_conc), ymax = (mean_conc + sd_conc), color = storage), 
                width = 1.5, alpha = 0.25) + 
#  facet_wrap(~factor(organism, levels = c("D. melanogaster", "Ae. aegypti")), scales = "free_x") +
  theme_few(base_size = 11) +
  theme_minimal(base_size = 11) +
  labs(x = "Weeks After Sample Storage", y = y_axis,
       color = "Sample Storage") +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        #        axis.title = element_text(face = "bold"),
        #        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) 

conc_fly
ggsave("plots/Figure1/Fly_RNA_concentrations.pdf", units = "in", width = 10, height = 8)
ggsave("plots/Figure1/Fly_RNA_concentrations.jpg", units = "in", width = 10, height = 8)
ggsave("plots/Figure1/Fly_RNA_concentrations.svg", units = "in", width = 10, height = 8)

# Mosquito Concentration
conc_mos <- ggplot(filter(conc2, organism == "Ae. aegypti")) +
  geom_point(aes(x = week, y = mean_conc, color = storage), size = 6, 
             alpha = 0.75) +
  scale_color_manual(values = c("firebrick", "navyblue")) +
  geom_errorbar(aes(x = week, ymin = (mean_conc - sd_conc), ymax = (mean_conc + sd_conc), color = storage), 
                width = 1.5, alpha = 0.25) + 
#  facet_wrap(~factor(organism, levels = c("D. melanogaster", "Ae. aegypti")), scales = "free_x") +
  theme_few(base_size = 11) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
#        axis.title = element_text(face = "bold"),
#        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  labs(x = "Weeks After Sample Storage", y = y_axis,
       color = "Sample Storage")

conc_mos
ggsave("plots/Supplemental2/Mos_RNA_concentrations.pdf", units = "in", width = 10, height = 8)
ggsave("plots/Supplemental2/Mos_RNA_concentrations.jpg", units = "in", width = 10, height = 8)
ggsave("plots/Supplemental2/Mos_RNA_concentrations.svg", units = "in", width = 10, height = 8)



