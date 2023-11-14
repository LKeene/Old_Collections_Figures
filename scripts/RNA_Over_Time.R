# load libraries
library(tidyverse)
library(readxl)
library(ggthemes)
library(writexl)
library(ggrepel)
library(ggpmisc)
library(rstatix)
library(ggpubr)

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

# delta ct fold change- three targets
ggplot(filter(fly_fc, target %in% c("Galbut virus", "Nora virus", "RpL32 mRNA")),
       aes(x = week)) +
  geom_point(aes(y = delta_ct, fill = group), shape = 21, size = 1.35, 
             stroke = 0.1, color = "black", alpha = 0.75) +
  stat_compare_means(aes(y= delta_ct, group = group), label = "p.signif", 
                     hide.ns = TRUE, label.y = 7, size = 6.5) +
  geom_line(aes(y = mean_dct, group = group, linetype = group, colour = group), 
            linewidth = 0.75, alpha = 0.75) +
  scale_fill_manual(values = c("firebrick", "navyblue")) +
  scale_colour_manual(values = c("firebrick", "navyblue")) +
  geom_errorbar(aes(ymin = (mean_dct - sd_dct), ymax = (mean_dct + sd_dct)), 
                width = 1, color = "grey50", alpha = 0.4) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  facet_wrap(~factor(target, levels = c("Galbut virus", "Nora virus",
                                        "RpL32 mRNA")), ncol = 1) +
  labs(x = 'Weeks After Collection', 
       y = "Log(2) Fold Change Relative to Time Point 0 Fresh FoCo-17", 
       fill = "Sample Storage", linetype = "Sample Storage", colour = "Sample Storage")

# remove # to save plot
ggsave("plots/Relative_fly_dct_3targets.pdf", units = "in", width = 10, height = 8)
ggsave("plots/Relative_fly_dct_3targets.svg", units = "in", width = 10, height = 8)

ggplot(filter(fly_fc, target %in% c("La Jolla virus", "Thika virus")), aes(x = week)) +
  geom_point(aes(y = delta_ct, fill = group), shape = 21, size = 1.35, 
             stroke = 0.1, color = "black", alpha = 0.75) +
  stat_compare_means(aes(y= delta_ct, group = group), label = "p.signif", hide.ns = TRUE, label.y = 7, size = 6.5) +
  geom_line(aes(y = mean_dct, group = group, linetype = group, colour = group), linewidth = 0.75, alpha = 0.75) +
  scale_fill_manual(values = c("firebrick", "navyblue")) +
  scale_colour_manual(values = c("firebrick", "navyblue")) +
  geom_errorbar(aes(ymin = (mean_dct - sd_dct), ymax = (mean_dct + sd_dct)), 
                width = 1, color = "grey50", alpha = 0.4) +  
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  facet_wrap(~factor(target, levels = c("La Jolla virus", "Thika virus")), ncol = 1, scales = "free_y") +
  labs(x = 'Weeks After Collection', 
       y = "Log(2) Fold Change Relative to Time Point 0 Fresh FoCo-17", 
       fill = "Sample Storage", linetype = "Sample Storage", colour = "Sample Storage")

ggsave("plots/Relative_fly_dct_SuppTargets.pdf", units = "in", width = 10, height = 8)
ggsave("plots/Relative_fly_dct_SuppTargets.svg", units = "in", width = 10, height = 8)
ggsave("plots/Relative_fly_dct_SuppTargets.jpg", units = "in", width = 10, height = 8)

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

short_v_long_wide$week <- as.factor(short_v_long_wide$week)

#Plot min and max
plot_min_x <- 16
plot_max_x <- 30
plot_min_y <- 16
plot_max_y <- 30

ggplot(short_v_long_wide) +
  geom_point(aes(x = mean_ct_long, y = mean_ct_short, fill = week), shape = 21, 
             size = 4, stroke = 0.25, alpha = 0.75) +
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
  facet_wrap(~ target) +
  stat_poly_line(aes(x = mean_ct_long, mean_ct_short), method = "lm", alpha = 0.5, 
                 se = FALSE, colour = "slategray3", linetype = "dotdash", linewidth = 0.5) +
  stat_poly_eq(aes(x = mean_ct_long, mean_ct_short, label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*"))) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  coord_fixed(xlim = c(plot_min_x, plot_max_x), ylim = c(plot_min_y, plot_max_y)) +
  geom_text_repel(aes(x = mean_ct_long, y = mean_ct_short, label = week), 
                  size = 4, colour = "grey30") +
  geom_abline(intercept = 0, slope = 1, color="grey40", alpha=0.5, linewidth=0.5,
              linetype=2) +
  labs(x = "Mean Ct of Each Time Point (long primer)", 
       y = "Mean Ct of Each Time Point (short primer)", fill = "Week")

# remove # to save plot
ggsave("plots/long_vs_short.pdf", units = "in", width = 10, height = 8)
ggsave("plots/long_vs_short.svg", units = "in", width = 10, height = 8)

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

dry_v_frozen_wide$week <- as.factor(dry_v_frozen_wide$week)

ggplot(dry_v_frozen_wide) +
  geom_point(aes(x = mean_ct_Dry, y = mean_ct_Frozen, fill = week), shape = 21, 
             size = 4, stroke = 0.25, alpha = 0.75) + 
  geom_errorbarh(aes(xmin = mean_ct_Dry - sd_ct_Dry, 
                     xmax = mean_ct_Dry + sd_ct_Dry, y = mean_ct_Frozen), 
                 height = 0.25, color = "slategray3", linewidth = 0.25, alpha = 0.75) +
  geom_errorbar(aes(ymin = mean_ct_Frozen - sd_ct_Frozen, 
                    ymax = mean_ct_Frozen + sd_ct_Frozen, x = mean_ct_Dry), 
                width = 0.25, color = "slategray3", linewidth = 0.25, alpha = 75) +
  facet_wrap(~factor(target, levels = c("Galbut virus", "La Jolla virus", 
                                        "Nora virus", "Thika virus", 
                                        "RpL32 mRNA")), nrow = 2) +
  stat_poly_line(aes(x = mean_ct_Dry, mean_ct_Frozen), method = "lm", alpha = 0.5, 
                 se = FALSE, colour = "slategray3", linetype = "dotdash", 
                 linewidth = 0.5) +
  stat_poly_eq(aes(x = mean_ct_Dry, mean_ct_Frozen, 
                   label = paste(after_stat(eq.label), 
                                 after_stat(rr.label), sep = "*\", \"*"))) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  coord_fixed(xlim = c(plot_min_x_df, plot_max_x_df), 
              ylim = c(plot_min_y_df, plot_max_y_df)) +
  geom_text_repel(aes(x = mean_ct_Dry, y = mean_ct_Frozen, label = week), 
                  size = 4, colour = "grey30") +
  geom_abline(intercept = 0, slope = 1, color="grey40", alpha=0.5, linewidth=0.5,
              linetype=2) +
  labs(x = "Mean Ct of Each Time Point (Dry)", 
       y = "Mean Ct of Each Time Point (Frozen)", fill = "Week")

# remove # to save plot
ggsave("plots/dry_vs_frozen.pdf", units = "in", width = 10, height = 8)
ggsave("plots/dry_vs_frozen.svg", units = "in", width = 10, height = 8)


# Mosquito figures

# Week 4 for normalizing
w4 <- mosquito %>% 
  filter(group == "Frozen",
         present == "y",
         week == 4) %>%
  select(week, target, ct) %>% 
  type_convert() %>% 
  group_by(week, target) %>% 
  summarise(mean_ct_fresh = mean(ct))

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
  mutate(delta_ct = mean_ct_fresh - ct) %>% 
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
ggplot(mos_data3, aes(x = week)) +
  geom_point(aes(y = delta_ct, fill = group), shape = 21, size = 1.35, 
             stroke = 0.1, color = "black", alpha = 0.75) +
  stat_compare_means(aes(y= delta_ct, group = group), label = "p.signif", hide.ns = TRUE, label.y = 6, size = 6.5) +
  geom_line(aes(y = mean_dct, group = group, linetype = group, colour = group), linewidth = 0.75, alpha = 0.75) +
  scale_fill_manual(values = c("firebrick", "navyblue")) +
  scale_colour_manual(values = c("firebrick", "navyblue")) +
  geom_errorbar(aes(ymin = (mean_dct - sd_dct), ymax = (mean_dct + sd_dct)), 
                width = 1, color = "grey50", alpha = 0.4) +  
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  facet_wrap(~factor(target, levels = c("Verdadero virus", "Guadeloupe mosquito virus", 
                                        "Actin mRNA")), ncol = 1) +
  labs(x = 'Weeks After Collection', 
       y = "Log(2) Fold Change Relative to Time Point \n4 Week Frozen Mosquito",
       linetype = "Sample Storage",
       fill = "Sample Storage",
       colour = "Sample Storage") 
  
# remove # to save plot
ggsave("plots/Relative_mosquito_dct.pdf", units = "in", width = 10, height = 8)
ggsave("plots/Relative_mosquito_dct.svg", units = "in", width = 10, height = 8)

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

dry_v_frozen_wide_mos$week <- as.factor(dry_v_frozen_wide_mos$week)

ggplot(dry_v_frozen_wide_mos) +
  geom_point(aes(x = mean_ct_Dry, y = mean_ct_Frozen, fill = week), shape = 21, 
             size = 4, stroke = 0.25, alpha = 0.75) +
  geom_errorbarh(aes(xmin = mean_ct_Dry - sd_ct_Dry, 
                     xmax = mean_ct_Dry + sd_ct_Dry, y = mean_ct_Frozen), 
                 height = 0.25, color = "slategray3", linewidth = 0.25, alpha = 0.75) +
  geom_errorbar(aes(ymin = mean_ct_Frozen - sd_ct_Frozen, 
                    ymax = mean_ct_Frozen + sd_ct_Frozen, 
                    x = mean_ct_Dry), 
                width = 0.25, color = "slategray3", linewidth = 0.25, alpha = 0.75) +
  facet_wrap(~factor(target, levels = c("Verdadero virus", "Guadeloupe mosquito virus", 
                                        "Actin mRNA")), ncol = 2) +
  stat_poly_line(aes(x = mean_ct_Dry, mean_ct_Frozen), method = "lm", alpha = 0.5, 
                 se = FALSE, colour = "slategray3", linetype = "dotdash", linewidth = 0.5) +
  stat_poly_eq(aes(x = mean_ct_Dry, mean_ct_Frozen, 
                   label = paste(after_stat(eq.label), 
                                 after_stat(rr.label), sep = "*\", \"*"))) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  geom_text_repel(aes(x = mean_ct_Dry, y = mean_ct_Frozen, label = week), 
                  size = 4, colour = "grey30") +
  geom_abline(intercept = 0, slope = 1, color="grey40", alpha=0.5, linewidth=0.5,
              linetype=2) +
  labs(x = "Mean Ct of Each Time Point (Dry)", 
       y = "Mean Ct of Each Time Point (Frozen)", fill = "Week")

# remove # to save plot
ggsave("plots/dry_vs_frozen_mos.pdf", units = "in", width = 10, height = 8)
ggsave("plots/dry_vs_frozen_mos.svg", units = "in", width = 10, height = 8)

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

# long vs short galbut & RpL32
ggplot(filter(fly_summary, target %in% c("Galbut virus", "RpL32 mRNA"))) +
  geom_point(aes(x = week, y = n, fill = group), shape = 21, 
              size = 3, stroke = 0.25, alpha = 0.65) +
  scale_fill_manual(values = c("firebrick", "navyblue")) +
  facet_grid(length ~ target) +
  theme_few(base_size = 11) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  labs(x = "Weeks After Collection", y = "Number of Flies Positive", fill = "Sample Storage") 

ggsave("plots/n_long_short.pdf", units = "in", width = 10, height = 8)
ggsave("plots/n_long_short.svg", units = "in", width = 10, height = 8)

# other fly targets
ggplot(filter(fly_summary, target %in% c("Thika virus", "Nora virus", "La Jolla virus"))) +
  geom_point(aes(x = week, y = n, fill = group), shape = 21, 
              size = 3, stroke = 0.25, alpha = 0.75) +
  scale_fill_manual(values = c("firebrick", "navyblue")) +
  facet_grid(group ~ factor(target, levels = c("Nora virus", "La Jolla virus", 
                                               "Thika virus"))) +
  theme_few(base_size = 11) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  labs(x = "Weeks After Collection", y = "Number of Flies Positive", fill = "Sample Storage")

ggsave("plots/n_fly_targets.pdf", units = "in", width = 10, height = 8)
ggsave("plots/n_fly_targets.svg", units = "in", width = 10, height = 8)

# mosquito targets
mos_summary <- mosquito %>% 
  group_by(target, group, week) %>% 
  count(present) %>% 
  mutate(target = str_replace(target, "Actin", "Actin mRNA"),
         target = str_replace(target, "Rennavirus", "Guadeloupe mosquito virus"),
         target = str_replace(target, "Verdadero", "Verdadero virus"))

mos_summary$week <- as.numeric(mos_summary$week)

ggplot(mos_summary) +
  geom_point(aes(x = week, y = n, fill = group), shape = 21, 
              size = 3, stroke = 0.25, alpha = 0.75) +
  scale_fill_manual(values = c("firebrick", "navyblue")) +
  facet_grid(group ~ factor(target, levels = c("Verdadero virus", "Guadeloupe mosquito virus", 
                                               "Actin mRNA"))) +
  theme_few(base_size = 11) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) +
  labs(x = "Weeks After Collection", y = "Number of Mosquitoes Positive", fill = "Sample Storage")

ggsave("plots/n_mos_targets.pdf", units = "in", width = 10, height = 8)
ggsave("plots/n_mos_targets.svg", units = "in", width = 10, height = 8)

# STATS
fly_stats <- fly_fc %>% 
  ungroup() %>% 
  group_by(week, target) %>% 
  t_test(delta_ct ~ group)

mos_stats <- mosquito_data3 %>% 
  ungroup() %>% 
  group_by(week, target) %>% 
  filter(target != "Guadeloupe mosquito virus") %>% 
  t_test(delta_ct ~ group)

mos_stats_gmv <- mosquito_data3 %>% 
  ungroup() %>% 
  group_by(week, target) %>% 
  filter(target == "Guadeloupe mosquito virus",
         week != 32) %>% 
  t_test(delta_ct ~ group)

mos_stats_joined <- full_join(mos_stats, mos_stats_gmv)

