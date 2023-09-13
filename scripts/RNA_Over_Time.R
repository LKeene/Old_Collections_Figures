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

# pre-processing 
#fly_data2 <- fly %>%
#  filter(present == "y") %>% 
#  type_convert() %>% 
#  group_by(target, group, week, length) %>%
#  mutate(mean_ct = mean(ct, na.rm = TRUE),
#         sd_ct = sd(ct, na.rm = TRUE),
#         sample_length = paste(target, length, sep = "_")) 

# determine delta ct, use average of time 2 week fresh as starting ct
#write_xlsx(fly_data2, "tidy_formats/fly_data3.xlsx")

#read back in the data
fly_data3 <- read_xlsx("tidy_formats/fly_data3.xlsx")

fly_data3 <- fly_data3 %>% 
  mutate(target = str_replace(target, "Galbut", "Galbut virus"),
         target = str_replace(target, "RpL32", "RpL32 mRNA"),
         target = str_replace(target, "La Jolla", "La Jolla virus"),
         target = str_replace(target, "Nora", "Nora virus"),
         target = str_replace(target, "Thika", "Thika virus"))


hline_max1 <- data.frame(group = c("Dry", "Frozen"), target = c("La Jolla virus", 
                                                                "Nora virus",
                                                                "Thika virus"), 
                        hline_1 = c(30.81469, 27.84886, 32.70908, 28.27032,
                                    27.49015, 30.87427))


hline_min1 <- data.frame(group = c("Dry", "Frozen"), target = c("La Jolla virus", 
                                                                "Nora virus", 
                                                                "Thika virus"),
                        hline_2 = c(16.96781, 17.54493, 26.47017, 22.34384,
                                    18.78465, 22.54356))


# Mean line bar of La Jolla, Nora & Thika targets
ggplot(filter(fly_data3, target %in% c("La Jolla virus", "Nora virus", 
                                       "Thika virus"), 
              group != "Fresh"), aes(x=week)) + 
  geom_hline(data = hline_max1, aes(yintercept = hline_1), alpha = 0.5) +
  geom_hline(data = hline_min1, aes(yintercept = hline_2), alpha = 0.5) +
  geom_point(aes(y=ct, fill=target), shape = 21, size = 1.5, stroke = 0.1, 
             color = "black", alpha = 0.75) + 
  scale_fill_manual(values = c("turquoise3", "navyblue", "blue")) + 
  geom_line(aes(y=mean_ct, group=sample_length)) + 
  geom_errorbar(aes(ymin = (mean_ct - sd_ct), ymax = (mean_ct + sd_ct)), 
                width = 0.2, color = "grey50", alpha = 0.25) +
  theme_few(base_size = 11) + 
  facet_grid(target ~group, scales = "free_y") +
  labs(x = "Time Since Sample Collection (weeks)", y = "Mean Ct", fill = "Target")


# remove # to save plot
#ggsave("plots/Mean_LaJollaNoraThika.pdf", units = "in", width = 10, height = 8)


hline_max2 <- data.frame(group = c("Dry", "Dry", "Dry", "Dry", "Frozen", "Frozen",
                                   "Frozen", "Frozen"), 
                         sample_length = c("Galbut_long", "Galbut_short",
                                           "RpL32_long", "RpL32_short"), 
                        hline_1 = c(27.11656, 21.71760, 29.948694, 25.24503, 
                                    25.20976, NA, 23.35923, NA))


hline_min2 <- data.frame(group = c("Dry", "Dry", "Dry", "Dry", "Frozen", "Frozen",
                                   "Frozen", "Frozen"), 
                         sample_length = c("Galbut_long", "Galbut_short",
                                           "RpL32_long", "RpL32_short"),
                        hline_2 = c(18.23504, 16.04849, 20.2549, 20.49891, 
                                    17.29496, NA, 17.70815, NA))

# Mean line bar of galbut and RpL32 short and long
ggplot(filter(fly_data3, sample_length %in% c("Galbut_long", "Galbut_short", 
                                             "RpL32_long", "RpL32_short"), 
              group != "Fresh"), aes(x=week)) + 
  geom_hline(data = hline_max2, aes(yintercept = hline_1), alpha = 0.5) +
  geom_hline(data = hline_min2, aes(yintercept = hline_2), alpha = 0.5) +
  geom_point(aes(y=ct, fill=sample_length), shape = 21, size = 1.5, stroke = 0.1, 
             color = "black", alpha = 0.75) + 
  scale_fill_manual(values = c("turquoise3", "navyblue", "blue", "purple")) + 
  geom_line(aes(y=mean_ct, group=sample_length)) + 
  geom_errorbar(aes(ymin = (mean_ct - sd_ct), ymax = (mean_ct + sd_ct)), 
                width = 0.2, color = "grey50", alpha = 0.25) +
  theme_few(base_size = 11) + 
  facet_grid(sample_length ~ group, scales = "free_y") +
  labs(x = "Time Since Sample Collection (weeks)", y = "Mean Ct" )


# remove # to save plot
#ggsave("plots/Mean_GalbutRpL32.pdf", units = "in", width = 10, height = 8)

# get mean fold change, remove fresh & short samples
fly_fc <- fly_data3 %>% 
  group_by(target, group, week) %>% 
  mutate(fold_change = 2^-(delta_ct),
         mean_fc = mean(fold_change, na.rm = TRUE),
         sd_fc = sd(fold_change, na.rm = TRUE),
         mean_dct = mean(delta_ct, na.rm = TRUE),
         sd_dct = sd(delta_ct, na.rm = TRUE)) %>% 
  filter(group != "Fresh") %>% 
  filter(length != "short")


# relative change for all fly targets
ggplot(fly_fc, aes(x = week)) +
  geom_point(aes(y = fold_change, fill = group), shape = 21, size = 1, 
             stroke = 0.1, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("turquoise3", "purple")) +
  geom_line(aes(y = mean_fc, group = group, linetype = group)) +
  scale_y_continuous(trans = "log2") +
  theme_few(base_size = 11) +
  facet_wrap(~factor(target, levels = c("Galbut virus", "La Jolla virus", 
                                        "Nora virus", "Thika virus", 
                                        "RpL32 mRNA")), ncol = 2) +
  labs(x = 'Weeks After Collection', y = "Relative change", fill = "Target", 
       linetype = "Group")

# remove # based on what transformation was used in scale_y_continuous
#ggsave("plots/Relative_fly.pdf", units = "in", width = 10, height = 8)
#ggsave("plots/Relative_fly_log10.pdf", units = "in", width = 10, height = 8)
#ggsave("plots/Relative_fly_log2.pdf", units = "in", width = 10, height = 8)
#ggsave("plots/Relative_fly_sqrt.pdf", units = "in", width = 10, height = 8)

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
        axis.text = element_text(face = "bold")) +
  facet_wrap(~factor(target, levels = c("Galbut virus", "Nora virus",
                                        "RpL32 mRNA")), ncol = 1) +
  labs(x = 'Weeks After Collection', 
       y = "Log(2) Fold Change Relative to Time Point 0 Fresh FoCo-17", 
       fill = "Sample Storage", linetype = "Sample Storage", colour = "Sample Storage")

# remove # to save plot
ggsave("plots/Relative_fly_dct_3targets.pdf", units = "in", width = 10, height = 8)

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
        axis.text = element_text(face = "bold")) +
  facet_wrap(~factor(target, levels = c("La Jolla virus", "Thika virus")), ncol = 1, scales = "free_y") +
  labs(x = 'Weeks After Collection', 
       y = "Log(2) Fold Change Relative to Time Point 0 Fresh FoCo-17", 
       fill = "Sample Storage", linetype = "Sample Storage", colour = "Sample Storage")

ggsave("plots/Relative_fly_dct_SuppTargets.pdf", units = "in", width = 10, height = 8)

# short vs long Galbut & Rpl
short_v_long <- fly_data3 %>% 
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
        axis.text = element_text(face = "bold")) +
  coord_fixed(xlim = c(plot_min_x, plot_max_x), ylim = c(plot_min_y, plot_max_y)) +
  geom_text_repel(aes(x = mean_ct_long, y = mean_ct_short, label = week), 
                  size = 4, colour = "grey30") +
  geom_abline(intercept = 0, slope = 1, color="grey40", alpha=0.5, linewidth=0.5,
              linetype=2) +
  labs(x = "Mean Ct of Each Time Point (long primer)", 
       y = "Mean Ct of Each Time Point (short primer)", fill = "Week")

# remove # to save plot
ggsave("plots/long_vs_short.pdf", units = "in", width = 10, height = 8)

# dry vs frozen fly
plot_min_x_df <- 16
plot_max_x_df <- 33
plot_min_y_df <- 16
plot_max_y_df <- 33

dry_v_frozen <- fly_data3 %>% 
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
        axis.text = element_text(face = "bold")) +
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

# pre-processing mosquito plots
#mosquito_data2 <- mosquito %>% 
#  filter(present == "y") %>% 
#  type_convert() %>% 
#  group_by(target, group, week) %>%
#  mutate(mean_ct = mean(ct, na.rm = TRUE),
#         sd_ct = sd(ct, na.rm = TRUE),
#         sample_group = paste0(target)) %>% 
#  ungroup()

# determine delta ct, use average of time 4 week frozen as starting ct
#write_xlsx(mosquito_data2, "tidy_formats/mosquito_data3.xlsx")

#read back in the data
mosquito_data3 <- read_xlsx("tidy_formats/mosquito_data3.xlsx")


# get mean fold change
mosquito_data3 <- mosquito_data3 %>% 
  group_by(target, group, week) %>% 
  mutate(fold_change = 2^-(delta_ct),
         mean_fc = mean(fold_change, na.rm = TRUE),
         sd_fc = sd(fold_change, na.rm = TRUE),
         mean_dct = mean(delta_ct, na.rm = TRUE),
         sd_dct = sd(delta_ct, na.rm = TRUE),
         target = str_replace(target, "Verdadero", "Verdadero virus"),
         target = str_replace(target, "Rennavirus", "Guadeloupe mosquito virus"),
         target = str_replace(target, "Actin", "Actin mRNA"))

hline_max_mos <- data.frame(group = c("Dry", "Dry", "Dry", "Frozen", "Frozen", 
                                      "Frozen"), 
                            target = c("Verdadero virus", "Rennavirus", 
                                       "Actin mRNA"), 
                            hline_1_mos = c(27.94278, 29.54662, 29.98109,
                                            27.91520, 24.06199, 27.85027))

hline_min_mos <- data.frame(group = c("Dry", "Dry", "Dry", "Frozen", "Frozen", 
                                      "Frozen"), 
                            target = c("Verdadero virus", "Guadeloupe mosquito virus", 
                                       "Actin mRNA"),
                            hline_2_mos = c(22.04536, 19.4869, 24.19849,
                                            22.04304, 14.35967, 23.48702))

# Mean line bar of all Ae. aegypti targets 
ggplot(mosquito_data3, aes(x = week)) + 
  geom_hline(data = hline_max_mos, aes(yintercept = hline_1_mos), alpha = 0.3) +
  geom_hline(data = hline_min_mos, aes(yintercept = hline_2_mos), alpha = 0.3) +
  geom_point(aes(y = ct, fill = target), shape=21, size=1, stroke=0.1, 
             color="black", alpha = 0.5) + 
  scale_fill_manual(values = c("turquoise3", "purple", "green")) + 
  geom_line(aes(y=mean_ct, group=sample_group)) + 
  geom_errorbar(aes(ymin = (mean_ct - sd_ct), ymax = (mean_ct + sd_ct)), 
                width = 0.2, color = "grey50", alpha = 0.25) +
  theme_few(base_size = 11) +
  facet_grid(factor(target, levels = c("Verdadero virus", "Guadeloupe mosquito virus", 
                                       "Actin mRNA")) ~ group, scales = "free_y") +
  labs(x = "Time Since Sample Collection (weeks)", y = "Ct value", 
       fill = "Target")

# remove # to save plot
#ggsave("plots/Mean_Mosquito.pdf", units = "in", width = 10, height = 8)

# Relative change for all Ae. aegypti targets
ggplot(mosquito_data3, aes(x = week)) +
  geom_point(aes(y = fold_change, fill = group), shape = 21, size = 1, 
             stroke = 0.1, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("turquoise3", "purple")) +
  geom_line(aes(y = mean_fc, group = group, linetype = group)) +
  scale_y_continuous(trans = "log2") +
  theme_few(base_size = 11) +
  facet_wrap(~factor(target, levels = c("Verdadero virus", "Guadeloupe mosquito virus", 
                                        "Actin mRNA")), scales = "free_y",
             ncol = 1) +
  labs(x = 'Weeks After Collection', y = "Relative change", fill = "Target")

# remove # based on what transformation was used in scale_y_continuous
#ggsave("plots/Relative_Mosquito_log10.pdf", units = "in", width = 10, height = 8)
#ggsave("plots/Relative_Mosquito_log2.pdf", units = "in", width = 10, height = 8)
#ggsave("plots/Relative_Mosquito_sqrt.pdf", units = "in", width = 10, height = 8)
#ggsave("plots/Relative_Mosquito.pdf", units = "in", width = 10, height = 8)

# delta ct fold change
ggplot(mosquito_data3, aes(x = week)) +
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
        axis.text = element_text(face = "bold")) +
  facet_wrap(~factor(target, levels = c("Verdadero virus", "Guadeloupe mosquito virus", 
                                        "Actin mRNA")), ncol = 1) +
  labs(x = 'Weeks After Collection', 
       y = "Log(2) Fold Change Relative to Time Point 4 Week Frozen Mosquito",
       linetype = "Sample Storage",
       fill = "Sample Storage",
       colour = "Sample Storage") 
  
# remove # to save plot
ggsave("plots/Relative_mosquito_dct.pdf", units = "in", width = 10, height = 8)

# dry vs frozen mos
plot_min_x_df_m <- 14
plot_max_x_df_m <- 30
plot_min_y_df_m <- 14
plot_max_y_df_m <- 30

dry_v_frozen_mos <- mosquito_data3 %>%
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
                                        "Actin mRNA")), ncol = 1) +
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
        axis.text = element_text(face = "bold")) +
  geom_text_repel(aes(x = mean_ct_Dry, y = mean_ct_Frozen, label = week), 
                  size = 4, colour = "grey30") +
  geom_abline(intercept = 0, slope = 1, color="grey40", alpha=0.5, linewidth=0.5,
              linetype=2) +
  labs(x = "Mean Ct of Each Time Point (Dry)", 
       y = "Mean Ct of Each Time Point (Frozen)", fill = "Week")

# remove # to save plot
ggsave("plots/dry_vs_frozen_mos.pdf", units = "in", width = 10, height = 8)

# Percent Positive
fly_summary <- fly %>% 
  group_by(target, group, week, length) %>% 
  count(present) %>% 
  mutate(percent = (n/6) * 100) %>% 
  filter(present == "y") %>% 
  filter(group != "Fresh") 

fly_summary$week <- as.numeric(fly_summary$week)

# long vs short galbut & RpL32
ggplot(filter(fly_summary, target %in% c("Galbut", "RpL32"))) +
  geom_jitter(aes(x = week, y = percent, fill = group), shape = 21, 
              size = 3, stroke = 0.25, alpha = 0.65, width = 0.5, height = 0.5) +
  scale_fill_manual(values = c("firebrick", "navyblue")) +
  facet_grid(length ~ target) +
  theme_few(base_size = 11) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold")) +
  labs(x = "Weeks After Collection", y = "Percent of Flies Positive", fill = "Sample Storage") 

ggsave("plots/percent_long_short.pdf", units = "in", width = 10, height = 8)

# other fly targets
ggplot(filter(fly_summary, target %in% c("Thika", "Nora", "La Jolla"))) +
  geom_jitter(aes(x = week, y = percent, fill = group), shape = 21, 
              size = 3, stroke = 0.25, alpha = 0.75, width = 0.5, height = 0.5) +
  scale_fill_manual(values = c("firebrick", "navyblue")) +
  facet_grid(group ~ target) +
  theme_few(base_size = 11) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold")) +
  labs(x = "Weeks After Collection", y = "Percent of Flies Positive", fill = "Sample Storage")

ggsave("plots/percent_fly_targets.pdf", units = "in", width = 10, height = 8)

# mosquito targets
mos_summary <- mosquito %>% 
  group_by(target, group, week) %>% 
  count(present) %>% 
  mutate(percent = (n/6) * 100) %>% 
  filter(present == "y") %>% 
  mutate(target = str_replace(target, "Actin", "Actin mRNA"),
         target = str_replace(target, "Rennavirus", "Gaudeloupe mosquito virus"))

mos_summary$week <- as.numeric(mos_summary$week)

ggplot(mos_summary) +
  geom_jitter(aes(x = week, y = percent, fill = group), shape = 21, 
              size = 3, stroke = 0.25, alpha = 0.75, width = 0.5, height = 0.5) +
  scale_fill_manual(values = c("firebrick", "navyblue")) +
  facet_grid(group ~ target) +
  theme_few(base_size = 11) +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold")) +
  labs(x = "Weeks After Collection", y = "Percent of Flies Positive", fill = "Sample Storage")

ggsave("plots/percent_mos_targets.pdf", units = "in", width = 10, height = 8)

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

