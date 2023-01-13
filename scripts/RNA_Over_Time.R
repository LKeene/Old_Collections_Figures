# load libraries
library(tidyverse)
library(readxl)
library(ggthemes)
library(writexl)

# read in data
fly <- read_xlsx("tidy_formats/fly_data2.xlsx")

mosquito <- read_xlsx("tidy_formats/mosquito_data2.xlsx")

# fly data figures

# pre-processing 
fly_data2 <- fly %>%
  filter(present == "y") %>% 
  type_convert() %>% 
  group_by(target, group, week, length) %>%
  mutate(mean_ct = mean(ct, na.rm = TRUE),
         sd_ct = sd(ct, na.rm = TRUE),
         sample_length = paste(target, length, sep = "_")) 

# determine delta ct, use average of time 2 week fresh as starting ct
write_xlsx(fly_data2, "tidy_formats/fly_data3.xlsx")

#read back in the data
fly_data3 <- read_xlsx("tidy_formats/fly_data3.xlsx")

hline_max1 <- data.frame(group = c("Dry", "Frozen"), target = c("La Jolla", "Nora", 
                                                               "Thika"), 
                        hline_1 = c(30.81469, 27.84886, 32.70908, 28.27032,
                                    27.49015, 30.87427))


hline_min1 <- data.frame(group = c("Dry", "Frozen"), target = c("La Jolla", "Nora",
                                                               "Thika"),
                        hline_2 = c(16.96781, 17.54493, 26.47017, 22.34384,
                                    18.78465, 22.54356))


# Mean line bar of La Jolla, Nora & Thika targets
ggplot(filter(fly_data3, target %in% c("La Jolla", "Nora", "Thika"), 
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

# get mean fold change, remove fresh samples & short samples
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
  facet_wrap(~factor(target, levels = c("Galbut", "La Jolla", "Nora", "Thika", "RpL32")),
             ncol = 2) +
  labs(x = 'Weeks After Collection', y = "Relative change", fill = "Target")

# remove # based on what transformation was used in scale_y_continuous
#ggsave("plots/Relative_fly.pdf", units = "in", width = 10, height = 8)
#ggsave("plots/Relative_fly_log10.pdf", units = "in", width = 10, height = 8)
#ggsave("plots/Relative_fly_log2.pdf", units = "in", width = 10, height = 8)
#ggsave("plots/Relative_fly_sqrt.pdf", units = "in", width = 10, height = 8)

# delta ct fold change
ggplot(fly_fc, aes(x = week)) +
  geom_point(aes(y = delta_ct, fill = group), shape = 21, size = 1, 
             stroke = 0.1, color = "black", alpha = 0.35) +
  scale_fill_manual(values = c("turquoise3", "purple")) +
  geom_line(aes(y = mean_dct, group = group, linetype = group)) +
  geom_errorbar(aes(ymin = (mean_dct - sd_dct), ymax = (mean_dct + sd_dct)), 
                width = 1, color = "grey50", alpha = 0.25) +  
  theme_few(base_size = 11) +
  facet_wrap(~factor(target, levels = c("Galbut", "La Jolla", "Nora", "Thika", 
                                        "RpL32")), ncol = 2) +
  labs(x = 'Weeks After Collection', 
       y = "Fold Change Relative to Time Point 0 Fresh FoCo-17", 
       fill = "Target", linetype = "Group")

# remove # to save plot
#ggsave("plots/Relative_fly_dct.pdf", units = "in", width = 10, height = 8)

# pre-processing 
mosquito_data2 <- mosquito %>% 
  filter(present == "y") %>% 
  type_convert() %>% 
  group_by(target, group, week) %>%
  mutate(mean_ct = mean(ct, na.rm = TRUE),
         sd_ct = sd(ct, na.rm = TRUE),
         sample_group = paste0(target)) %>% 
  ungroup()

# determine delta ct, use average of time 4 week frozen as starting ct
write_xlsx(mosquito_data2, "tidy_formats/mosquito_data3.xlsx")

#read back in the data
mosquito_data3 <- read_xlsx("tidy_formats/mosquito_data3.xlsx")

# get mean fold change
mosquito_data3 <- mosquito_data3 %>% 
  group_by(target, group, week) %>% 
  mutate(fold_change = 2^-(delta_ct),
         mean_fc = mean(fold_change, na.rm = TRUE),
         sd_fc = sd(fold_change, na.rm = TRUE),
         mean_dct = mean(delta_ct, na.rm = TRUE),
         sd_dct = sd(delta_ct, na.rm = TRUE))

hline_max_mos <- data.frame(group = c("Dry", "Dry", "Dry", "Frozen", "Frozen", 
                                      "Frozen"), 
                            target = c("Verdadero", "Rennavirus", "Actin"), 
                            hline_1_mos = c(27.94278, 29.54662, 29.98109,
                                            27.91520, 24.06199, 27.85027))

hline_min_mos <- data.frame(group = c("Dry", "Dry", "Dry", "Frozen", "Frozen", 
                                      "Frozen"), 
                            target = c("Verdadero", "Rennavirus", "Actin"),
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
  facet_grid(factor(target, levels = c("Verdadero", "Rennavirus", "Actin"))
             ~ group, scales = "free_y") +
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
  facet_wrap(~factor(target, levels = c("Verdadero", "Rennavirus", "Actin")), 
             scales = "free_y", ncol = 1) +
  labs(x = 'Weeks After Collection', y = "Relative change", fill = "Target")

# remove # based on what transformation was used in scale_y_continuous
#ggsave("plots/Relative_Mosquito_log10.pdf", units = "in", width = 10, height = 8)
#ggsave("plots/Relative_Mosquito_log2.pdf", units = "in", width = 10, height = 8)
#ggsave("plots/Relative_Mosquito_sqrt.pdf", units = "in", width = 10, height = 8)
#ggsave("plots/Relative_Mosquito.pdf", units = "in", width = 10, height = 8)

# delta ct fold change
ggplot(mosquito_data3, aes(x = week)) +
  geom_point(aes(y = delta_ct, fill = group), shape = 21, size = 1, 
             stroke = 0.1, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("turquoise3", "purple")) +
  geom_line(aes(y = mean_dct, group = group, linetype = group)) +
  geom_errorbar(aes(ymin = (mean_dct - sd_dct), ymax = (mean_dct + sd_dct)), 
                width = 1, color = "grey50", alpha = 0.25) + 
  theme_few(base_size = 11) +
  facet_wrap(~factor(target, levels = c("Verdadero", "Rennavirus", "Actin")),
             ncol = 1) +
  labs(x = 'Weeks After Collection', 
       y = "Fold Change Relative to Time Point 4 Week Mosquito", 
       fill = "Target", linetype = "Group")

# remove # to save plot
ggsave("plots/Relative_mosquito_dct.pdf", units = "in", width = 10, height = 8)
