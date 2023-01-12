library(tidyverse)
library(readxl)
library(writexl)
library(ggthemes)

df <- read_xlsx("CombinedAllTargets_ForFigure_LK_20220815.xlsx")

df <- rename(.data = df,
             sample_name = Sample,
             sex = Sex,
             target = Target,
             ct = Ct,
             group = Group,
             present = Present,
             weeks = Weeks) %>% 
  filter(present == "Y")

write_xlsx(df, "Mean_AllTargets2.xlsx")

df <- df %>%
  group_by(target, group, weeks) %>%
  mutate(mean_ct = mean(ct, na.rm = TRUE),
         sd_ct = sd(ct, na.rm = TRUE),
         sample_group = paste0(target)) 


hline_max <- data.frame(group = c("Dry", "Frozen"), target = c("Galbut",
                                                               "La Jolla", "Nora", 
                                                               "RpL32", "Thika"), 
                        hline_1 = c(27.1165, 28.68144, 27.49, 23.359, 32.70908, 
                                    22.546, 29.2265, 27.84886, 27.414, 29.43864))


hline_min <- data.frame(group = c("Dry", "Frozen"), target = c("Galbut",
                                                               "La Jolla", "Nora",
                                                               "RpL32", "Thika"),
                        hline_2 = c(18.23504, 22.3484, 18.972817, 17.70815, 
                                    26.47017, 16.07851, 16.96781, 17.54493, 
                                    20.2565, 24.26979))


# Mean line bar  
ggplot(filter(df, group != "Fresh"), aes(x=weeks)) + 
  geom_hline(data = hline_max, aes(yintercept = hline_1), alpha = 0.5) +
  geom_hline(data = hline_min, aes(yintercept = hline_2), alpha = 0.5) +
  geom_point(aes(y=mean_ct, fill=target), shape=21, size=2, stroke=0.1, color="black") + 
  scale_fill_manual(values = c("turquoise3", "navyblue", "orange", "yellow", "green")) + 
  geom_line(aes(y=mean_ct, group=sample_group)) + 
  geom_errorbar(aes(ymin = (mean_ct - sd_ct), ymax = (mean_ct + sd_ct)), width = 0.2, color = "grey50", alpha = 0.25) +
  theme_few(base_size = 11) + 
  facet_grid(factor(target, levels = c("Galbut", "La Jolla", "Nora", "Thika", "RpL32")) ~group, scales = "free_y") +
  labs(x = "Time Since Sample Collection (weeks)", y = "Mean Ct" )



ggsave("Mean_AllTargetsLK_few.pdf", units = "in", width = 10, height = 8)


# df2 <- read_xlsx("Mean_AllTargets.xlsx")
# 
# df2 <- mutate(.data = df2, 
#               relative = (2^-delta_ct))
# 
# df2 <- df2 %>%
#   group_by(target, group, weeks) %>%
#   mutate(mean_relative = mean(relative, na.rm = TRUE),
#          sd_relative = sd(relative, na.rm = TRUE))
# 
# # Relative RNA line graph
# ggplot(filter(df2, group != "Fresh"), aes(x=weeks)) + 
#   geom_point(aes(y=mean_relative, fill=target), shape=21, size=2, stroke=0.1, color="black") + 
#   scale_fill_manual(values = c("turquoise3", "navyblue", "orange", "yellow", "green")) + 
#   geom_line(aes(y=mean_relative, group=sample_group)) + 
#   geom_errorbar(aes(ymin = (mean_relative - sd_relative), ymax = (mean_relative + sd_relative)), width = 0.2, color = "grey50", alpha = 0.5) +
#   theme_bw(base_size = 11) + 
#   facet_grid(target~group, scales = "free_y") +
#   scale_y_continuous(trans = "log10") +
#   labs(x = "Time Since Sample Collection (weeks)", y = "Relative Ct to Time 0 Samples" )
# 
# ggsave("Relative_AllTargets_log10_LK.pdf")

# Mosquito Data
df2 <- read_xlsx("Combined_PozaRica_Tidy_LK.xlsx")

# Group by with sex
df_sex <- df2 %>% 
  rename(sample_name = Sample,
         sex = Sex,
         target = Target,
         ct = Ct,
         group = Group,
         present = Present,
         weeks = Weeks) %>% 
  filter(present == "Y") %>% 
  group_by(target, group, weeks, sex) %>%
  mutate(mean_ct = mean(ct, na.rm = TRUE),
         sd_ct = sd(ct, na.rm = TRUE),
         sample_group = paste0(target)) 


# Mean line bar  
ggplot(df_sex, aes(x=weeks)) + 
  geom_point(aes(y=mean_ct, fill = target), shape=21, size=2, stroke=0.1, color="black") + 
  scale_fill_manual(values = c("turquoise3", "purple", "green")) + 
  geom_line(aes(y=mean_ct, group=sample_group)) + 
  geom_errorbar(aes(ymin = (mean_ct - sd_ct), ymax = (mean_ct + sd_ct)), width = 0.2, color = "grey50", alpha = 0.25) +
  theme_few(base_size = 11) +
  facet_grid(sex~group, scales = "free_y") +
  labs(x = "Time Since Sample Collection (weeks)", y = "Mean Ct" )

# Group by without sex
df_nosex <- df2 %>% 
  rename(sample_name = Sample,
         sex = Sex,
         target = Target,
         ct = Ct,
         group = Group,
         present = Present,
         weeks = Weeks) %>% 
  filter(present == "Y") %>% 
  group_by(target, group, weeks) %>%
  mutate(mean_ct = mean(ct, na.rm = TRUE),
         sd_ct = sd(ct, na.rm = TRUE),
         sample_group = paste0(target)) 

write_xlsx(df_nosex, "Mean_Mosquito.xlsx")

hline_max_mos <- data.frame(group = c("Dry", "Dry", "Dry", "Frozen", "Frozen", 
                                      "Frozen"), 
                            target = c("Actin", "Verdadero", "Rennavirus"), 
                         hline_1_mos = c(29.98109404, 27.94278399, 27.48395252, 
                                         27.48395252, 27.91520468, 24.062))

hline_min_mos <- data.frame(group = c("Dry", "Dry", "Dry", "Frozen", "Frozen", 
                                      "Frozen"), 
                            target = c("Verdadero", "Actin", "Rennavirus"),
                         hline_2_mos = c(22.33625, 24.19849205, 19.487, 
                                         22.04304091, 23.48701763, 14.35966667))


# Mean line bar  
ggplot(df_nosex, aes(x=weeks)) + 
  geom_hline(data = hline_max_mos, aes(yintercept = hline_1_mos), alpha = 0.3) +
  geom_hline(data = hline_min_mos, aes(yintercept = hline_2_mos), alpha = 0.3) +
  geom_point(aes(y=mean_ct, fill = target), shape=21, size=2, stroke=0.1, 
             color="black") + 
  scale_fill_manual(values = c("cyan2", "purple", "green")) + 
  geom_line(aes(y=mean_ct, group=sample_group)) + 
  geom_errorbar(aes(ymin = (mean_ct - sd_ct), ymax = (mean_ct + sd_ct)),
                width = 0.1, color = "grey50", alpha = 0.1) +
  theme_few(base_size = 11) +
  facet_grid(factor(target, levels = c("Verdadero", "Rennavirus", "Actin")) 
             ~group, scales = "free_y") +
  labs(x = "Time Since Sample Collection (weeks)", y = "Mean Ct", fill = "Target" )



ggsave("Mean_Mosquito.pdf", units = "in", width = 10, height = 8)


