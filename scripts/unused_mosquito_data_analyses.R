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

# Calc. ct differences
changes2 <- mos_data3 %>% 
  filter(week %in% c(4, 52),
         sex == "M",
         rep == 1) %>% 
  select(week, group, target, mean_dct)

changes2_wide <- changes2 %>% 
  pivot_wider(names_from = week, values_from = mean_dct) %>% 
  rename(`4_week` = `4`,
         `52_week` = `52`) %>% 
  mutate(diff = `4_week` - `52_week`)

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


# MLRs
mos_length <- read_csv("metadata/mos_length_data.csv")

# Mosquito concentration summary table
conc %>% 
  filter(organism == "mosquito") %>% 
  select(week, storage, concentration) %>% 
  group_by(week, storage) %>% 
  summarize(n = n(),
            Mean = mean(concentration),
            SD = sd(concentration),
            SE_mean = SD/sqrt(n))

mos_length <- mos_length %>% 
  select(weeks, group, mean_length, sample_name) %>% 
  group_by(weeks, group) %>% 
  mutate(group = as.factor(group))

# Mosquito length summary table
mos_length %>% 
  select(weeks, group, mean_length) %>% 
  group_by(weeks, group) %>% 
  summarize(n = n(),
            Mean = mean(mean_length),
            SD = sd(mean_length),
            SE_mean = SD/sqrt(n)) 

df_mean_mos <- mos_length %>% 
  group_by(group, weeks) %>% 
  mutate(mean_length_g = mean(mean_length),
         sd_mean_length_g = sd(mean_length))

## Mosquito Concentration
mos_conc <- conc_mean %>% 
  filter(organism == "mosquito")

# Mosquito concentration One-Way Model
conc_mos_mlr <- lm(concentration~week+storage, data = mos_conc)
tidy(conc_mos_mlr)
Anova(conc_mos_mlr)

check_model(conc_mos_mlr, check = c("qq", "linearity", "homogeneity"))

emm_mos_conc <- emmeans(conc_mos_mlr, ~week)
pwpp(emm_mos_conc)


# Mosquito length One-Way Model
df_mean_mos <- df_mean_mos %>% 
  filter(group != "Fresh")
length_mos_mlr <- lm(mean_length_g~weeks+group, data = df_mean_mos)
tidy(length_mos_mlr)
Anova(length_mos_mlr)

check_model(length_mos_mlr, check = c("qq", "linearity", "homogeneity"))

emm_mos_leng <- emmeans(length_mos_mlr, ~weeks)
pwpp(emm_mos_leng)

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
