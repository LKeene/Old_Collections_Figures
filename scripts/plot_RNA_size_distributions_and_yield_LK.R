# This script takes raw Tapestation output files to generate figures depicting 
# the average length of RNA obtained from experimentally dried and frozen specimens
# as well as museum specimens.

library(bioanalyzeR)
library(tidyverse)
library(readxl)
library(patchwork)
library(rstatix)
library(performance)
library(gt)
library(webshot2)

# Load all dried fly data
dry_fly_rna1 <- read.electrophoresis("tapestation/2022-08-05 - 11-17-40-HSRNA.xml")
dry_fly_rna2 <- read.electrophoresis("tapestation/2022-08-05 - 11-49-24-HSRNA.xml")
dry_fly_rna3 <- read.electrophoresis("tapestation/2022-10-23 - 11-08-02-HSRNA.xml")
dry_fly_rna4 <- read.electrophoresis("tapestation/2023-03-10 - 13-39-28-HSRNA.xml")
dry_fly_rna5 <- read.electrophoresis("tapestation/2023-03-10 - 15-24-10-HSRNA.xml")

# Load all frozen fly data
frozen_fly_rna1 <- read.electrophoresis("tapestation/2022-10-11 - 10-04-51-HSRNA.xml")
frozen_fly_rna2 <- read.electrophoresis("tapestation/2022-08-04 - 11-37-28-HSRNA.xml")
frozen_fly_rna3 <- read.electrophoresis("tapestation/2023-03-10 - 14-32-39-HSRNA.xml")
frozen_fly_rna4 <- read.electrophoresis("tapestation/2023-03-10 - 15-24-10-HSRNA.xml")

# Load fresh fly data
rna_fresh <- read.electrophoresis("tapestation/2022-08-25 - 10-03-38[14518]-HSRNA.xml")

# Load all dried mosquito data
dry_mos_rna1 <- read.electrophoresis("tapestation/2022-10-07 - 12-36-06-HSRNA.xml")
dry_mos_rna2 <- read.electrophoresis("tapestation/2022-10-10 - 15-01-37-HSRNA.xml")
dry_mos_rna3 <- read.electrophoresis("tapestation/2022-10-10 - 15-12-34-HSRNA.xml")
dry_mos_rna4 <- read.electrophoresis("tapestation/2023-03-13 - 14-23-40-HSRNA.xml")
dry_mos_rna5 <- read.electrophoresis("tapestation/2023-03-24 - 13-19-06-HSRNA.xml")
dry_mos_rna6 <- read.electrophoresis("tapestation/2023-03-24 - 14-03-01-HSRNA.xml")

# Load all frozen mosquito data
frozen_mos_rna1 <- read.electrophoresis("tapestation/2023-02-02 - 12-27-58-HSRNA.xml")
frozen_mos_rna2 <- read.electrophoresis("tapestation/2023-03-13 - 15-05-19-HSRNA.xml")
frozen_mos_rna3 <- read.electrophoresis("tapestation/2023-03-24 - 13-19-06-HSRNA.xml")
frozen_mos_rna4 <- read.electrophoresis("tapestation/2023-03-24 - 14-03-01-HSRNA.xml")
frozen_mos_rna5 <- read.electrophoresis("tapestation/2023-03-24 - 14-35-30-HSRNA.xml")

# Load fresh mosquito data
fresh_mos <- read.electrophoresis("tapestation/2024-01-24 - 15-14-22-HSRNA.xml")

# Load Old Collections data
oc_samples <- read.electrophoresis("tapestation/2022-08-04 - 11-37-28-HSRNA.xml")

# read in sample metadata for each tape
metadata_fly <- read_excel("tapestation/tape_metadata_tidy.xlsx")
metadata_mos <- read_excel("tapestation/mos_metadata_tidy.xlsx")

# pull out the actual data
df1 <- dry_fly_rna1$data
df2 <- dry_fly_rna2$data
df3 <- dry_fly_rna3$data 
df4 <- dry_fly_rna4$data
df5 <- dry_fly_rna5$data

df6 <- frozen_fly_rna1$data
df7 <- frozen_fly_rna2$data
df8 <- frozen_fly_rna3$data
df9 <- frozen_fly_rna4$data

df10 <- rna_fresh$data

df11 <- dry_mos_rna1$data
df12 <- dry_mos_rna2$data
df13 <- dry_mos_rna3$data
df14 <- dry_mos_rna4$data
df15 <- dry_mos_rna5$data
df16 <- dry_mos_rna6$data

df17 <- frozen_mos_rna1$data
df18 <- frozen_mos_rna2$data
df19 <- frozen_mos_rna3$data
df20 <- frozen_mos_rna4$data
df21 <- frozen_mos_rna5$data

df22 <- fresh_mos$data

df_oc <- oc_samples$data

# these values match tape_id column in metadata excel
df1$tape_id <- "f1"
df2$tape_id <- "f2"
df3$tape_id <- "f3"
df4$tape_id <- "f6"
df5$tape_id <- "f8"

df6$tape_id <- "f4"
df7$tape_id <- "f5"
df8$tape_id <- "f7"
df9$tape_id <- "f8"

df10$tape_id <- "fresh"

df11$tape_id <- "m1"
df12$tape_id <- "m2"
df13$tape_id <- "m3"
df14$tape_id <- "m5"
df15$tape_id <- "m7"
df16$tape_id <- "m8"

df17$tape_id <- "m4"
df18$tape_id <- "m6"
df19$tape_id <- "m7"
df20$tape_id <- "m8"
df21$tape_id <- "m9"
df22$tape_id <- "m10"

df_oc$tape_id <- "oc"

# merge the data from different tapes
fly_data <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df_oc) 

mos_data <- rbind(df11, df12, df13, df14, df15, df16, df17, df18, df19, df20, df21, df22)

# make a joint ID that combines tape_id and sample.index 
# this will uniquely identify each lane and allow us to link
# tape data to metadata
fly_data <- fly_data %>% mutate(id = str_c(tape_id, sample.index))
metadata_fly <- metadata_fly %>% mutate(id = str_c(tape_id, sample.index))
mos_data <- mos_data %>% mutate(id = str_c(tape_id, sample.index))
metadata_mos <- metadata_mos %>% mutate(id = str_c(tape_id, sample.index))

# get rid of lanes for which no metadata provided
fly_data <- filter(fly_data, id %in% metadata_fly$id)
mos_data <- filter(mos_data, id %in% metadata_mos$id)

# check we have all the tapes/samples represented as expected
fly_data %>% 
  group_by(tape_id,sample.index) %>% 
  summarize()

# determine cutoff for where the lower marker ends - then check visually
lower_marker_length_max <- 33
                      
#ggplot(fly_data) +       
#  geom_line(aes(x=length, y=fluorescence, group=sample.index)) +
  # add in coloring under the lines: 
  # length < lower_marker_mlength_max will produce a T/F value, which we can color with scale_fill_manual
#  geom_area(aes(x=length, y=fluorescence, fill = length <= lower_marker_length_max)) +
#  scale_fill_manual(values=c("grey50", "lightsteelblue")) +
  # get rid of the ugly legend
#  theme_classic(base_size=10) +
#  scale_x_log10(limits=c(NA,1000)) +
#  xlab("RNA Length (nt)") +
#  ylab("Fluorescence (arbitrary units)") +
#  facet_wrap(~id, ncol=1) +
#  theme(legend.position="none",
#        strip.text.x = element_blank(),
#        axis.ticks.y = element_blank(),
#        axis.text.y  = element_blank()) 


# calculate mean from what is essentially histogram data.
# follow the strategy described here:
# https://www.statology.org/histogram-mean-median/

# get rid of rows with no assigned length or below lower marker cutoff
fly_data <- fly_data %>% filter(!is.na(length) & length > lower_marker_length_max)

# only keep fluorescence values if above an arbitrary threshold 
# this is because legit fluorescence signals are much larger
# and we want to avoid fluorescence signal noise from impacting mean, which it was doing
min_fluorescence_cutoff <- 10
fly_data <- fly_data %>% filter(fluorescence > min_fluorescence_cutoff)

# calculate fluorescence * length for each bin
fly_data <- fly_data %>% mutate(mn = fluorescence * length)

# calculate mean length from each TS trace data.
# and put in a new dataframe with just one value per sample
df_mean_fly <- fly_data %>% 
  group_by(id) %>% 
  summarize(mean_length = sum(mn) / sum(fluorescence))

# merge in metadata
df_mean_fly <- left_join(df_mean_fly, metadata_fly, by="id")
#write_csv(df_mean_fly, "metadata/fly_length_data.csv")

# investigate filtering impact
#ggplot(filter(df, tape_id=="oc" & sample.index < 8)) +
#  facet_wrap(~id, scales="free_y", ncol=1) +
#  scale_x_log10() + 
#  geom_point(aes(x=length, y=fluorescence)) 


# plot oc samples: we'll show these separately from new/experimental samples
p_oc <- ggplot(filter(df_mean_fly, group == "old")) +
  # geom_smooth(aes(x=year, y=mean_length), alpha=0.25, color="slateblue", size=0.5) +
  geom_point(aes(x = year, y = mean_length), shape = 21, size = 8, 
             fill = "orchid4", color = "black", stroke = 0.1, alpha = 0.95) +
  scale_x_reverse(breaks = seq(from=1880, to=2020, by=20))+
  ylim(c(0,500)) +
  xlab("Year of sample collection") +
  ylab("Mean length of RNA (nt)") +
  theme_minimal(base_size = 16) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        #        axis.title = element_text(face = "bold"),
        #        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) 

p_oc
ggsave("plots/Figure2/OC_RNA_length_vs_time.pdf", width=10, height=7, units="in")
ggsave("plots/Figure2/OC_RNA_length_vs_time.svg", width=10, height=7, units="in")
ggsave("plots/Figure2/OC_RNA_length_vs_time.jpg", width=10, height=7, units="in")

# calculate average lengths
oc_only <- df_mean_fly %>% 
  filter(tape_id == "oc") %>% 
  select(mean_length, year)

oc_1960_up <- oc_only %>% 
  filter(year > 1960) %>% 
  summarise(mean(mean_length))

oc_1960_below <- oc_only %>% 
  filter(year < 1960) %>% 
  summarise(mean(mean_length))

# get mean and sd of triplicates
df_mean_fly_grouped <- df_mean_fly %>% 
  filter(group != "old") %>% 
  group_by(group, weeks) %>% 
  mutate(mean_length_g = mean(mean_length),
         sd_mean_length_g = sd(mean_length))

# plot new experimental samples
p_new <- ggplot(df_mean_fly_grouped, aes(as.numeric(weeks))) +
  geom_point(aes(y = mean_length_g, fill = group), show.legend = FALSE,
             shape = 21, size = 6, color = "black", stroke = 0.1, alpha = 0.75) +
# stat_compare_means(aes(y= mean_length_g, group = ), label = "p.signif", 
#                    hide.ns = TRUE, label.y = 7, size = 3.5, alpha = 0.75, symnum.args = 
#                      list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
#                           symbols = c("D", "C", "B", "A", "ns"))) +
  scale_fill_manual(values = c("firebrick3", "gray30", "navyblue")) +
  geom_errorbar(aes(ymin = (mean_length_g - sd_mean_length_g), 
                    ymax = (mean_length_g + sd_mean_length_g), color = group), 
                width = 1.5, alpha = 0.25, show.legend = FALSE) +
  ylim(c(0,1000)) +
  labs(x = "Weeks After Sample Storage", 
       y = "Mean length of RNA (nt)", fill = "Group") +
  theme_minimal(base_size = 16) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        #        axis.title = element_text(face = "bold"),
        #        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) 

p_new
ggsave("plots/Figure1/Fly_RNA_length_vs_time.pdf", width=10, height=7, units="in")
ggsave("plots/Figure1/Fly_RNA_length_vs_time.svg", width=10, height=7, units="in")
ggsave("plots/Figure1/Fly_RNA_length_vs_time.jpg", width=10, height=7, units="in")

# plot the two plots side by side
# lay out combined plot using patchwork library
p_new + p_oc +  plot_layout(widths = c(1.5, 1.5))


ggsave("plots/RNA_length_vs_time.pdf", width=10, height=7, units="in")

# plot RNA length vs time of mosquito samples

# get rid of rows with no assigned length or below lower marker cutoff
mos_data <- mos_data %>% filter(!is.na(length) & length > lower_marker_length_max)

# only keep fluorescence values if above an arbitrary threshold 
# this is because legit fluorescence signals are much larger
# and we want to avoid fluorescence signal noise from impacting mean, which it was doing
mos_data <- mos_data %>% filter(fluorescence > min_fluorescence_cutoff)

# calculate fluorescence * length for each bin
mos_data <- mos_data %>% mutate(mn = fluorescence * length)

# calculate mean length from each TS trace data.
# and put in a new dataframe with just one value per sample
df_mean_mos <- mos_data %>% 
  group_by(id) %>% 
  summarize(mean_length = sum(mn) / sum(fluorescence))

# merge in metadata
df_mean_mos <- left_join(df_mean_mos, metadata_mos, by="id")
#write_csv(df_mean_mos, "metadata/mos_length_data.csv")

# get mean and sd of triplicates
df_mean_mos_grouped <- df_mean_mos %>% 
  group_by(group, weeks) %>% 
  mutate(mean_length_g = mean(mean_length),
         sd_mean_length_g = sd(mean_length))

p_mos <- ggplot(df_mean_mos_grouped, aes(x = as.numeric(weeks))) +
  # geom_smooth(aes(x=as.numeric(weeks), y=mean_length), alpha=0.25, color="slateblue", size=0.5) +
  geom_point(aes(y = mean_length_g, fill = group), 
             shape = 21, size = 6, color = "black", stroke = 0.5, alpha = 0.75) +
  scale_fill_manual(values = c("firebrick", "gray30", "navyblue")) +
  geom_errorbar(aes(ymin = (mean_length_g - sd_mean_length_g), 
                    ymax = (mean_length_g + sd_mean_length_g), color = group), 
                width = 1.5, alpha = 0.25, show.legend = FALSE) +
  ylim(c(0,1400)) +
  labs(x = "weeks After Sample Storage", 
       y = "Mean length of RNA (nt)", fill = "Group") +
  theme_minimal(base_size = 16) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        #        axis.title = element_text(face = "bold"),
        #        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        text = element_text(size = 20)) 

p_mos
ggsave("plots/Supplemental2/Mos_RNA_length_vs_time.pdf", width=10, height=7, units="in")
ggsave("plots/Supplemental2//Mos_RNA_length_vs_time.svg", width=10, height=7, units="in")
ggsave("plots/Supplemental2//Mos_RNA_length_vs_time.jpg", width=10, height=7, units="in")


