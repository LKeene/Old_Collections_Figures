library(bioanalyzeR)
library(tidyverse)
library(svglite)

rna_time_1 <- read.electrophoresis("tapestation/2022-08-05 - 11-17-40-HSRNA.xml")
rna_time_2 <- read.electrophoresis("tapestation/2022-08-05 - 11-49-24-HSRNA.xml")
rna_time_3 <- read.electrophoresis("tapestation/2022-10-23 - 11-08-02-HSRNA.xml")
rna_fresh <- read.electrophoresis("tapestation/2022-08-25 - 10-03-38[14518]-HSRNA.xml")
dna_pool <- read.electrophoresis("tapestation/2023-05-25 - 13-27-33-HSD1000.xml")
rna_OC <- read.electrophoresis("tapestation/2022-08-04 - 11-37-28-HSRNA.xml")

# the actual data is in the data element
rna_time_1$data

# make some plots with the plots built-in to bioanalyzeR
?qplot.electrophoresis
qplot.electrophoresis(rna_time_1)  
qplot.electrophoresis(rna_time_1,  log="x", scales="free_y", y = "fluorescence", x="relative.distance")
qplot.electrophoresis(rna_time_1,  log="x", y = "fluorescence", show.peaks = "markers", normalize = "total")
qplot.electrophoresis(rna_time_1,  log="x", scales="free_y", y = "fluorescence", show.peaks = "markers", normalize = "total")

# make our own plots in ggplot2
df1 <- rna_time_1$data
df2 <- rna_time_2$data
df3 <- rna_fresh$data
df4 <- rna_time_3$data
df5 <- dna_pool$data
df6 <- rna_OC$data
  
# remove the ladder from df1 & df2 and faulty data from df1
df1 <- filter(.data = df1, sample.index != 1 & sample.index != 8 & sample.index != 9)
df2 <- filter(.data = df2, sample.index !=1)
df4 <- filter(.data = df4, sample.index == 3)

# remove ladder & fresh male sample from df3
df3 <- filter(.data = df3, sample.index == 3)

# remove ladder from pool sample from df5
df5 <- filter(.data = df5, sample.index == 2) 

# get a couple of the OC samples from df5
df6 <- filter(.data = df6, sample.index == c(9, 11, 7))
  
# rename sample index from 2 and 3 to 8 and 9 respectively
df2["sample.index"][df2["sample.index"] == 2] <- 8
df2["sample.index"][df2["sample.index"] == 3] <- 9
df4["sample.index"][df4["sample.index"] == 3] <- 11

# rename sample index in df3 from 3 to 1
df3["sample.index"][df3["sample.index"]==3] <- 1

# merge the three data frames to get the complete data set
comb_df <- rbind(df1, df2, df4, df3)

# reorder so that samples are in appropriate order
comb_df <- comb_df %>%  arrange(sample.index)

# 30 seems to be where the lower marker ends
lower_marker_length_max <- 33
                      
over_time <- ggplot(filter(comb_df, length > lower_marker_length_max)) +       
  geom_line(aes(x=length, y=fluorescence, group=sample.index)) +
  # add in coloring under the lines: 
  # length < lower_marker_mlength_max will produce a T/F value, which we can color with scale_fill_manual
  geom_area(aes(x=length, y=fluorescence, fill = length <= lower_marker_length_max)) +
  scale_fill_manual(values=c("grey60", "lightsteelblue")) +
  # get rid of the ugly legend
  theme_classic(base_size=10) +
  # scale_x_log10(lim = c(0,1000)) +
  scale_x_log10() +
  xlab("RNA Length (nt)") +
  ylab("Fluorescence (arbitrary units)") +
  facet_wrap(~sample.index, ncol=1, scales = "free_y") +
  theme(legend.position="none",
        strip.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_blank()) 

over_time
ggsave("plots/RNAVsTime_Bioanalyzer.pdf", height=10, width=7, units="in")
ggsave("plots/RNAVsTime_Bioanalyzer.svg", width=10, height=7, units="in")
ggsave("plots/RNAVsTime_Bioanalyzer.jpg", width=10, height=7, units="in")

# Comparison of OC samples and final pool

# rename the sample indexes
df5["sample.index"][df5["sample.index"] == 2] <- 1
df6["sample.index"][df6["sample.index"] == 7] <- 4
df6["sample.index"][df6["sample.index"] == 11] <- 3
df6["sample.index"][df6["sample.index"] == 9] <- 2

# Combine OC samples with pool sample
comb_df2 <- rbind(df5, df6)

pool_vs_old <- ggplot(filter(comb_df2, length > lower_marker_length_max)) +       
  geom_line(aes(x=length, y=fluorescence, group=sample.index)) +
  # add in coloring under the lines: 
  # length < lower_marker_mlength_max will produce a T/F value, which we can color with scale_fill_manual
  geom_area(aes(x=length, y=fluorescence, fill = length <= lower_marker_length_max)) +
  scale_fill_manual(values=c("grey60", "lightsteelblue")) +
  # get rid of the ugly legend
  theme_classic(base_size=10) +
  # scale_x_log10(lim = c(0,1000)) +
  scale_x_log10(lim = c(32,1400)) +
  xlab("RNA Length (nt)") +
  ylab("Fluorescence (arbitrary units)") +
  facet_wrap(~sample.index, ncol=1, scales = "free_y") +
  theme(legend.position="none",
        strip.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_blank()) 

pool_vs_old
ggsave("plots/PoolVsSample_Bioanalyzer.pdf", height=7, width=10, units="in")
ggsave("plots/PoolVsSample_Bioanalyzer.svg", width=10, height=7, units="in")
ggsave("plots/PoolVsSample_Bioanalyzer.jpg", width=10, height=7, units="in")
