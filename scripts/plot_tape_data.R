library(bioanalyzeR)
library(tidyverse)

rna_time_1 <- read.electrophoresis("2022-08-05 - 11-17-40-HSRNA.xml")
rna_time_2 <- read.electrophoresis("2022-08-05 - 11-49-24-HSRNA.xml")
rna_fresh <- read.electrophoresis("2022-08-25 - 10-03-38[14518]-HSRNA.xml")

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

# remove the ladder from df1 & df2 and faulty data from df1
df1 <- filter(.data = df1, sample.index != 1 & sample.index != 8 & sample.index != 9)
df2 <- filter(.data = df2, sample.index !=1)

# remove ladder & fresh male sample from df3
df3 <- filter(.data = df3, sample.index == 3)

# rename sample index from 2 and 3 to 8 and 9 respectively
df2["sample.index"][df2["sample.index"] == 2] <- 8
df2["sample.index"][df2["sample.index"] == 3] <- 9

# rename sample index in df3 from 3 to 1
df3["sample.index"][df3["sample.index"]==3] <- 1

# merge the three data frames to get the complete data set
comb_df <- rbind(df1, df2, df3)

# reorder so that samples are in appropriate order
comb_df <- comb_df %>%  arrange(sample.index)

# 30 seems to be where the lower marker ends
lower_marker_length_max <- 33
                      
ggplot(filter(comb_df, length > lower_marker_length_max)) +       
  geom_line(aes(x=length, y=fluorescence, group=sample.index)) +
  # add in coloring under the lines: 
  # length < lower_marker_mlength_max will produce a T/F value, which we can color with scale_fill_manual
  geom_area(aes(x=length, y=fluorescence, fill = length <= lower_marker_length_max)) +
  scale_fill_manual(values=c("grey50", "lightsteelblue")) +
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

ggsave("RNAVsTime_Bioanalyzer.pdf", height=10, width=7, units="in")

# do the same for the Old Collection Samples
oc_samples <- read.electrophoresis("2022-08-04 - 11-37-28-HSRNA.xml")

qplot.electrophoresis(oc_samples)  
qplot.electrophoresis(oc_samples,  log="x", scales="free_y", y = "fluorescence", x="relative.distance")
qplot.electrophoresis(oc_samples,  log="x", y = "fluorescence", show.peaks = "markers", normalize = "total")
qplot.electrophoresis(oc_samples,  log="x", scales="free_y", y = "fluorescence", show.peaks = "markers", normalize = "total")

df3 <- oc_samples$data

df3 <- filter(.data = df3, sample.index != 1)

ggplot(df3) +       
  geom_line(aes(x=length, y=fluorescence, group=sample.index)) +
  # add in coloring under the lines: 
  # length < lower_marker_mlength_max will produce a T/F value, which we can color with scale_fill_manual
  geom_area(aes(x=length, y=fluorescence, fill = length <= lower_marker_length_max)) +
  scale_fill_manual(values=c("grey50", "lightsteelblue")) +
  # get rid of the ugly legend
  theme_classic(base_size=10) +
  # scale_x_log10(lim = c(0,1000)) +
  scale_x_log10(limits=c(NA,1000)) +
  xlab("RNA Length (nt)") +
  ylab("Fluorescence (arbitrary units)") +
  facet_wrap(~sample.index, ncol=1) +
  theme(legend.position="none",
        strip.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_blank()) 

ggsave("OC_Bioanalyzer.pdf", height=10, width=7, units="in")
