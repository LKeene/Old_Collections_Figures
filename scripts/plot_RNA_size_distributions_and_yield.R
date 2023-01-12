library(bioanalyzeR)
library(tidyverse)
library(readxl)
library(patchwork)

rna_time_1 <- read.electrophoresis("2022-08-05 - 11-17-40-HSRNA.xml")
rna_time_2 <- read.electrophoresis("2022-08-05 - 11-49-24-HSRNA.xml")
rna_fresh <- read.electrophoresis("2022-08-25 - 10-03-38[14518]-HSRNA.xml")
oc_samples <- read.electrophoresis("2022-08-04 - 11-37-28-HSRNA.xml")

# read in sample metadata for each tape
metadata <- read_excel("tape_metadata_tidy.xlsx")

# pull out the actual data
df1 <- rna_time_1$data
df2 <- rna_time_2$data
df3 <- rna_fresh$data
df_oc <- oc_samples$data

# these values match tape_id column in metadtata excel
df1$tape_id <- "1"
df2$tape_id <- "2"
df3$tape_id <- "fresh"
df_oc$tape_id <- "oc"

# merge the data from different tapes
df <- rbind(df1, df2, df3, df_oc)

# make a joint ID that combines tape_id and sample.index 
# this will uniquely identify each lane and allow us to link
# tape data to metadata
df <- df %>% mutate(id = str_c(tape_id, sample.index))
metadata <- metadata %>% mutate(id = str_c(tape_id, sample.index))

# get rid of lanes for which no metadata provided
df <- filter(df, id %in% metadata$id)

# check we have all the tapes/samples represented as expected
df %>% group_by(tape_id,sample.index) %>% summarize()

# determine cutoff for where the lower marker ends - then check visually
lower_marker_length_max <- 33
                      
ggplot(filter(df, tape_id == "1")) +       
  geom_line(aes(x=length, y=fluorescence, group=sample.index)) +
  # add in coloring under the lines: 
  # length < lower_marker_mlength_max will produce a T/F value, which we can color with scale_fill_manual
  geom_area(aes(x=length, y=fluorescence, fill = length <= lower_marker_length_max)) +
  scale_fill_manual(values=c("grey50", "lightsteelblue")) +
  # get rid of the ugly legend
  theme_classic(base_size=10) +
  scale_x_log10(limits=c(NA,1000)) +
  xlab("RNA Length (nt)") +
  ylab("Fluorescence (arbitrary units)") +
  facet_wrap(~sample.index, ncol=1) +
  theme(legend.position="none",
        strip.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_blank()) 


# calculate mean from what is essentially histogram data.
# follow the strategy described here:
# https://www.statology.org/histogram-mean-median/

# get rid of rows with no assigned length or below lower marker cutoff
df <- df %>% filter(!is.na(length) & length > lower_marker_length_max)

# only keep fluorescence values if above an arbitrary threshold 
# this is because legit fluoreseence signals are much larger
# and we want to avoid fluorescence signal noise from impacting mean, which it was doing
min_fluorescence_cutoff <- 10
df <- df %>% filter(fluorescence > min_fluorescence_cutoff)

# calculate fluorescence * length for each bin
df <- df %>% mutate(mn = fluorescence * length)

# calculate mean length from each TS trace data.
# and put in a new dataframe with just one value per sample
df_mean <- df %>% group_by(id) %>% summarize(mean_length = sum(mn) / sum(fluorescence))

# merge in metadata
df_mean <- left_join(df_mean, metadata, by="id")

# investigate filtering impact
ggplot(filter(df, tape_id=="oc" & sample.index < 8)) +
  facet_wrap(~id, scales="free_y", ncol=1) +
  scale_x_log10() + 
  geom_point(aes(x=length, y=fluorescence)) 


# plot oc samples: we'll show these separately from new/experimental samples
p_oc <- ggplot(filter(df_mean, new_or_old == "old")) +
  # geom_smooth(aes(x=year, y=mean_length), alpha=0.25, color="slateblue", size=0.5) +
  geom_point(aes(x=year, y=mean_length), shape=21, size=4, fill="slateblue", color="black", stroke=0.25) +
  scale_x_reverse(breaks = seq(from=1880, to=2020, by=20))+
  ylim(c(0,1000)) +
  xlab("Year of sample collection") +
  # we don't need a y-axis label because it's repeated in p_new below
  # ylab("Mean length of RNA on Tapestation (nt)") +
  ylab("") +
  theme_minimal(base_size = 16) 

p_oc
  
# plot new experimental samples
p_new <- ggplot(filter(df_mean, new_or_old=="new")) +
  # geom_smooth(aes(x=as.numeric(weeks), y=mean_length), alpha=0.25, color="slateblue", size=0.5) +
  geom_point(aes(x=as.numeric(weeks), y=mean_length), shape=21, size=4, fill="slateblue", color="black", stroke=0.25) +
  ylim(c(0,1000)) +
  xlab("Weeks since sample pinning") +
  ylab("Mean length of RNA on Tapestation (nt)") +
  theme_minimal(base_size = 16) 

p_new

# plot the two plots side by side
# lay out combined plot using patchwork library
p_new + p_oc +  plot_layout(widths = c(1, 2))

ggsave("RNA_length_vs_time.pdf", width=10, height=7, units="in")


# -----------------------------
# plot some of the sample data
# -----------------------------
# things like RNA recovery vs. time or extraction method...

# read in csv
sample_data <- read.delim("LocationData.csv", sep=",", header=T)

# yield vs.  year
# a color-blind palette
# from: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/ 
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# compare yield of extraction samples: destructive vs. not 
dried_sample_data <- filter(sample_data, Storage.Type =="Dried")

p_method <- ggplot(dried_sample_data) +
  geom_jitter(aes(x=Extraction.Method, y=Concentration, fill=Extraction.Method), 
              width=0.1, height=0, shape=21, size=2.5, color="black", stroke=0.25) +
  geom_boxplot(aes(x=Extraction.Method, y=Concentration), 
               color="grey40", outlier.shape = NA, fill=NA, width=0.25, size=0.25) +
  scale_fill_manual(values=c(cbPalette[6], cbPalette[7])) +
  ylab("Concentration of RNA extract (ng/µL)") +
  xlab("Extraction method") + 
  theme_minimal(base_size = 14)  +
  theme(legend.position = "none")

p_method
ggsave("RNA_yield_vs_extraction_method.pdf", width=10, height=7, units="in")

# compare yield of extraction samples: dried vs. ethanol
p_storage <- ggplot(sample_data) +
  geom_jitter(aes(x=Storage.Type, y=Concentration, fill=Storage.Type), 
              width=0.1, height=0, shape=21, size=2.5, color="black", stroke=0.25) +
  geom_boxplot(aes(x=Storage.Type, y=Concentration), 
               color="grey40", outlier.shape = NA, fill=NA, width=0.25, size=0.25) +
  scale_fill_manual(values=c(cbPalette[2], cbPalette[3])) +
  ylab("Concentration of RNA extract (ng/µL)") +
  xlab("") +
  theme_minimal(base_size = 14)  +
  theme(legend.position = "none")

p_storage
ggsave("RNA_yield_vs_storage_method.pdf", width=10, height=7, units="in")

# plot RNA Yield as a function of sample age
p_year <- ggplot(sample_data) +
  geom_point(aes(x=Date.Collected, y=Concentration), 
              shape=21, size=3, fill="slateblue", color="black", stroke=0.25) +
  # scale_fill_manual(values=c(cbPalette[6], cbPalette[7])) +
  scale_x_reverse() +
  ylab("Concentration of RNA extract (ng/µL)") +
  xlab("Year of sample collection") +
  theme_minimal(base_size = 14) 

p_year
ggsave("RNA_yield_vs_sample_age.pdf", width=10, height=7, units="in")



