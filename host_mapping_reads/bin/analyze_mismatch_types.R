#!/usr/bin/env Rscript 

library(tidyverse)

# This code block sets up input arguments to either come from the command line
# (if running from the pipeline, so not interactively) or to use expected default values
# (if running in interactive mode, for example in RStudio for troubleshooting
# or development).
#
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  # lib_dir=args[1]
  tsv_input       = args[1]
  metadata_input  = args[2]
  minimum_total_count  = args[3]
} else {
  # if running via RStudio
  # r_lib_dir = "../lib/R/"
  tsv_input      = "../results/process/all_mismatch_counts.txt"
  metadata_input = "../refseq/metadata.csv"
  minimum_total_count  = 10000
}


# read in metadata
metadata <- read.delim(metadata_input, sep=",",
                       header=T, stringsAsFactors = F)


# read in the file with the mismatch data for all datasets
datasets <- read.delim(tsv_input, header=F, sep="\t")
colnames(datasets) <- c("sample_id", "refseq", "position", "ref_base", "read_base", "count", "succeeding_base", "preceeding_base")

# confirm that metadata exists for all datasets

dataset_names <- datasets %>%
  group_by(sample_id) %>% summarize() %>% pull(sample_id)
missing_metadata <- !(dataset_names %in% metadata$sample_id)

if (any(missing_metadata)) {
  message("ERROR: missing metadata for the following datasets:")
  cat(dataset_names[missing_metadata])
}

datasets <- datasets %>% select(-succeeding_base, preceeding_base)

# filter out datasets with not enough total counts
dataset_counts <- datasets %>% group_by(sample_id) %>% summarize(total_counts = sum(count))

# plot # of counts per dataset
ggplot(dataset_counts) + geom_histogram(aes(x=total_counts), bins=60) + scale_x_log10() 

# filter out datasets with fewer than 1e7 total counts
too_few_counts <- dataset_counts %>% filter(total_counts < 1e7) %>% pull(sample_id)
datasets <- datasets %>% filter(!sample_id %in% too_few_counts)

# make a composite name for the type of mutation (e.g. "G->A")
datasets <- datasets %>% mutate(mismatch = str_c(ref_base, "->", read_base))

# replace T with U because RNA
datasets <- datasets %>% mutate(mismatch = str_replace_all(mismatch, "T", "U"))

# totals for each mismatch type
datasets_by_type <- datasets %>% group_by(sample_id, mismatch) %>% summarize(count = sum(count), .groups="drop")

# totals counts for each dataset
datasets_by_type <- datasets_by_type %>% group_by(sample_id) %>% mutate(total_count = sum(count),
                                                                        mismatch_freq = count / total_count)

# remove non-mutations (e.g. A->A)
# datasets_by_type <- datasets_by_type %>% filter(ref_base != read_base)

# only keep positions with sufficient total coverage
# datasets <- datasets %>% filter(total_count > minimum_total_count)

# what is the range?
# min_freq <- min(df$mismatch_freq)
# max_freq <- max(df$mismatch_freq)
# mean_freq <- mean(df$mismatch_freq)
# median_freq <- median(df$mismatch_freq)

# merge metadata into counts table
datasets_by_type <- left_join(datasets_by_type, metadata, by="sample_id")

# filter out simulans samples since possible issues with mapping to D. mel genome
# plenty of D. mel samples
datasets_by_type <- datasets_by_type %>% filter(species != "simulans")

# reorder sample factors for display
datasets_by_type$sample_id <- fct_reorder(datasets_by_type$sample_id, datasets_by_type$date_collected, min)

# don't plot negative control dataset (very few reads)
datasets_by_type <- datasets_by_type %>% filter(sample_type != "Negative")

# reorder for display
datasets_by_type$sample_type <- fct_relevel(datasets_by_type$sample_type, "Old_Collection", "Experimental_dried", "Fresh_frozen")

# fancier text
datasets_by_type$sample_type <- 
  recode(datasets_by_type$sample_type, 
         Old_Collection     = "Old\ncollections", 
         Experimental_dried = "Experimental\ndried",
         Fresh_frozen       = "Fresh\nfrozen")

#RColorBrewer Dark2 colors
dark2_lightblue   <- rgb(166/255,206/255,227/255)
dark2_blue        <- rgb(31/255,120/255,180/255)
dark2_lightgreen  <- rgb(178/255,223/255,138/255)
dark2_green       <- rgb(51/255,160/255,44/255)
dark2_pink        <- rgb(251/255,154/255,153/255)
dark2_red         <- rgb(227/255,26/255,28/255)
dark2_lightorange <- rgb(253/255,191/255,111/255)
dark2_orange      <- rgb(255/255,127/255,0/255)
dark2_lightpurple <- rgb(202/255,178/255,214/255)
dark2_purple      <- rgb(106/255,61/255,154/255)
dark2_yellow      <- rgb(255/255,255/255,153/255)
dark2_brown       <- rgb(177/255,89/255,40/255)

fancy_color_scale <- c(dark2_orange, dark2_blue, dark2_green)

# plot raw mismatch frequencies 
ggplot(datasets_by_type) +
  geom_jitter (aes(x=sample_type, y=mismatch_freq, fill=sample_type),  shape=21, color="black", stroke=0.2, size=1.5, height=0, width=0.25, alpha=0.5) +
  geom_boxplot(aes(x=sample_type, y=mismatch_freq, color=sample_type), size=0.25, outlier.shape = NA, fill=NA) +
  theme_classic(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5),) +
  theme(legend.position = "none") +
  scale_y_log10() + 
  scale_color_manual(values = fancy_color_scale) +
  scale_fill_manual (values = fancy_color_scale) +
  xlab("") +
  ylab("Frequency of indicated mutation")  +
  facet_wrap(~mismatch) 

# save a pdf of plot
ggsave("mismatch_frequencies.pdf", height=7.5, width=10, units="in")

# just mismatches involving Cs
c_mismatches <- datasets_by_type %>% filter(str_detect(mismatch, "C->"))

# plot raw mismatch frequencies 
ggplot(c_mismatches) +
  geom_jitter (aes(x=sample_type, y=mismatch_freq, fill=sample_type),  shape=21, color="black", stroke=0.2, size=1.5, height=0, width=0.25, alpha=0.5) +
  geom_boxplot(aes(x=sample_type, y=mismatch_freq, color=sample_type), size=0.25, outlier.shape = NA, fill=NA) +
  theme_classic(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5),) +
  # theme(axis.text.x = element_text(size = 10)) +
  theme(legend.position = "none") +
  scale_y_log10() + 
  scale_color_manual(values = fancy_color_scale) +
  scale_fill_manual (values = fancy_color_scale) +
  xlab("") +
  ylab("Mismatch frequency in fly-mapped reads")  +
  facet_wrap(~mismatch, nrow=1) 

ggsave("c_mismatches.pdf", units="in", height=4, width=7)

# medians of groups
dataset_medians <- datasets_by_type %>% 
  group_by(sample_type, mismatch) %>% 
  summarize(median_freq = median(mismatch_freq),
                                                sd_freq     = sd(mismatch_freq))


