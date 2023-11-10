#!/usr/bin/env Rscript

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
} else {
  # if running via RStudio
  # r_lib_dir = "../lib/R/"
  tsv_input      = "../results/process/collected_alfa_counts.tsv"
  metadata_input = "../refseq/metadata.csv"
}

library(tidyverse)
library(patchwork)

# read in metadata
metadata <- read.delim(metadata_input, sep=",", 
                       header=T, stringsAsFactors = F)
# convert to character type in case an integer-like sample ID (e.g. 1004277)
metadata$sample_id <- as.character(metadata$sample_id)

datasets <- read.delim(file = tsv_input, sep="\t")
colnames(datasets) <- c("sample_id", "type", "counts", "size_in_genome")
# convert to character type in case an integer-like sample ID (e.g. 1004277)
datasets$sample_id <- as.character(datasets$sample_id)

# plot count totals per dataset
datasets_per_sample <- datasets %>% group_by(sample_id) %>% summarize(total_counts = sum(counts))
ggplot(datasets_per_sample) +
  geom_histogram(aes(x=total_counts), bins=60) + 
  scale_x_log10() +
  xlab("# of fly genome mapping reads") +
  ylab("# datasets") +
  theme_bw()

ggsave("fly-mapping_reads_per_dataset.pdf", units="in", height=6.5, width=9)

# filter out datasets with fewer than 1e7 total counts
too_few_counts <- datasets_per_sample %>% filter(total_counts < 1e7) %>% pull(sample_id)

datasets <- datasets %>% filter(!sample_id %in% too_few_counts)

# confirm that metadata exists for all datasets

dataset_names <- datasets %>% 
  group_by(sample_id) %>% summarize() %>% pull(sample_id)
missing_metadata <- !(dataset_names %in% metadata$sample_id)

if (any(missing_metadata)) {
  message("ERROR: missing metadata for the following datasets:")
  cat(dataset_names[missing_metadata])
}

# merge metadata into counts table
datasets <- left_join(datasets, metadata, by="sample_id")

# filter out simulans samples since possible issues with mapping to D. mel genome
# plenty of D. mel samples
datasets <- datasets %>% filter(species != "simulans")

# parse out category,biotype separated by a comma in column #2
type_split <- str_match(datasets$type, "(.*),(.*)")
datasets$category <- type_split[,2]
datasets$biotype  <- type_split[,3]

datasets <- datasets %>% select(-type)

# combine counts for biotypes with >1 category
cats <- datasets %>% group_by(biotype, category) %>% summarize()
summed_datasets <- datasets %>% group_by(sample_id, biotype) %>% summarize(counts = sum(counts))

# use summed values for plots, etc
datasets <- left_join(summed_datasets, metadata, by="sample_id")

# calculate normalized counts 
datasets <- datasets %>% 
  group_by(sample_id) %>% 
  mutate(fractional_counts = counts / sum(counts)) 

# reorder sample factors for display
datasets$sample_id <- fct_reorder(datasets$sample_id, datasets$date_collected, min)

datasets$sample_type <- fct_relevel(datasets$sample_type, "Old_Collection", "Experimental_dried", "Fresh_frozen")


# don't plot negative control dataset (very few reads)
datasets <- datasets %>% filter(sample_type != "negative")

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

datasets$sample_type <- 
  recode(datasets$sample_type, 
         Old_Collection     = "Old\ncollections", 
         Experimental_dried = "Experimental\ndried",
         Fresh_frozen       = "Fresh\nfrozen")

# plot all 
all_p <- ggplot(datasets) +
  geom_jitter (aes(x=sample_type, y=fractional_counts, fill=sample_type), shape=21, color="black", stroke=0.2, size=1.5, height=0, width=0.25, alpha=0.5) +
  geom_boxplot(aes(x=sample_type, y=fractional_counts, color=sample_type), size=0.25, outlier.shape = NA, fill=NA) +
  theme_classic(base_size = 12) + 
  xlab("") +
  ylab("Percent reads mapping to indicated feature type") +
  scale_color_manual(values = fancy_color_scale) +
  scale_fill_manual (values = fancy_color_scale) +
  facet_wrap(~biotype) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.position = "none")

all_p

all_log_p <- all_p + scale_y_log10() 

# all_log_p 

all_both_p <- all_p + all_log_p + plot_layout(ncol = 1)
# all_both_p

ggsave("alfa_all_plot.pdf", width=6.5, height=9, units="in")

dataset_medians <- datasets %>% group_by(sample_type, biotype) %>% summarize(median_fraction = median(fractional_counts),
                                                                             sd_fraction = sd(fractional_counts))


# plot some 

some_sample_types <- c("opposite_strand", "protein_coding", "rRNA", "tRNA")
some_datasets <- datasets %>% filter(biotype %in% some_sample_types)

some_p <- ggplot(some_datasets) +
  geom_jitter (aes(x=sample_type, y=fractional_counts, fill=sample_type), shape=21, color="black", stroke=0.2, size=1.5, height=0, width=0.25, alpha=0.5) +
  geom_boxplot(aes(x=sample_type, y=fractional_counts, color=sample_type), size=0.25, outlier.shape = NA, fill=NA) +
  theme_classic(base_size = 12) + 
  xlab("") +
  ylab("Percent reads mapping to indicated feature type") +
  scale_color_manual(values = fancy_color_scale) +
  scale_fill_manual (values = fancy_color_scale) +
  facet_wrap(~biotype, nrow=1) +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.position = "none")

some_p

ggsave("alfa_some_plot.pdf", width=9, height=4, units="in")
