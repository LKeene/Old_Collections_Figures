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
  R_lib_dir       = args[3]
  R_script_dir    = args[4]
  output_dir      = "./"
} else {
  # if running via RStudio
  # r_lib_dir = "../lib/R/"
  tsv_input      = "../results/process/all_mismatch_counts.txt"
  metadata_input = "../results/collect/collected_metadata.csv"
  R_lib_dir       = "../lib/R/"
  R_script_dir    = "../../scripts/"
  output_dir      = "../results/process/"
}

# these libraries are part of the tidyverse, so will be availabile in the
# tidyverse singularity image we are using (or analogous conda env)
library(tidyverse)

# these libraries are not part of the standard tidyverse, so may have to load it 
# from a specified path
# either from pipeline's R lib dir or from R environment
if (R_lib_dir != "NA") {
  library(rstatix, lib.loc=R_lib_dir)
  library(ggpubr, lib.loc=R_lib_dir)
  library(patchwork, lib.loc=R_lib_dir)
  
} else {
  # in this case assuming these will be installed
  library(rstatix)
  library(ggpubr)
  library(patchwork)
}

# create an output file and text variable to capture text output from this analysis
output_file <-file(paste0(output_dir, "mismatch_text.txt"))
output_text <- ""

# read in common color definitions
source(paste0(R_script_dir, "/plot_colors.R"))
fancy_color_scale <- c(fresh_frozen_color, experimental_dried_color, old_collection_color)


# read in metadata
metadata <- read.delim(metadata_input, sep=",",
                       header=T, stringsAsFactors = F)
# make sure sample_id is character type 
metadata$sample_id <- as.character(metadata$sample_id)

# define subsets of fresh-frozen datasets
metadata <- metadata %>% mutate(ff_type = case_when(
  str_detect(sample_id, "SRX")  ~ "sra_dataset",
  str_detect(sample_id, "\\d+Wk_[FM]") ~ "crisitunity",
  str_detect(sample_id, "Col30") ~ "colony30",
  .default = "pos_ctrl"))


# read in the file with the mismatch data for all datasets
datasets <- read.delim(tsv_input, header=F, sep="\t")
colnames(datasets) <- c("sample_id", "refseq", "position", "ref_base", "read_base", "count", "succeeding_base", "preceeding_base")
# make sure sample_id is character type 
datasets$sample_id <- as.character(datasets$sample_id)

# confirm that metadata exists for all datasets
dataset_names <- datasets %>%
  group_by(sample_id) %>% summarize() %>% pull(sample_id)
missing_metadata <- !(dataset_names %in% metadata$sample_id)

if (any(missing_metadata)) {
  message("WARNING: missing metadata for the following datasets:")
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

# make a composite name for the type of mutation (e.g. "G to A")
datasets <- datasets %>% mutate(mismatch = str_c(ref_base, " to ", read_base))

# replace T with U because RNA
datasets <- datasets %>% mutate(mismatch = str_replace_all(mismatch, "T", "U"))

# totals for each mismatch type
datasets_by_type <- datasets %>% group_by(sample_id, mismatch) %>% summarize(count = sum(count), .groups="drop")

# totals counts for each dataset
datasets_by_type <- datasets_by_type %>% group_by(sample_id) %>% mutate(total_count = sum(count),
                                                                        mismatch_freq = count / total_count)

# remove non-mutations (e.g. A to A)
# datasets_by_type <- datasets_by_type %>% filter(ref_base != read_base)

# merge metadata into counts table
datasets_by_type <- left_join(datasets_by_type, metadata, by="sample_id")

# filter out simulans samples since possible issues with mapping to D. mel genome
# plenty of D. mel samples
datasets_by_type <- datasets_by_type %>% filter(str_detect(species, "melano"))

# reorder sample factors for display
datasets_by_type$sample_id <- fct_reorder(datasets_by_type$sample_id, datasets_by_type$date_collected, min)

# don't plot negative control dataset (very few reads)
datasets_by_type <- datasets_by_type %>% filter(sample_type != "Negative")

# reorder for display
datasets_by_type$sample_type <- fct_relevel(datasets_by_type$sample_type, "Fresh_frozen", "Experimental_dried", "Old_Collection")

# fancier labels for display
datasets_by_type$sample_type <- 
  recode(datasets_by_type$sample_type, 
         Old_Collection     = "Museum\nsamples",
         Experimental_dried = "Experimental\ndried",
         Fresh_frozen       = "Fresh\nfrozen")

# calculate mismatch frequencies relative to median of fresh frozen samples
# first calculate medians
median_mismatch_freqs <- datasets_by_type %>% 
  filter(sample_type == "Fresh\nfrozen") %>% 
  group_by(mismatch) %>% 
  summarize(median_ff_freq = median(mismatch_freq))

# calculate mismatch frequencies relative to median of fresh frozen samples
datasets_by_type <- left_join(datasets_by_type, median_mismatch_freqs)
datasets_by_type <- datasets_by_type %>% mutate(relative_mismatch_freq = mismatch_freq / median_ff_freq)

# --------
# Do stats
# --------

# are data normally distributed?


# this Shapiro test indicates frequencies not normally distributed
df_shapiro     <- datasets_by_type %>% group_by(mismatch, sample_type) %>% shapiro_test(mismatch_freq)
filter(df_shapiro, p<0.05)

# what about log-normally distributed?
datasets_by_type <- datasets_by_type %>% mutate(log_mismatch_freq = log10(mismatch_freq))
df_shapiro_log <- datasets_by_type %>% group_by(mismatch, sample_type) %>% shapiro_test(log_mismatch_freq)
filter(df_shapiro_log, p<0.05)

# what about relative frequencies?
df_shapiro_rel <- datasets_by_type %>% group_by(mismatch, sample_type) %>% shapiro_test(relative_mismatch_freq)
filter(df_shapiro_rel, p<0.05)

# OK - neither normally nor log-normally distributed according to Shapiro test
# will have to use non-parametric tests

# run stats

# Wilcoxon test 
df_wilcox <- datasets_by_type %>% group_by(mismatch) %>% wilcox_test(mismatch_freq ~ sample_type)
df_wilcox_plot <- df_wilcox %>% add_xy_position(x = "sample_type")

# Wilcoxon test on log10 transformed data to match log10 plot
df_wilcox_log <- datasets_by_type %>% group_by(mismatch) %>% wilcox_test(log_mismatch_freq ~ sample_type)
# df_wilcox_log_plot <- df_wilcox_log %>% add_xy_position(x = "sample_type", scales="free_y")
df_wilcox_log_plot <- df_wilcox_log %>% add_xy_position(x = "sample_type")

# Wilcoxon test on relative frequencies 
df_wilcox_rel <- datasets_by_type %>% group_by(mismatch) %>% wilcox_test(relative_mismatch_freq ~ sample_type)
# df_wilcox_log_plot <- df_wilcox_log %>% add_xy_position(x = "sample_type", scales="free_y")
df_wilcox_rel_plot <- df_wilcox_rel %>% add_xy_position(x = "sample_type")


# plot raw mismatch frequencies 
all_p <- ggplot(datasets_by_type) +
  geom_jitter (aes(x=sample_type, y=mismatch_freq, fill=sample_type),  shape=21, color="black", stroke=0.2, size=1.5, height=0, width=0.25, alpha=0.5) +
  geom_boxplot(aes(x=sample_type, y=mismatch_freq, color=sample_type), size=0.25, outlier.shape = NA, fill=NA) +
  theme_this_paper(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5),) +
  scale_y_log10() + 
  scale_color_manual(values = fancy_color_scale) +
  scale_fill_manual (values = fancy_color_scale) +
  xlab("") +
  ylab("Frequency of indicated mutation")  +
  facet_wrap(~mismatch) 

all_p

y_positions <- rep(c(0.15,0.45,0.15), 16)
all_p_stats <- all_p + 
  stat_pvalue_manual(df_wilcox_plot, y.position = y_positions, tip.length = 0.01, bracket.shorten=0.1, size=3, bracket.nudget.y=0.1)
  
all_p_stats

ggsave(paste0(output_dir, "/Fig_SX_raw_mismatch_frequencies.pdf"), height=7.5, width=10, units="in")

# plot relative mismatch frequencies 
all_rel_p <- ggplot(datasets_by_type) +
  geom_jitter (aes(x=sample_type, y=relative_mismatch_freq, fill=sample_type),  shape=21, color="black", stroke=0.2, size=1.5, height=0, width=0.25, alpha=0.5) +
  geom_boxplot(aes(x=sample_type, y=relative_mismatch_freq, color=sample_type), size=0.25, outlier.shape = NA, fill=NA) +
  theme_this_paper(base_size = 12) + 
  scale_y_log10() +
  scale_color_manual(values = fancy_color_scale) +
  scale_fill_manual (values = fancy_color_scale) +
  xlab("") +
  ylab("Frequency of indicated mismatch\nrelative to average frequency in fresh-frozen samples")  +
  facet_wrap(~mismatch) 

all_rel_p

y_positions <- rep(c(2.3,2.5,2.3), 16)
all_rel_p_stats <- all_rel_p + 
  stat_pvalue_manual(df_wilcox_rel_plot, y.position = y_positions, tip.length = 0.01, bracket.shorten=0.1, size=3, bracket.nudget.y=0.1)
  
all_rel_p_stats

ggsave(paste0(output_dir, "/Fig_10_relative_mismatch_frequencies.pdf"), height=7.5, width=10, units="in")


# median frequency of C to U mismatches in old samples vs fresh-frozen
medians_by_mismatch <- datasets_by_type %>%
  group_by(mismatch, sample_type) %>%
  summarize(
    median_freq          = median(mismatch_freq), 
    median_relative_freq = median(relative_mismatch_freq), .groups="drop")

# text for paper
oc_c_to_u_fold_increase <- medians_by_mismatch %>% 
  filter(mismatch == "C to U" & sample_type == "Museum\nsamples") %>% 
  pull(median_relative_freq)

oc_c_to_u_freq <- medians_by_mismatch %>%
  filter(mismatch == "C to U" & sample_type == "Museum\nsamples") %>% 
  pull(median_freq) 

ff_c_to_u_freq <- medians_by_mismatch %>%
  filter(mismatch == "C to U" & sample_type == "Fresh\nfrozen") %>% 
  pull(median_freq)

oc_v_ff_pval <- df_wilcox_rel %>%
  filter(mismatch == "C to U" & group1 == "Museum\nsamples" & group2 == "Fresh\nfrozen")  %>%
  pull(p.adj)

output_text <- 
  paste0(output_text,
         "The median level of C to U mismatches in datasets from museum samples was ",
         sprintf("%0.1f", oc_c_to_u_fold_increase),
         " times higher than in fresh-frozen datasets (",
         sprintf("%0.1e", oc_c_to_u_freq),
         " vs. ",
         sprintf("%0.1e", ff_c_to_u_freq),
         " ; p=",
         sprintf("%0.1e", oc_v_ff_pval),
         ").")


# A-to-G mutations
oc_a_to_g_fold_increase <- medians_by_mismatch %>% 
  filter(mismatch == "A to G" & sample_type == "Museum\nsamples") %>% 
  pull(median_relative_freq)

oc_a_to_g_freq <- medians_by_mismatch %>%
  filter(mismatch == "A to G" & sample_type == "Museum\nsamples") %>% 
  pull(median_freq) 

ff_a_to_g_freq <- medians_by_mismatch %>%
  filter(mismatch == "A to G" & sample_type == "Fresh\nfrozen") %>% 
  pull(median_freq)

oc_v_ff_a_to_g_pval <- df_wilcox_rel %>%
  filter(mismatch == "A to G" & group1 == "Museum\nsamples" & group2 == "Fresh\nfrozen")  %>%
  pull(p.adj)

output_text <- 
  paste0(output_text,
         "The median level of A to G mismatches in museum datasets from museum samples was ",
         sprintf("%0.1f", oc_a_to_g_fold_increase),
         " fold higher than in fresh datasets (",
         sprintf("%0.1e", oc_a_to_g_freq),
         " vs. ",
         sprintf("%0.1e", ff_a_to_g_freq),
         " ; p=",
         sprintf("%0.1e", oc_v_ff_a_to_g_pval),
         ").")
       


# write out text for paper to a file
print(output_text)
writeLines(output_text, output_file)

# close file 
close(output_file)
