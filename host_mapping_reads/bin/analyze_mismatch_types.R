library(tidyverse)

sessionInfo()
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
  tsv_input      = "../results/all_mismatch_counts.txt"
  metadata_input = "../refseq/metadata.csv"
  minimum_total_count  = 10000
}


# read in metadata
metadata <- read.delim(metadata_input, sep=",",
                       header=T, stringsAsFactors = F)


# read in the file with the mismatch data for all datasets
datasets <- read.delim(tsv_input, header=F, sep="\t")
colnames(datasets) <- c("sample_id", "position", "ref_base", "read_base", "count", "succeeding_base", "preceeding_base")

# confirm that metadata exists for all datasets

dataset_names <- datasets %>%
  group_by(sample_id) %>% summarize() %>% pull(sample_id)
missing_metadata <- !(dataset_names %in% metadata$sample_id)

if (any(missing_metadata)) {
  message("ERROR: missing metadata for the following datasets:")
  cat(dataset_names[missing_metadata])
}

datasets <- datasets %>% select(-succeeding_base, preceeding_base)

# make a composite name for the type of mutation (e.g. "G->A")
datasets <- datasets %>% mutate(mismatch = str_c(ref_base, "->", read_base))

# replace T with U because RNA
datasets <- datasets %>% mutate(mismatch = str_replace(mismatch, "T", "U"))

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

# reorder sample factors for display
datasets_by_type$sample_id <- fct_reorder(datasets_by_type$sample_id, datasets_by_type$date_collected, min)


# don't plot negative control dataset (very few reads)
datasets_by_type <- datasets_by_type %>% filter(control_type != "Negative")

# plot raw mismatch frequencies 
ggplot(datasets_by_type) +
# ggplot(filter(datasets_by_type, mismatch == "C->G")) +
  geom_point(aes(x=sample_id, y=mismatch_freq, color=control_type, group=sample_id)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5)) +
  scale_fill_manual(values=c("cornflowerblue", "orange", "grey")) +
  facet_wrap(~mismatch) +
  scale_y_log10() + 
  xlab("") +
  ylab("Frequency of indicated mutation") 

# save a pdf of it
ggsave("mismatch_frequencies.pdf", height=7.5, width=10, units="in")


ggplot(datasets_by_type) +
  geom_boxplot(aes(x=control_type, y=mismatch_freq, group=control_type)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5)) +
  scale_fill_manual(values=c("cornflowerblue", "orange", "grey")) +
  facet_wrap(~mismatch) +
  scale_y_log10() + 
  xlab("") +
  ylab("Frequency of indicated mutation") 
