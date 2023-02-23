


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
  tsv_input      = "../results/collected_rRNA_locus_coverage_individual_positions.tsv"
  metadata_input = "../refseq/metadata.csv"
}

library(tidyverse)
library(patchwork)

# read in metadata
metadata <- read.delim(metadata_input, sep=",", 
                       header=T, stringsAsFactors = F)

coverage <- read.delim(file = tsv_input, sep="\t", header=F)
coverage_columns <- c("sample_id", "refseq", 
                      "position", 
                      "fwd_coverage", "rev_coverage")
colnames(coverage) <- coverage_columns

# flesh out coverage df, replacing intervals (with same coverage), with individual values
min_x <- min(coverage %>% group_by(position) %>% pull(position))
max_x <- max(coverage %>% group_by(position) %>% pull(position))

# confirm that metadata exists for all datasets

dataset_names <- coverage %>% 
  group_by(sample_id) %>% summarize() %>% pull(sample_id)
missing_metadata <- !(dataset_names %in% metadata$sample_id)

if (any(missing_metadata)) {
  message("ERROR: missing metadata for the following datasets:")
  cat(dataset_names[missing_metadata])
}

# merge metadata into counts table
coverage <- left_join(coverage, metadata, by="sample_id")


# calculate ratio of Fwd/Rev coverage
coverage <- coverage %>% mutate(cov_ratio = fwd_coverage/rev_coverage)


# reorder sample factors for display
coverage$sample_id <- fct_reorder(coverage$sample_id, coverage$date_collected, min)

# don't plot negative control dataset (very few reads)
coverage <- coverage %>% filter(control_type != "Negative")

# create coverage windows 
window_size <- 20
coverage <- coverage %>% mutate(window_start = window_size * (position %/% window_size),
                         window_mid   = window_start + (window_size / 2))

# calculate % base paired in each window
coverage_windows <- coverage %>% group_by(sample_id, window_mid) %>% 
  summarize(mean_fwd_cov = mean(fwd_coverage), .groups="drop")

# merge metadata into counts table
coverage_windows <- left_join(coverage_windows, metadata, by="sample_id")

# read in base-pairing data: a file that has whether individual bases
# in the rRNA reference sequence are predicted to be involved in base pairing
bp <- read.delim("../refseq/rdna_bp.txt", sep="\t")
colnames(bp) <- c("position", "in_bp")

bp <- bp %>% mutate(window_start = window_size * (position %/% window_size),
                    window_mid   = window_start + (window_size / 2))

# calculate % base paired in each window
bp_window_values <- bp %>% group_by(window_mid) %>% summarize(pct_bp = sum(in_bp)/window_size)

# flesh out values in bp_windows: need NA values for windows that have no data
possible_windows <- seq(from = 0 + (window_size / 2), to = max_x, by=window_size)

# create a blank DF with all possible window midpoints
bp_windows <- data.frame(matrix(ncol = 2, nrow = length(possible_windows)))
colnames(bp_windows) <- c("window_mid", "percent_base_paired")
bp_windows$window_mid <- possible_windows
bp_windows$percent_base_paired <- NA_real_

bp_windows <- left_join(bp_windows, bp_window_values, by="window_mid")

# use NA or the real value if available
bp_windows <- bp_windows %>% mutate(percent_base_paired = 
                        if_else(is.na(pct_bp), 
                                percent_base_paired, 
                                pct_bp)) %>% select(-pct_bp)


# merge bp & cov
cov_bp_windows <- left_join(coverage_windows, bp_windows, by="window_mid")


oldest <- filter(coverage, sample_id =="1004284")

ggplot(coverage) +
  geom_line(aes(x=interval_start, y=fwd_coverage, color=sample_id)) +
  scale_y_log10() +
  facet_wrap(~sample_id) +
  theme_bw() 

# some_dataset_ids <- c("Albany1902_1", "1004279", "1004283", "Davidson2006_1", "FoCo17_Pos", "PosCtrl_Pool1")
some_dataset_ids <- c("Albany1902_1", "1004283","PosCtrl_Pool1")
some_datasets <- filter(coverage, sample_id %in% some_dataset_ids)
cov_plot <- ggplot(some_datasets) +
  geom_line(aes(x=position, y=fwd_coverage, color=sample_id)) +
  scale_y_log10() +
  facet_wrap(~sample_id, ncol=1) +
  # facet_wrap(~sample_id) +
  xlab("") + 
  ylab ("Coverage of sense mapping reads") +
  theme_bw(base_size=14) +
  theme(legend.position = "none")
  
cov_plot
  
bp_plot <- ggplot (bp_windows) +
  geom_line(aes(x=window_mid, y=percent_base_paired)) +
  theme_bw(base_size = 14) +
  ylab(paste0("Fraction of ", window_size, " nt windows\nin structure\nor base-paired")) +
  xlab("Position in pre-rRNA (nt)") 

bp_plot

cov_plot + bp_plot  + plot_layout(ncol=1, heights=c(8,2))

ggsave("coverage_bp_plot.pdf", units="in", height=7.5, width=10)

# scatter plot of cov vs. %basepaired
some_datasets <- filter(cov_bp_windows, sample_id %in% some_dataset_ids)

ggplot(some_datasets) +
  geom_point(aes(x=percent_base_paired, y=mean_fwd_cov, fill=sample_id), 
             shape=21, size=1, color="black", stroke=0.25) +
  geom_smooth(aes(x=percent_base_paired, y=mean_fwd_cov)) +
  theme_bw() +
  # scale_x_log10() +
  scale_y_log10() +
  xlab("Percent of window base paired or in structure") +
  ylab("Mean sense coverage(x)") +
  facet_wrap(~sample_id)

ggplot(coverage) +
  geom_line(aes(x=interval_start, y=rev_coverage, color=sample_id)) +
  scale_y_log10() +
  facet_wrap(~sample_id)
  theme_bw() 

ggplot(coverage) +
  geom_line(aes(x=interval_start, y=cov_ratio, color=sample_id)) +
  scale_y_log10() +
  facet_wrap(~sample_id)
  theme_bw() 

ggplot(filter(coverage, is.finite(cov_ratio))) +
  # geom_jitter(aes(x=sample_id, y=cov_ratio, fill=control_type), shape=21, color="black", stroke=0.25) +
  geom_violin(aes(x=sample_id, y=cov_ratio, fill=control_type), color="black", stroke=0.25) +
  scale_y_log10() +
  ylab("Ratio of +strand mapping to -strand mapping rRNA reads")
  theme_bw() 


# individual base data
individual_positions <- left_join(coverage, bp)
individual_positions <- individual_positions %>% filter(!is.na(in_bp))

some_individual_positions <- filter(individual_positions, sample_id %in% some_dataset_ids)

ggplot(some_individual_positions) +
  geom_point(aes(x=in_bp, y=fwd_coverage)) +
  geom_boxplot(aes(x=in_bp, y=fwd_coverage, group=in_bp)) +
  theme_bw() + 
  scale_y_log10() +
  facet_wrap(~sample_id)
