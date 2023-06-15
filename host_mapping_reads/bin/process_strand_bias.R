#!/usr/bin/env Rscript 


library(tidyverse)

# sessionInfo()

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
} else {
  # if running via RStudio
  # r_lib_dir = "../lib/R/"
  tsv_input      = "../results/collected_strand_bias.tsv"
}

# read in the file with the mismatch data for all datasets
datasets <- read.delim(tsv_input, header=F, sep="\t")
colnames(datasets) <- c("sample_id", "sam_file", "refseq", "fraction_negative_strand")

datasets <- datasets %>% mutate(fraction_positive_strand = 1-fraction_negative_strand)

# extract metadata out of sample_id values
metadata <- str_match(datasets$sample_id, "(\\d+)WkDry_([MF])_([\\d+]_)")

datasets$week <- as.numeric(metadata[,2])
datasets$replicate <- metadata[,4]

# average replicates
datasets_avg <- datasets %>% 
  group_by(week, refseq) %>% 
  summarize(mean_frac_pos = mean(fraction_positive_strand),
            sd_frac_pos = sd(fraction_positive_strand),
            .groups="drop")

# plot raw mismatch frequencies 
ggplot(datasets_avg) +
  geom_point(aes(x=week, y=mean_frac_pos, fill=refseq), shape=21) + 
  geom_errorbar(aes(x=week, ymin=mean_frac_pos-sd_frac_pos, ymax=mean_frac_pos+sd_frac_pos, color=refseq), width=0.5, size=0.5) + 
  theme_classic() +
  # theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5)) +
  facet_wrap(~refseq) +
  # scale_y_log10() + 
  # ylim(c(0,1)) +
  xlab("Weeks samples dried") +
  ylab("Fraction of reads +strand mapping")

ggsave("strand_bias_in_experimental_samples.pdf", width=10, height=7.5, units="in")
