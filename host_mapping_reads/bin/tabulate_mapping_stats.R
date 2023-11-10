#!/usr/bin/env Rscript

# This code block sets up input arguments to either come from the command line
# (if running from the pipeline, so not interactively) or to use expected default values 
# (if running in interactive mode, for example in RStudio for troubleshooting
# or development).  
#
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  tsv_input        = args[1]
  output_dir       = "./"
} else {
  # if running via RStudio
  args = commandArgs(trailingOnly=TRUE)
  tsv_input        = "../results/all_sam_stats.txt"
  output_dir       = "../results/"
}

library(tidyverse)

# read in metadata
df <- read.delim(file = tsv_input, sep="\t", header=F)
colnames(df) <- c("sample_id", "refseq", "pct_id")

df_stats <- df %>% group_by(sample_id, refseq) %>% summarize(count=n(), mean_pct_id = mean(pct_id, na.rm = T))

write.table(df_stats, file=paste0(output_dir, "/summarized_mapping_stats.txt"), quote=F, sep="\t", row.names=F, col.names=F)

head(df_stats)
