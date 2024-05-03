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
  output_dir      = "./"
} else {
  # if running via RStudio
  # r_lib_dir = "../lib/R/"
  tsv_input      = "../results/process/transcriptome_tallies.txt"
  metadata_input = "../results/collected_metadata.csv"
  output_dir     = "../results/"
}

#
#
# MDS 5/6/2020
#

# read in the data
df <- read.delim(tsv_input, sep="\t", header=F)


