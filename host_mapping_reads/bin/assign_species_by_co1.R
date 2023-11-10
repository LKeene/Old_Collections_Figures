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
  tsv_input      = "../results_species_id/process/summarized_mapping_stats.txt"
  metadata_input = "../refseq/metadata.csv"
  output_dir     = "../results_species_id/"
}

#
# assign Drosophila species based on counts of CO1-mapping reads
# from a competitive mapping analysis using all available
# Drosophila CO1 (cytochrome oxidase 1) sequencs
#
# MDS 5/6/2020
#

# read in the data
df <- read.delim(tsv_input, sep="\t", header=F)

# name columns
colnames(df) <- c("dataset", "species", "count", "pct_id")

# calculate CO1-mapping count totals for each dataset
# then determine the fractional counts for CO1-mapping reads
# then only keep any counts that represent >1% of the datasets
df <- df %>% 
  group_by(dataset) %>% 
  mutate(total_count = sum(count), 
         fractional_count = count / total_count) %>% 
  filter(fractional_count > 0.1) %>%
  mutate(accession_species = species,  species = str_match(species, "_([A-Z].*)")[,2])


# top hits for each datset
top_hits <- 
  df %>% 
  group_by(dataset) %>% 
  arrange(-fractional_count) %>% 
  filter(row_number()==1) 

# what species were observed?
observed_species <- top_hits %>% 
  group_by(species) %>% 
  summarize(n=n()) %>%
  arrange(-n)

# output observed spp
write.table(observed_species, 
            file=paste0(output_dir, "co1_observed_species.txt"), 
            quote=F, sep="\t", row.names=F)

# what about Lexi's flies from FoCo20/21
Lexi_FoCos <- df %>% filter(str_detect(dataset, "^FoCo2"))

# actually make assignments based on highest fractional count
df_assigned <- df %>% 
  group_by(dataset) %>% 
  arrange(-fractional_count) %>% 
  filter(row_number()==1) %>%
  mutate(assigned_species = if_else(fractional_count > 0.50, species, "undetermined")) %>%
  select(dataset, assigned_species, count, fractional_count, pct_id) %>%
  arrange(dataset)

  # summarize(fraction_co1 = max(fractional_count), 
            # count=count[1], 
            # pct_id = pct_id[1],
            # assigned_species = if_else(fraction_co1 > 0.50, species[1], "undetermined"))

# write output to tsv file
write.table(df_assigned, 
            file = paste0(output_dir, "co1_based_species_assignments.txt"), 
            quote=F, sep="\t", row.names=F)
