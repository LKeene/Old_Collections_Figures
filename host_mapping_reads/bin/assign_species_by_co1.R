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
} else {
  # if running via RStudio
  # r_lib_dir = "../lib/R/"
  tsv_input      = "../results/process/all_species_tallies.txt"
  metadata_input = "../refseq/metadata.csv"
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
colnames(df) <- c("dataset", "species", "count")

# calculate CO1-mapping count totals for each dataset
# then determine the fractional counts for CO1-mapping reads
# then only keep any counts that represent >1% of the datasets
df <- df %>% 
  group_by(dataset) %>% 
  mutate(total_count = sum(count), 
         fractional_count = count / total_count) %>% 
  filter(fractional_count > 0.1) %>%
  mutate(accession_species = species,  species = str_match(species, "Drosophila_.*")[,1])


# top hits for each datset
top_hits <- 
  df %>% 
  group_by(dataset) %>% 
  arrange(-fractional_count) %>% 
  filter(row_number()==1) 

ggplot(top_hits) +
  geom_jitter(aes(x=species, y=fractional_count), height=0, width=0.25) +
  theme_bw()

# what species were observed?
observed_species <- top_hits %>% 
  group_by(species) %>% 
  summarize(n=n())

observed_species %>% arrange(-n)

# what about Lexi's flies from FoCo20/21
Lexi_FoCos <- df %>% filter(str_detect(dataset, "^FoCo2"))

df %>% filter(str_match(dataskos))
  
# what species were observed with high confidence?
observed_species <- top_hits %>% 
  group_by(species) %>% 
  summarize(n=n())
  
# actually make assignments if possible
df_assigned <- df %>% 
  group_by(dataset) %>% 
  arrange(-fractional_count) %>% 
  summarize(fraction_co1 = max(fractional_count), 
            count=count[1], 
            assigned_species = if_else(fraction_co1 > 0.85, species[1], "undetermined"))

# write assignments to an excel file
# openxlsx::write.xlsx(df_assigned, file = "co1_based_species_assignments.xlsx")

# write output to tsv file
write.table(df_assigned, file = "co1_based_species_assignments.txt", quote=F, sep="\t", row.names=F)
