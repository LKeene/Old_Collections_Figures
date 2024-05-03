#!/usr/bin/env Rscript

# This script reads in sample metadata from several sources
# and writes out a collected unified metadata file
#
# Some cleanup of metadata is required because of non-tidy column names in 
# LocationData.csv (which also has the OC sample metadata)
#
# and because of columns missing in fresh frozen and experimental dried 
# metadata files (because some columns are only relevant to museum samples,
# and so aren't included in those metadata files)
# 
# MDS 3/1/2024

# if running from command line
if (!interactive()) {
  # access metadata directory via commandArgs
  args = commandArgs(trailingOnly=TRUE)
  metadata_dir <- args[1]
} else {
  # if running via RStudio
  metadata_dir   = "../metadata/"
}

library(tidyverse)

# read in old collection samples (museum collections) metadata
oc_md <- read.delim(paste0(metadata_dir, "/LocationData.csv"), sep=",", header=T)

# rename non-tidy column names
oc_md <- oc_md %>% rename(species = Species,
                          date_collected = Date.Collected,
                          location = Location,
                          donating_institution = Donating.Institution,
                          biosample = BioSample,
                          sample_id_in_paper = Sample.ID,
                          museum_id = Museum.ID,
                          sample_id = Fastq.ID, 
                          storage_type = Storage.Type,
                          extraction_method = Extraction.Method,
                          RNA_concentration = Concentration..ng.μl.,
                          X260_280 = X260.280,
                          X260_230 = X260.230,
                          galbut_virus_positive_qPCR = Galbut)


# R1_strand refers to the orientation of R1 reads following strand-specific library prep 
# different adapters result in different R1 orientations
# see: https://www.idtdna.com/pages/support/faqs/can-the-xgen-unique-dual-index-umi-adapters-be-used-for-rna-seq
# all the OC sequences were prepped with IDT UDI/UMI adapters
oc_md$R1_strand = "reverse"
oc_md$weeks_dried = NA_integer_
oc_md$sample_type = "Old_Collection"

# read in metadata for experimental dried samples 
ed_md <- read.delim(paste0(metadata_dir, "/experimental_dried_metadata.csv"), sep=",", header=T)

# read in metadata for fresh frozen samples 
ff_md <- read.delim(paste0(metadata_dir, "/fresh_frozen_metadata.csv"), sep=",", header=T)

ed_md

# rbind dried & fresh frozen metadata
fe_md <- rbind(ed_md, ff_md)

# add in columns missing from dried/frozen metadata
fe_md$biosample = NA_character_
fe_md$sample_id_in_paper = NA_character_
fe_md$museum_id = NA_character_
fe_md$lat  = NA_integer_
fe_md$long = NA_integer_

# reorder columns in alphabetical order so we can rbind 
oc_md <- oc_md %>% select(order(colnames(oc_md)))
fe_md <- fe_md %>% select(order(colnames(fe_md)))

# double check our column names are OK before binding
stopifnot(all(colnames(oc_md) == colnames(fe_md)))

# this is the final combined metadata
metadata <- rbind(oc_md, fe_md)

# sample type should be a factor
metadata$sample_type <- as.factor(metadata$sample_type)

# write out the collected metadata file
write.table(metadata, "collected_metadata.csv", sep=",", row.names=F, col.names=T)

