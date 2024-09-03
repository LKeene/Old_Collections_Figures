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
  output_dir   <- "."
} else {
  # if running via RStudio
  metadata_dir   <- "../metadata/"
  output_dir     <- metadata_dir
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
                          museum_id = Museum.Accession,
                          sample_id = Fastq.ID, 
                          storage_type = Storage.Type,
                          extraction_method = Extraction.Method,
                          RNA_concentration = Concentration..ng.μl.,
                          X260_280 = X260.280,
                          X260_230 = X260.230,
                          galbut_virus_positive_qPCR = Galbut.virus.RT.qPCR.Result)


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

# read in metadata for positive and negative control samples 
control_md  <- read.delim(paste0(metadata_dir, "/control_metadata.csv"), sep=",", header=T)

# rbind dried & fresh frozen & control metadata
fe_md <- rbind(ed_md, ff_md, control_md)

# add in columns missing from dried/frozen metadata
fe_md$biosample = NA_character_
fe_md$museum_id = NA_character_
fe_md$lat  = NA_integer_
fe_md$long = NA_integer_

# handle sample IDs for paper: if defined, used that
fe_md <- fe_md %>% mutate(sample_id_in_paper = if_else(is.na(sample_id_in_paper), sample_id, sample_id_in_paper))

# reorder columns in alphabetical order so we can rbind 
oc_md <- oc_md %>% select(order(colnames(oc_md)))
fe_md <- fe_md %>% select(order(colnames(fe_md)))

# double check our column names are OK before binding
stopifnot(all(colnames(oc_md) == colnames(fe_md)))

# this is the final combined metadata
metadata <- rbind(oc_md, fe_md)

# sample type should be a factor
metadata$sample_type <- as.factor(metadata$sample_type)

metadata %>% group_by(sample_type) %>% summarize()

# reorder metadata sample types so they display in desired order
metadata$sample_type <- fct_relevel(metadata$sample_type, 
                                    "Fresh_frozen", 
                                    "Experimental_dried", 
                                    "Old_Collection", 
                                    "Positive_control", 
                                    "Negative_control")

# change category labels for plotting 
metadata$sample_type <- 
  recode(metadata$sample_type, 
         Old_Collection     = "Museum\nsamples", 
         Experimental_dried = "Experimental\ndried",
         Fresh_frozen       = "Fresh\nfrozen",
         Positive_control   = "Positive\ncontrol",
         Negative_control   = "Negative\ncontrol")


# make a sample order variable that combines collection data and weeks dried
# this will sort both old samples and new samples in chronological order
metadata <- metadata %>% mutate(sample_order = date_collected + if_else(is.na(weeks_dried), 0, weeks_dried))

# reorder sample factors for display
metadata$sample_id <- fct_reorder(metadata$sample_id, metadata$sample_order, min)
metadata$sample_id_in_paper <- fct_reorder(metadata$sample_id_in_paper, metadata$sample_order, min)

# write out the collected metadata file
write.table(metadata, paste0(output_dir,"/collected_metadata.csv"), sep=",", row.names=F, col.names=T)

# write out the collected metadata as an RDS object
saveRDS(metadata, paste0(output_dir, "/collected_metadata.rds"))
