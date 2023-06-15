library(tidyverse)
library(readxl)
library(gt)
library(webshot2)

# Old Collections meta data table
oc_metadata <- read_csv("metadata/LocationData.csv")

oc_metadata %>% gt(groupname_col = "Species") %>% 
  tab_header(title = "Museum Collection Specimen Metadata") %>% 
  cols_align(align = "center") %>% 
  gtsave("plots/MetaData.png", vwidth = 1500)


# other viruses- known
others <- read_csv("metadata/OtherViruses_known.csv")

others <- others %>% 
  select(year, Drosophila_species, sample_location, virus_name, percent_covered, depth,
         percent_nucleotide_similarity) %>% 
  arrange(year) %>% 
  rename("Year" = "year",
         "Species" = "Drosophila_species",
         "Location" = "sample_location",
         "Virus" = "virus_name",
         "% Coverage" = "percent_covered",
         "Depth" = "depth",
         "% Similarity to Reference" = "percent_nucleotide_similarity") 

others %>% gt(groupname_col = "Species") %>% 
  tab_header(title = "Other Known Virus Sequences") %>% 
  cols_align(align = "center") %>% 
  gtsave("plots/known_otherviruses.png", vwidth = 1500)

# other viruses- new

