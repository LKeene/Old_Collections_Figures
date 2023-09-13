library(tidyverse)
library(gt)
library(webshot2)

# Old Collections meta data table
oc_metadata <- read_csv("metadata/LocationData.csv") %>% 
  select(-c(long, lat))

oc_metadata %>% gt(groupname_col = "Species") %>% 
  tab_header(title = "Museum Collection Specimen Metadata") %>% 
  cols_align(align = "center") %>% 
  gtsave("plots/MetaData.png", vwidth = 3000)


# other viruses- known
others <- read_csv("metadata/OtherViruses_known.csv")

others <- others %>%
  select(year, Drosophila_species, sample_location, taxon_id, percent_coverage,
         region, average_coverage, reference_accession, percent_nucleotide_identity) %>% 
  arrange(year) %>% 
  rename("Year" = "year",
         "Species" = "Drosophila_species",
         "Location" = "sample_location",
         "Taxon ID" = "taxon_id",
         "Nearest Genbank Sequence" = "reference_accession",
         "% Query Coverage" = "percent_coverage",
         "Average Coverage" = "average_coverage",
         "% Identity" = "percent_nucleotide_identity",
         "Region" = "region") 

others %>% gt(groupname_col = "Species") %>% 
  tab_header(title = "Other Known Virus Sequences") %>% 
  cols_align(align = "center") %>% 
  cols_move(columns = `Nearest Genbank Sequence`, after = `Taxon ID`) %>% 
  cols_move(columns = `% Identity`, after = `% Query Coverage`) %>% 
  cols_move(columns = `Average Coverage`, after = Region) %>% 
  tab_footnote(footnote = "nt, percentage nucleotide identity calculated for 
               genome length sequences.", 
               locations = cells_column_labels(columns = `% Identity`)) %>% 
  tab_footnote(footnote = "Average mapped read coverage across contigs as 
               reported in Geneious", 
               locations = cells_column_labels(columns = `Average Coverage`)) %>% 
  gtsave("plots/known_otherviruses.png", vwidth = 1500)


# other viruses- new

