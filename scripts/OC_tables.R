library(tidyverse)
library(gt)
library(webshot2)

# Old Collections meta data table- Supplemental
oc_metadata <- read_csv("metadata/LocationData.csv") %>% 
  select(-c(Galbut, lat, long))

oc_metadata %>% gt(groupname_col = "Species") %>% 
  tab_style(style = list(cell_text(style = "italic")), locations = cells_row_groups()) %>% 
  tab_header(title = "Museum Collection Specimen Metadata") %>% 
  cols_align(align = "center") %>% 
  tab_footnote(footnote = "All specimens from USA unless otherwise noted", 
               locations = cells_column_labels(columns = `Location`)) %>% 
  tab_footnote(footnote = "Specimens that did not yield sufficient sequencing 
               library were not sequenced and thus have no BioSample ID", 
               locations = cells_column_labels(columns = `BioSample`)) %>% 
  tab_footnote(footnote = "RNA concentrations below the limit of detection on
               the Qubit were assigned a concentration of 0", 
               locations = cells_column_labels(columns = `Concentration (ng/μl)`)) %>% 
  tab_footnote(footnote = "Our study sample ID", 
               locations = cells_column_labels(columns = `Sample ID`)) %>%
  tab_footnote(footnote = "RNA quality was assessed using nanodrop", 
               locations = cells_column_labels(columns = `260/280`)) %>% 
  tab_footnote(footnote = "RNA quality was assessed using nanodrop", 
               locations = cells_column_labels(columns = `260/230`)) %>%
  gtsave("plots/MetaData_Supp.png", vwidth = 3000)

# Old Collections meta data table- Main
oc_metadata2 <- read_csv("metadata/LocationData.csv") %>% 
  select(-c(Galbut, lat, long, "260/280", "260/230"))

oc_metadata2 %>% gt(groupname_col = "Species") %>% 
  tab_style(style = list(cell_text(style = "italic")), locations = cells_row_groups()) %>% 
  tab_header(title = "Museum Collection Specimen Metadata") %>% 
  cols_align(align = "center") %>% 
  tab_footnote(footnote = "All specimens from USA unless otherwise noted", 
               locations = cells_column_labels(columns = `Location`)) %>% 
  tab_footnote(footnote = "Specimens that did not yield sufficient sequencing 
               library were not sequenced and thus have no BioSample ID", 
               locations = cells_column_labels(columns = `BioSample`)) %>% 
  tab_footnote(footnote = "Our study sample ID", 
               locations = cells_column_labels(columns = `Sample ID`)) %>%
  tab_footnote(footnote = "RNA concentrations below the limit of detection on
               the Qubit were assigned a concentration of 0", 
               locations = cells_column_labels(columns = `Concentration (ng/μl)`)) %>% 
  gtsave("plots/MetaData.png", vwidth = 3000)

# other viruses- known
others <- read_csv("metadata/OtherViruses_known.csv")

others <- others %>%
  select(known, year, sample_location, virus_name, percent_covered,
         region, depth, reference_accession, percent_nucleotide_similarity, sample_id, accession, rate) %>% 
  arrange(year) %>% 
  rename("Year" = "year",
         "Location" = "sample_location",
         "Taxon ID" = "virus_name",
         "Nearest GenBank Sequence" = "reference_accession",
         "% Query Coverage" = "percent_covered",
         "Average Coverage" = "depth",
         "% Identity" = "percent_nucleotide_similarity",
         "Region" = "region",
         "Sample ID" = "sample_id",
         "GenBank Sequence Accession" = "accession",
         "Rate (s/n/y)" = "rate") 

others %>% gt(groupname_col = "known") %>% 
  tab_header(title = "Virus Sequences Identified in Old Collection Specimens") %>% 
  cols_align(align = "center") %>% 
  cols_move(columns = `Sample ID`, after = `Location`) %>% 
  cols_move(columns = `% Identity`, after = `% Query Coverage`) %>% 
  cols_move(columns = `Average Coverage`, after = `Region`) %>% 
  cols_move(columns = `GenBank Sequence Accession`, after = `Taxon ID`) %>% 
  cols_move(columns = `Nearest GenBank Sequence`, after = `GenBank Sequence Accession`) %>% 
  tab_footnote(footnote = "Nearest GenBank sequence is provided for RdRp of 
               Drosophila-associated sobemo-like virus and the L and M segment of
               Puslinch virus.",
               locations = cells_column_labels(columns = `Nearest GenBank Sequence`)) %>% 
  tab_footnote(footnote = "% Query coverage was identified from the Query Cover 
               column for the nearest GenBank sequence.",
               locations = cells_column_labels(columns = `% Query Coverage`)) %>%
  tab_footnote(footnote = "% Query coverage for novel virus sequences was 
               determined based on nearest GenBank sequence.",
               locations = cells_column_labels(columns = `% Query Coverage`)) %>% 
  tab_footnote(footnote = "nt, percentage nucleotide identity calculated for 
               genome length sequences.", 
               locations = cells_column_labels(columns = `% Identity`)) %>% 
  tab_footnote(footnote = "For novel virus sequences, % nt identity to the 
               reference sequence is shown.", 
               locations = cells_column_labels(columns = `% Identity`)) %>% 
  tab_footnote(footnote = "Average mapped read coverage across contigs as 
               reported in Geneious", 
               locations = cells_column_labels(columns = `Average Coverage`)) %>% 
  tab_footnote(footnote = "Calculated using the nearest GenBank sequence with 
               date collected data", 
               locations = cells_column_labels(columns = `Rate (s/n/y)`)) %>%
  gtsave("plots/OtherViruses.png", vwidth = 1500)

