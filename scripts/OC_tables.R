# The script generates .png tables for tables in the main text of the manuscript.

library(tidyverse)
library(gt)
library(webshot2)

# Old Collections meta data table- Supplemental- full data
oc_metadata <- read_csv("metadata/LocationData.csv") %>% 
  select(-c(lat, long, "Fastq ID")) 
  

oc_metadata %>% gt(groupname_col = "Species") %>% 
  tab_style(style = list(cell_text(style = "italic")), locations = cells_row_groups()) %>% 
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
  tab_footnote(footnote = "RNA quality was assessed using a spectrophotometer", 
               locations = cells_column_labels(columns = `260/280`)) %>% 
  tab_footnote(footnote = "RNA quality was assessed using a spectrophotometer", 
               locations = cells_column_labels(columns = `260/230`)) %>%
  tab_footnote(footnote = "Galbut virus RNA 1 positive using primer numbers 1948 & 1949 (supplemental table 2)", 
               locations = cells_column_labels(columns = `Galbut virus RT-qPCR Result`)) %>% 
  gtsave("plots/Tables/MetaData_Supp.png", vwidth = 3000)

# Old Collections meta data table- Main
oc_metadata2 <- read_csv("metadata/LocationData.csv") %>% 
  select(-c("Galbut virus RT-qPCR Result", "Fastq ID", lat, long, "260/280", "260/230"))

oc_metadata2 %>% gt(groupname_col = "Species") %>% 
  tab_style(style = list(cell_text(style = "italic")), locations = cells_row_groups()) %>% 
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
  gtsave("plots/Tables/MetaData.png", vwidth = 3000)

# other viruses- known
others <- read_csv("metadata/OtherViruses_known.csv")

others <- others %>%
  select(known, year, sample_location, virus_name, percent_covered,
         region, depth, reference_accession, percent_nucleotide_similarity, sample_id, accession, rate) %>% 
  arrange(year) %>% 
  rename("Date Collected" = "year",
         "Location" = "sample_location",
         "Virus" = "virus_name",
         "Nearest GenBank" = "reference_accession",
         "% Query Coverage" = "percent_covered",
         "Average Coverage" = "depth",
         "% Identity" = "percent_nucleotide_similarity",
         "Completeness" = "region",
         "Sample ID" = "sample_id",
         "GenBank Sequence Accession" = "accession",
         "Estimated Evolutionary Rate" = "rate") 

others %>% gt(groupname_col = "known") %>% 
  cols_align(align = "center") %>% 
  cols_move(columns = `Sample ID`, after = `Location`) %>% 
  cols_move(columns = `% Identity`, after = `% Query Coverage`) %>% 
  cols_move(columns = `Average Coverage`, after = `Completeness`) %>% 
  cols_move(columns = `GenBank Sequence Accession`, after = `Virus`) %>% 
  cols_move(columns = `Nearest GenBank`, after = `GenBank Sequence Accession`) %>% 
  tab_footnote(footnote = "Nearest GenBank sequence is provided for RdRp of 
               Drosophila-associated sobemo-like virus and the L, M and S segment of
               Puslinch virus.",
               locations = cells_column_labels(columns = `Nearest GenBank`)) %>% 
  tab_footnote(footnote = "% Query coverage from BLASTN alignment",
               locations = cells_column_labels(columns = `% Query Coverage`)) %>%
  tab_footnote(footnote = "% Query coverage from BLASTN alignment for novel virus sequences was 
               determined based on nearest GenBank sequence.",
               locations = cells_column_labels(columns = `% Query Coverage`)) %>% 
  tab_footnote(footnote = "nt, percentage nucleotide identity to closest GenBank sequence identified via BLASTn.", 
               locations = cells_column_labels(columns = `% Identity`)) %>% 
  tab_footnote(footnote = "For novel virus sequences, % nt identity to the 
               reference sequence is shown.", 
               locations = cells_column_labels(columns = `% Identity`)) %>% 
  tab_footnote(footnote = "Average mapped read coverage across contigs as 
               reported in Geneious.", 
               locations = cells_column_labels(columns = `Average Coverage`)) %>% 
  tab_footnote(footnote = "Calculated using the nearest GenBank sequence with 
               date collected data.", 
               locations = cells_column_labels(columns = `Estimated Evolutionary Rate`)) %>%
  gtsave("plots/Tables/OtherViruses.png", vwidth = 1500)

# Old Collections Dates- Supplemental Methods
oc_dates <- read_csv("metadata/Extraction_Prep_Dates_Table.csv")


oc_dates %>% gt() %>% 
  tab_style(style = list(cell_text(style = "italic")), locations = cells_row_groups()) %>% 
  cols_align(align = "center") %>% 
  tab_footnote(footnote = "Select specimens sequenced a second time to increase
               coverage of virus sequecnes", 
               locations = cells_column_labels(columns = `Secondary Library Preparation Date`)) %>% 
  tab_footnote(footnote = "Select specimens sequenced a second or third time to
               increase coverage of virus sequences", 
               locations = cells_column_labels(columns = `Tertiary Library Preparation Date`)) %>% 
  gtsave("plots/Tables/Date_Supp.png", vwidth = 3000)

