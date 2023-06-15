library(tidyverse)

matrix <- read.delim("metadata/OC_All_taxa_matrix.txt", sep = "\t", skip = 1, row.names = 1)

matrix <- matrix %>% 
  slice(2:38)

matrix <- t(matrix) 
matrix <- as.data.frame(matrix)

matrix <- rownames_to_column(matrix, var = "organisms")

# filter down to viruses
attempt <- matrix %>% 
  filter(str_detect(organisms, "virus|phage|viridae|virales|viriform|viria|env59"))
