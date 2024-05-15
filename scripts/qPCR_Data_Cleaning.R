# This script uses the SangerTools r package to read in raw qPCR data and then 
# uses tidy functions to join and cleanup the data into a single tidy format for
# downstream analysis.

# Packages
library(tidyverse)
library(readxl)
library(SangerTools)
library(ggthemes)
library(writexl)

# There were some manual changes made to ensure that all the sample names followed
# the appropriate format..

# Load in data
all_data <- multiple_excel_reader("qPCR_data", pattern = "*.xls", 
                                      sheet = 3, rows_to_skip = 43, col_names = TRUE)

# Select relevant columns
all_data <- all_data %>% 
  select(sample_name, ct, tm1)

# remove control data
controls <- all_data %>% 
  filter(str_detect(sample_name, "Ctrl|cDNA|qPCR|ET|cDAN"))

# move short data
short <- all_data %>% 
  filter(str_detect(sample_name, "short")) %>% 
  filter(!str_detect(sample_name, "Ctrl|cDNA|qPCR"))

# remove misc. samples (OC/Tests/Fresh)
cleaned_noshort <- all_data %>% 
  filter(!str_detect(sample_name, paste("Ctrl|cDNA|qPCR|1919|short|",
                                        "|NoDest|ET|FoCo17|cDAN|Frozen_M|", 
                                        "|Frozen_F|EtOH|Vera|1004")))

# separate sample name into week, group, target, sex and replicate 
cleaned_noshort <- cleaned_noshort %>% 
  mutate(sample = sample_name) %>% 
  separate(col = sample_name, into = c("weekgroup", "target", "sex", "rep")) %>% 
  separate(col = weekgroup, into = c("week", "group"), sep = "wk|Wk|W")

cleaned_short <- short %>% 
  mutate(sample = sample_name) %>% 
  filter(!str_detect(sample_name, "Frozen")) %>% # short primers used on frozen 72 wk only
  separate(col = sample_name, into = c("weekgroup", "target", "length", "sex",
                                       "rep")) %>% 
  separate(col = weekgroup, into = c("week", "group"), sep = "Wk") %>% 
  mutate(target = str_replace(target, "Rpl", "RpL32"))

# separate into fly & mosquito data  
fly_data1 <- cleaned_noshort %>% 
  filter(target %in% c("Galbut", "Gal", "Galbutbut", "galbut", "Rpl", "RpL", 
                       "Thika", "LaJolla", "Nora")) %>% 
  mutate(length = "long") %>% 
  mutate(target = str_replace(target, "Rpl|RpL", "RpL32")) %>% 
  mutate(target = str_replace(target, "^Gal$|^Galbutbut$|galbut", "Galbut")) %>% 
  mutate(target = str_replace(target, "LaJolla", "La Jolla"))
# fly data w/ long primers is ready!

mosquito_data1 <- cleaned_noshort %>% 
  filter(target %in% c("Verdadero", "Act", "Actin", "Renna")) %>% 
  mutate(target = str_replace(target, "^Act$", "Actin")) %>% 
  mutate(target = str_replace(target, "Renna", "Rennavirus"))
# mosquito df is ready!

# add back short data
fly_data1 <- rbind(fly_data1, cleaned_short)
# fly df is ready!

# write both df's to manually determine presence of absence. Tm values followed:
# galbut: 84.1
# short galbut: 78.4
# rpl: 85.1
# short rpl: 80.8
# Thika: 80.67
# La Jolla: 81.1
# Nora: 81.7
# Verdadero: 83.0
# Actin: 86.4
# Rennavirus: 83.5
# samples will be determined to be present (y) by being within 1 Tm of these controls.

write_xlsx(fly_data1, "tidy_formats/fly_data1.xlsx")

write_xlsx(mosquito_data1, "tidy_formats/mosquito_data1.xlsx")
