library(tidyverse)
library(readxl)

reseq <- read_xlsx("metadata/potential_reseq_2.xlsx")

reseq <- reseq %>% 
  mutate(count_type = factor(count_type, 
                                levels = c("initial", "post_trimming", 
                                           "post_collapse", "post_host_filtered")))

reseq_plot <- ggplot(reseq, aes(x = count_type, y = count)) +
  geom_point(aes(color = sample_id)) +
  scale_y_log10()

reseq_plot
