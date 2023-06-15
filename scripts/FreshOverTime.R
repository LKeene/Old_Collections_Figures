library(tidyverse)
library(readxl)
library(ggthemes)

fly_data3 <- read_xlsx("tidy_formats/fly_data3.xlsx")
fresh <- fly_data3 %>% 
  filter(group == "Fresh") %>% 
  filter(target %in% c("Galbut", "RpL32")) 

ggplot(fresh, aes(x = week)) + 
  geom_smooth(aes(y = mean_ct, group = sample_length), alpha = 0.2, 
              colour = "darkgrey", span = 0.3, linewidth = 0.5)+
  geom_point(aes(y = ct, fill = sex), shape = 21, size = 2, stroke = 0.1, 
             color = "black", alpha = 0.75) + 
  scale_fill_manual(values = c("blue", "blueviolet")) + 
  theme_few(base_size = 11) + 
  facet_grid(.~target, scales = "free_y") +
  labs(x = "Time (weeks)", y = "Mean Ct", fill = "Sex")

ggsave("plots/fresh.pdf", units = "in", width = 10, height = 8)  
  