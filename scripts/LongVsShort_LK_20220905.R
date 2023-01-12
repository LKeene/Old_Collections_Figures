library(tidyverse)
library(readxl)

df <- read_xlsx("LongVsShort_LK_20220905.xlsx")

df <- rename(.data = df,
             sample_name = Sample,
             sex = Sex,
             target = Target,
             ct = CT,
             present = Present,
             weeks = Weeks,
             replicate = Replicate)
df <- filter(.data = df, 
       present == "Y")

df <- df %>%
  group_by(target, weeks, Primer) %>%
  mutate(mean_ct = mean(ct, na.rm = TRUE),
         sd_ct = sd(ct, na.rm = TRUE),
         sample_group = paste0(target))

# Mean line bar  
ggplot(df, aes(x=weeks)) + 
  geom_point(aes(y=mean_ct, fill=Primer, shape = Primer), size=2.5, stroke=0.1, color="blue") + 
  scale_fill_manual(values = c("indianred1", "blueviolet")) + 
  scale_shape_manual(values = c(21, 22)) +
  geom_line(aes(y=mean_ct, group=Primer)) + 
  geom_errorbar(aes(ymin = (mean_ct - sd_ct), ymax = (mean_ct + sd_ct)), width = 0.1, color = "grey50", alpha = 0.5) +
  theme_bw(base_size = 11) + 
  facet_grid(~target) +
  labs(x = "Time Since Sample Collection (weeks)", y = "Mean Ct", fill = "Primer Set", shape = "Primer Set")

ggsave("Mean_LongVsShort_LK.pdf", units = "in", width = 10, height = 8)



