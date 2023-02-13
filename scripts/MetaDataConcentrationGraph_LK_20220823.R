library(tidyverse)
df <- read.csv("metadata/LocationData.csv")

df2 <- sapply(df,function(x) {x <- gsub("Too Low",0,x)})
sapply(df2, class)
df2 <- transform(df2, Concentration..ng.ul. = as.numeric(Concentration..ng.ul.))
df2 <- transform(df2, Date.Collected = as.numeric(Date.Collected))
sapply(df2, class)

ggplot(df2, aes(x=Date.Collected)) + 
  geom_point(aes(y=Concentration..ng.ul., fill=Extraction.Method), shape=21, size=2, stroke=0.1, color="black") +
  scale_fill_manual(values = c("blueviolet", "darkgreen")) +
  #geom_line(aes(y=Concentration..ng.ul.) + 
  theme_bw(base_size = 11) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_rect(colour = "black", fill = "white")) +
  labs(x = "Year Collected", y = "Concentration ng/ul", fill = 'Extraction Method')
ggsave("CollectionVsConcentration.pdf", units = "in", width = 10, height = 8)
