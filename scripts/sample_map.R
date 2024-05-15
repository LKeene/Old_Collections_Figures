library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(dplyr)
library(readxl)
library(gridExtra)
library(ggpubr)
library(ggthemes)
library(viridis)

locations <- read.csv("metadata/LocationData.csv")

# Get world map from map_data package
usa <- map_data("state")
canada <- map_data("worldHires", "canada")
mexico <- map_data("worldHires", "mexico")

# geom_polygon -> makes base of map 
# geom_point -> create points colored by host
# Add geom_points (latitude and longitude) 
# do 2 geom_point calls: 1 for colored circles and one for black outlines
# coord_map clips the map 
map <- ggplot() + 
 geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = "snow2",
              color = "grey70", size=0.2) + 
 geom_polygon(data = canada, aes(x = long, y = lat, group = group), 
              fill = "snow3", color = "grey70", size = 0.2) +
 geom_polygon(data = mexico, aes(x = long, y = lat, group = group), 
              fill = "snow4", color = "grey70", size = 0.2) +
 coord_fixed(xlim = c(-160, -60), ylim = c(20,49), ratio = 1.2) +
 theme_few() +
 scale_fill_viridis() +
 geom_jitter(data = locations, aes(x = long, y = lat, fill = Date.Collected), 
             shape = 21, stroke = 0.25, size = 4, alpha = 0.6, width = 0.5,
             height = 0.5) +
 labs(fill = "Date Collected") +
 theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
      legend.text=element_text(family="Helvetica", size=10),
      axis.ticks=element_blank(),
      axis.title.y=element_blank(),
      axis.title.x=element_blank(),
      text = element_text(size = 20))
                                                    

map <- map + geom_text() +
  annotate("text", label = "Hawaii", x = -158, y = 24.5)
map
ggsave("plots/Figure2/location_map.pdf", units = "in", width = 10, height = 3)
ggsave("plots/Figure2/location_map.jpg", units = "in", width = 10, height = 3)







