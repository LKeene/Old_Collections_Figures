library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
# library(plyr)
library(dplyr)
# library(googlesheets)
library(readxl)
library(gridExtra)
library(ggpubr)
library(RColorBrewer)

# can read from google sheets?

# gs_auth()
# listing<-gs_ls()
# listing$sheet_title

# open google spreadsheet
# bunya_gs <- gs_title("Bunyavirus genome assembly status")
# gs_ws_ls(bunya_gs)
# bunya_meta_data <- gs_read(ss=bunya_gs, ws = "Virus_metadata")
# bunya_meta_data

bunya_meta_data <- read_excel("./Bunyavirus_metadata.xlsx")

bunya_meta_data %>% filter(Abbreviation == 'PICV') %>% select(Host_category)

# convert empty values in table to NA values for consistency
# make a vector of the Host_Name colum
# convert empty -> NA
is.na(bunya_meta_data$Host_category) <- bunya_meta_data$Host_category == ""


# Get world map from map_data package
world <- map_data("world")

# this defines limits of the world map in long/lat that will be plotted
xlim = c(-180,180)
ylim = c(-55,75)


# Create a caegorical color scale for hosts
# see: http://www.datavis.ca/sasmac/brewerpal.html
# myColors <- brewer.pal(8,"Dark2")
# myColors <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
myColors <- brewer.pal(8,"Set1")
# names(myColors) <- levels(bunya_meta_data$Host_category)
fillScale <- scale_fill_manual(values = myColors, na.value="black")
# unique(bunya_meta_data$Host_category)

# geom_polygon -> makes base of map 
# geom_point -> create points colored by host
# Add geom_points (latitude and longitude) 
# do 2 geom_point calls: 1 for colored circles and one for black outlines
# coord_map clips the map 
world_map <- ggplot() + 
 geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "grey90", color = "grey70", size=0.2) + theme_bw() +
 # geom_point(data = filter(bunya_meta_data, Virus=="Caimito"), aes(x = Longitude, y = Latitude, fill=Host_category), color="black", shape=21, size=1.5, alpha=0.90, stroke=0.2) +
 geom_point(data = bunya_meta_data, aes(x = Longitude, y = Latitude, fill=Host_category), color="black", shape=21, size=1.5, alpha=0.90, stroke=0.2) +
 # geom_point(data = bunya_meta_data, aes(x = Longitude, y = Latitude), size=1.5, color="black", shape = 21, alpha=0.35, stroke=0.2) +
 coord_map(xlim=xlim, ylim=ylim) + fillScale + theme(axis.text.x=element_blank(),
                                                    axis.text.y=element_blank(),
                                                    legend.text=element_text(family="Helvetica", size=10),
                                                    axis.ticks=element_blank(),
                                                    axis.title.y=element_blank(),
                                                    axis.title.x=element_blank())
                                                    


world_map

# how many countries? 
countries <- unique(bunya_meta_data$country)



# make plot of samples by years
years <- as.numeric(format(as.Date(bunya_meta_data$`collection-date`, origin = "1970-01-01", format="%m/%d/%Y"),"%Y"))
bunya_meta_data$years <- years

newest_sample <- max(years, na.rm = TRUE)
oldest_sample <- min(years, na.rm = TRUE)
year_span <- newest_sample - oldest_sample + 1

# extract decades from years
decades <- sapply(years, function(x) x - x %% 10 )
decades
bunya_meta_data$decades <- as.factor(decades)

# plot years as a histogram
year_plot <- ggplot(bunya_meta_data, aes(years)) + 
  geom_histogram(bins=year_span, fill="slategray", color="black", size=0.2) +
  labs(y = "Isolates") +
  scale_x_continuous(breaks=seq(oldest_sample, newest_sample, by = 5)) + 
  scale_y_continuous(breaks=c(2,4,6,8,10,12)) +
  
  theme(axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text.x = element_text(family="Helvetica", size=8, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(family="Helvetica", size=8),
        axis.title.y = element_text(family="Helvetica", size=9))
#          axis.line = element_line(colour = "grey"))
year_plot

# plot countries
countries <- bunya_meta_data$country

country_counts <- count(bunya_meta_data, country)
country_counts <- country_counts[complete.cases(country_counts), ]
country_counts


country_plot <- ggplot(filter(country_counts, n>1), aes(x=reorder(country, -n), y=n)) + 
  geom_bar(stat="identity", fill="slategray", color="black", size=0.2) +
  labs(y = "Isolates") +
  # scale_x_continuous(breaks=seq(oldest_sample, newest_sample, by = 5)) + 
  # scale_y_continuous(breaks=c(2,4,6,8,10,12)) +
  theme(axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text.x = element_text(family="Helvetica", size=8, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(family="Helvetica", size=8),
        axis.title.y = element_text(family="Helvetica", size=9))
country_plot

# plot host freq chart
host_counts <- count(bunya_meta_data, host_category_medium)
# host_counts <- host_counts[complete.cases(host_counts), ]
host_colors <- select(bunya_meta_data, host_category_medium, Host_category)
host_counts <- inner_join(host_counts, host_colors)
host_counts <- unique(host_counts)

genus_counts <- count(bunya_meta_data, Host_genus)


# host_plot <- ggplot(filter(host_counts, !is.na(host_category_medium) & host_category_medium != "Not known"), aes(x=reorder(host_category_medium, -n), y=n, fill=Host_category)) + 
host_plot <- ggplot(filter(host_counts, !is.na(host_category_medium)), aes(x=reorder(host_category_medium, -n), y=n, fill=Host_category)) + 
  # geom_bar(stat="identity", color="black", size=0.2) +
  geom_col(color="black", size=0.2) +
  fillScale + 
  labs(y = "Isolates") +
  # scale_x_continuous(breaks=seq(oldest_sample, newest_sample, by = 5)) + 
  # scale_y_continuous(breaks=c(2,4,6,8,10,12)) +
  theme(axis.title.x=element_blank(),
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text.x = element_text(family="Helvetica", size=8, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(family="Helvetica", size=8),
        axis.title.y = element_text(family="Helvetica", size=9))
host_plot

# layout as a multi-panel plot
# do this in 3 rows because:
# - 2nd row will have aligned plots but with no labels (labels mess up alignment)
# - 3rd row will have unaligned plots, but will have labels
# will have to merge 2nd & 3rd rows in illustrator to make a final version with aligned axes & labels
ggarrange(world_map,                                         # First row with scatter plot
          ggarrange(year_plot + theme(axis.text.x=element_blank()), country_plot + theme(axis.text.x=element_blank()), host_plot  + theme(axis.text.x=element_blank()), 
                    ncol = 3, 
                    labels = c("B", "C", "D"), 
                    font.label = list(size = 10, face = "plain"),
                    align="v"), # Second row with box and dot plots
          ggarrange(year_plot, country_plot, host_plot, 
                    ncol = 3, 
                    labels = c("B", "C", "D"), 
                    font.label = list(size = 10, face = "plain"),
                    align="v"), # Second row with box and dot plots
          nrow = 3, 
          font.label = list(size = 10, face = "plain"),
          labels = "A"                                        # Labels of the scatter plot
) 

ggsave("fig_1_extra_row.pdf", width=7, height=7, units="in")

# create a decades color scale
myColors <- brewer.pal(8,"Reds")
names(myColors) <- levels(as.factor(decades))
decadesColScale <- scale_colour_manual(name = "Decade of original isolation",values = myColors, na.value="black")
gg1 + geom_point(data = plot_data, aes(x = long, y = lat, color = plot_data$decades), size=3, alpha=0.9) +
  geom_point(data = plot_data, aes(x = long, y = lat), size=3, color="black", shape = 21, alpha=0.25, stroke=0.25) +
  coord_map(xlim=xlim, ylim=ylim) + decadesColScale






