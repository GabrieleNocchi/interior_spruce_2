library(ggplot2)
library(cowplot)
library(sp)
library(raster)

data <- read.table("metadata_coord.txt", sep ="\t")
data_name <- cbind(paste(data$V3,data$V4),data$V4,as.numeric(data$V6),as.numeric(data$V5))
data_name <- unique(data_name)
data_name <- as.data.frame(data_name)
data <- cbind(as.numeric(data_name$V3), as.numeric(data_name$V4))
data <- as.data.frame(data)
colnames(data) <- c( "Long", "Lat")
coordinates(data) <- c("Long", "Lat")
proj4string(data) <- CRS("+proj=longlat +datum=WGS84")
files <- list.files(pattern = "\\.tif$")
files <- files[!grepl("wc2.1_30s_elev\\.tif", files)]
env_data <- stack(files)
my_env <- extract(env_data,data)
my_env <- as.data.frame(my_env)
colnames(data_name) <- c( "Site","Topography","Long", "Lat")
my_env <- cbind(data_name, my_env)
write.table(my_env, file = "env.txt", sep = "\t", quote = F, row.names = F)

# Assuming your data is stored in a variable called 'final'
final <- read.table("env.txt",sep = "\t", h = T)





# First ordering

final <- final[order(final$Site), ]

# Create a list to store plots
plots <- list()

for (i in colnames(final[, 5:8])) {
 # Create the bar plot for each variable
 plot <- ggplot(final, aes(x = Site, y = !!sym(i))) +
 geom_bar(stat = "identity", fill = c("turquoise","royalblue","blue","turquoise","royalblue","blue","turquoise","royalblue","blue","turquoise","royalblue","blue","turquoise","royalblue","blue","turquoise","royalblue","blue","turquoise","royalblue","blue","turquoise","royalblue","blue")) +
 labs(y = i) +
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

 # Add the plot to the list
 plots[[i]] <- plot
}

for (i in colnames(final[, 9:12])) {
 # Create the bar plot for each variable
 plot <- ggplot(final, aes(x = Site, y = !!sym(i))) +
 geom_bar(stat = "identity", fill = c("red3","orange","gold","red3","orange","gold","red3","orange","gold","red3","orange","gold","red3","orange","gold","red3","orange","gold","red3","orange","gold","red3","orange","gold")) +
 labs(y = i) +
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

 # Add the plot to the list
 plots[[i]] <- plot
}

for (i in colnames(final[, 13:16])) {
 # Create the bar plot for each variable
 plot <- ggplot(final, aes(x = Site, y = !!sym(i))) +
 geom_bar(stat = "identity", fill = c("lightcyan","grey70", "grey16","lightcyan","grey70", "grey16","lightcyan","grey70", "grey16","lightcyan","grey70", "grey16","lightcyan","grey70", "grey16","lightcyan","grey70", "grey16","lightcyan","grey70", "grey16","lightcyan","grey70", "grey16")) +
 labs(y = i) +
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

 # Add the plot to the list
 plots[[i]] <- plot
}

# Arrange the plots in a grid with 4 plots per row and 3 rows
combined_plot <- plot_grid(plotlist = plots, ncol = 4, nrow = 3)

# Print the combined plot
print(combined_plot)





# Second ordering
final <- final[order(final$Topography,final$Site), ]
final$Site <- factor(final$Site,levels = unique(final$Site))
# Create a list to store plots
plots <- list()

for (i in colnames(final[, 5:8])) {
 # Create the bar plot for each variable
 plot <- ggplot(final, aes(x = Site, y = !!sym(i))) +
 geom_bar(stat = "identity", fill = c("turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","royalblue","royalblue","royalblue","royalblue","royalblue","royalblue","royalblue","royalblue","blue","blue","blue","blue","blue","blue","blue","blue")) +
 labs(y = i) +
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

 # Add the plot to the list
 plots[[i]] <- plot
}

for (i in colnames(final[, 9:12])) {
 # Create the bar plot for each variable
 plot <- ggplot(final, aes(x = Site, y = !!sym(i))) +
 geom_bar(stat = "identity", fill = c("red3","red3","red3","red3","red3","red3","red3","red3","orange","orange","orange","orange","orange","orange","orange","orange","gold","gold","gold","gold","gold","gold","gold","gold")) +
 labs(y = i) +
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

 # Add the plot to the list
 plots[[i]] <- plot
}

for (i in colnames(final[, 13:16])) {
 # Create the bar plot for each variable
 plot <- ggplot(final, aes(x = Site, y = !!sym(i))) +
 geom_bar(stat = "identity", fill = c("lightcyan","lightcyan","lightcyan","lightcyan","lightcyan","lightcyan","lightcyan","lightcyan","grey70", "grey70", "grey70", "grey70", "grey70", "grey70", "grey70", "grey70", "grey16","grey16","grey16","grey16","grey16","grey16","grey16","grey16")) +
 labs(y = i) +
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

 # Add the plot to the list
 plots[[i]] <- plot
}

# Arrange the plots in a grid with 4 plots per row and 3 rows
combined_plot <- plot_grid(plotlist = plots, ncol = 4, nrow = 3)

# Print the combined plot
print(combined_plot)









# MAP
library(ggplot2)
library(raster)
library(rasterVis)
library(rgdal)
library(grid)
library(scales)
library(viridis)
library(ggthemes)
library(ggrepel)
library(raster)
library(marmap)


coo <- read.table("env.txt", h = T, sep = "\t")

# Raster file altitude
alt <- raster::raster("./wc2.1_30s_elev.tif")

# Load shapefile of Ugandan waterbodies
wat <-  rgdal::readOGR("./Base Waterbody Polygon.shp")
# roads <-  rgdal::readOGR("./alberta_highway.shp")
# roads_2 <-  rgdal::readOGR("./british_columbia_highway.shp")
border <-  rgdal::readOGR("./CAN_adm1.shp")

ll <- c(51.42,-116.17, "Lake Louise")
canmore <-c(51.08,-115.34, "Canmore")
banff <- c(51.17,-115.57,"Banff")
kana <- c(50.91,-115.14,"Kananaskis Village")

cities <- rbind(ll,canmore,banff, kana)
cities <- as.data.frame(cities)
colnames(cities) <- c("Lat","Long","Site")

cities$Long <- as.numeric(cities$Long)
cities$Lat <- as.numeric(cities$Lat)

coo_df <- rbind(cities, coo[, c("Lat", "Long", "Site")])

# Assuming 'alt' and 'wat' are your raster objects
ext <- extent(
  min(coo$Long) - 0.1,
  max(coo$Long) + 0.1,
  min(coo$Lat) - 0.1,
  max(coo$Lat) + 0.1
)



# Crop 'alt' and 'wat' to the defined extent
alt_cropped <- crop(alt, ext)
wat_cropped <- crop(wat, ext)
# roads_cropped <- crop(roads,ext)
# roads_cropped_2 <- crop(roads_2,ext)
border_cropped <- crop(border,ext)


# GGPLOT
test_spdf <- as(alt_cropped, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")


brown_palette <- c("darkolivegreen4","darkgreen","wheat2","sienna", "#CD853F","#8B4513", "#654321")

ggplot() +
geom_tile(data = test_df, aes(x = x, y = y, fill = value), alpha = 0.8) +
geom_polygon(data = wat_cropped, aes(x = long, y = lat, group = group),color = "#00a0a0",fill = "azure", size = 0.2) +
geom_polygon(data = border_cropped, aes(x = long, y = lat, group = group),color ="black" ,fill = "NA", size = 0.75) +
#geom_polygon(data = roads_cropped, aes(x = long, y = lat, group = group),color = "white", fill = "NA", size = 0.1) +
#geom_polygon(data = roads_cropped_2, aes(x = long, y = lat, group = group),color = "white", fill = "NA", size = 0.1) +
scale_fill_gradientn(colors = brown_palette) +
# scale_fill_etopo() +
coord_equal() +
theme_map() +
theme(legend.position = "bottom", legend.key.width = unit(0.5, "cm")) +

geom_point(data = coo_df, aes(x = Long, y = Lat),color = "black", fill = c("black","black","black", "black","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red"
), size = c(3,3,3,3,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5), shape = c(22,22,22,22,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21
)) +
geom_label_repel(data = coo_df, aes(x = Long, y = Lat, label = Site), size = c(3.5,3.5,3.5,3.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5), color = "black", fill = c("white","white","white","white","lightpink","lightpink","lightpink","lightyellow","lightyellow","lightyellow","lightblue","lightblue","lightblue","darkolivegreen3","darkolivegreen3","darkolivegreen3","sienna1","sienna1","sienna1","lightseagreen","lightseagreen","plum1","plum1","lightgoldenrod","plum1","lightgoldenrod","lightgoldenrod","lightseagreen"

)) +
ggtitle("Spruce Sampling") +
ggspatial::annotation_scale(location = "br",pad_x = unit(1.2, "cm"),pad_y = unit(0.1, "cm")) +
ggspatial::annotation_north_arrow(location = "tr", which_north = "true", pad_x = unit(1.8, "cm"), pad_y = unit(1, "cm"),style = ggspatial::north_arrow_fancy_orienteering(fill = c("black", "white"), line_col = "black")) +
theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 30), legend.position = c(0.07,0.1),legend.box.background = element_rect(color = "black", fill = NA)) + labs(fill = "Elevation (m)")
