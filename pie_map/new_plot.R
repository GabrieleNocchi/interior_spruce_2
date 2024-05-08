library(mapplots)
library(rworldmap)
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

brown_palette <- c("darkolivegreen4","wheat2", "#CD853F","#8B4513", "#654321")

grey_palette <- gray.colors(10, start = 0.9, end = 0.1)

data(countriesLow)
plot(countriesLow)
coo <- read.table("ordered_metadata_with_coord.txt", h = F, stringsAsFactors = F, row.names = NULL, sep = "\t")
my_dir <- getwd()
file <- list.files(path=my_dir, pattern="*meanQ", full.names=TRUE, recursive=FALSE)
a <- read.table(file)
gabri <- cbind(a$V1,a$V2)
gabri <- as.numeric(gabri)
n <- length(a$V1)
gabri_2 <- c(rep("pop1",n), rep("pop2",n))
gabri_2 <- as.factor(gabri_2)
lon <- rep(coo$V5,2)
lat <- rep(coo$V6,2)
# added extra digit to duplicate coordinates to avoid removal of records
xyz <- make.xyz(lon, lat, gabri,gabri_2)
#pdf("map.pdf")
plot(countriesLow, xlim= c(min(na.omit(coo$V5)-0.05),max(na.omit(coo$V5) + 0.05)), ylim = c(min(na.omit(coo$V6)-0.05),max(na.omit(coo$V6)+0.05)))
COL <- c("darkgreen", "lightgrey")

# Raster file altitude
alt <- raster::raster("./wc2.1_30s_elev.tif")

# Load shapefile of Ugandan waterbodies
wat <-  rgdal::readOGR("./Base Waterbody Polygon.shp")

border <-  rgdal::readOGR("./CAN_adm1.shp")



border <- crop(border, extent(c(min(na.omit(coo$V5) - 0.05), max(na.omit(coo$V5) + 0.05),
                                min(na.omit(coo$V6) - 0.05), max(na.omit(coo$V6) + 0.05))))


wat <- crop(wat, extent(c(min(na.omit(coo$V5) - 0.05), max(na.omit(coo$V5) + 0.05),
                                min(na.omit(coo$V6) - 0.05), max(na.omit(coo$V6) + 0.05))))
								
alt <- crop(alt, extent(c(min(na.omit(coo$V5) - 0.05), max(na.omit(coo$V5) + 0.05),
                                min(na.omit(coo$V6) - 0.05), max(na.omit(coo$V6) + 0.05))))




ll <- c(51.42,-116.16, "Lake Louise")
canmore <-c(51.08,-115.34, "Canmore")
banff <- c(51.17,-115.51,"Banff")
kana <- c(50.91,-115.14,"Kananaskis Village")
emerald <- c(51.455322,-116.540603, "1. Emerald")
cathedral <- c(51.445664,-116.327133, "2. Cathedral")
bell <- c(51.305969,-116.024544, "3. Bell")
norquay <- c(51.204481,-115.667814, "4. Norquay")
cascade <- c(51.248633,-115.554167, "5. Cascade")
canmore_2 <- c(51.050556,-115.329208, "6. Canmore")
oldbaldy <- c(50.842411,-115.131753, "8. Old Baldy")
yates <- c(51.038728,-115.053825, "7. Yates")


plot(alt, add = TRUE, col = grey_palette)
plot(wat, add = TRUE, col ="azure", border = "#00a0a0")
plot(border, add = TRUE, lwd = 2)



points_df_2 <- data.frame(
  lat = as.numeric(c(emerald[1], cathedral[1], bell[1], norquay[1], cascade[1], canmore_2[1], oldbaldy[1], yates[1])),
  lon = as.numeric(c(emerald[2], cathedral[2], bell[2], norquay[2], cascade[2], canmore_2[2], oldbaldy[2], yates[2])),
  name = c(emerald[3], cathedral[3], bell[3], norquay[3], cascade[3], canmore_2[3], oldbaldy[3], yates[3])
)


# Create a data frame with points information
points_df <- data.frame(
  lat = as.numeric(c(ll[1],  banff[1], kana[1])),
  lon = as.numeric(c(ll[2],  banff[2], kana[2])),
  name = c(ll[3], banff[3], kana[3])
)

# Add points to the plot
#points(points_df$lon, points_df$lat, pch = 16, col = "red")

# Add text labels for the points with an offset

# Add a white square background around text labels
# Add a white square background around text labels
#for (i in 1:nrow(points_df)) {
#  text_label <- points_df$name[i]
#  text_x <- points_df$lon[i]
#  text_y <- points_df$lat[i]
#  text_width <- strwidth(text_label, cex = 1, font = 2)
#  text_height <- strheight(text_label, cex = 1, font = 2)
  
 # rect(
 #   text_x - text_width / 2 - 0.005 ,  # Adjust the padding as needed
  #  text_y + text_height / 2 + 0.035,  # Adjust the padding as needed
   # text_x + text_width / 2  + 0.005,  # Adjust the padding as needed
   # text_y - text_height / 2 + 0.025,  # Adjust the padding as needed
   # col = "white", border = "black"
 # )
#}


for (i in 1:nrow(points_df_2)) {
  text_label <- points_df_2$name[i]
  text_x <- points_df_2$lon[i]
  text_y <- points_df_2$lat[i]
  text_width <- strwidth(text_label, cex = 1, font = 2)
  text_height <- strheight(text_label, cex = 1, font = 2)
  
  rect(
    text_x - text_width / 2 - 0.005 ,  # Adjust the padding as needed
    text_y + text_height / 2 + 0.035,  # Adjust the padding as needed
    text_x + text_width / 2  + 0.005,  # Adjust the padding as needed
    text_y - text_height / 2 + 0.025,  # Adjust the padding as needed
    col = "white", border = "black"
  )
}

text(points_df_2$lon, points_df_2$lat, labels = points_df_2$name, pos = 3, offset = 1, col = "black", cex = 0.8, font = 2)

# Add text labels for the points with an offset
text(points_df$lon, points_df$lat, labels = points_df$name, pos = 3, offset = 1, col = "black", cex = 1, font = 2)


draw.pie(xyz$x, xyz$y, xyz$z, radius = 0.01, col = COL)
