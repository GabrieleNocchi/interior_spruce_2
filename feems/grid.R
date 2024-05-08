# Working but high RAM
library(dggridR)
dggs <- dgconstruct(res=11, topology = "TRIANGLE", aperture = 4)
dgearthgrid(dggs, savegrid = "/Users/gnocc/Desktop/")




# My Solution
library(sf)
library(dplyr)
library(dggridR)

# Define the coordinates of the rectangle
area <- data.frame(
  x = c(-113.0538, -113.0538, -118.5683, -118.5683, -113.0538),  # x coordinates of vertices
  y = c(53.45184, 48.89241, 48.89241, 53.45184, 53.45184)   # y coordinates of vertices
)

# Create a simple feature (sf) object representing the rectangle
rectangle_sf <- st_polygon(list(as.matrix(area))) %>%
  st_sfc(crs = 4326)  # Assign a coordinate reference system (e.g., EPSG:4326)

# Print the object
print(rectangle_sf)


st_write(rectangle_sf, "rectangle_shapefile.shp")

library(dggridR)

dggs <- dgconstruct(topology = "TRIANGLE", res = 11, aperture = 4)
dgshptogrid(dggs, shpfname = "rectangle_shapefile.shp", savegrid = ".")
lp_sh<-st_read("dggrid.shp")
st_write(lp_sh, "final.shp", SHPT = "polygon")









## Tom Solution
dggs <- dgconstruct(topology = "TRIANGLE", res = 11, aperture = 4)
dgshptogrid(dggs, shpfname = "world-administrative-boundaries.shp", savegrid = ".")

range_wide_grid <- maptools::readShapeLines ("dggrid.shp")

# Check class of this (important for compatibility with FEEMs)
class(range_wide_grid)

lp_sh<-st_read("dggrid.shp")

spdf<- as(range_wide_grid,"SpatialPointsDataFrame")

## Grab the coords of the lodgepole samples...
lodgepole_pine_samples_all <- read.table("gab_superscaffold_heterozygosity_filtered.coord", sep = "\t")
colnames(lodgepole_pine_samples_all) <- c("LON","LAT")

lp_sample_coords <- matrix(c(lodgepole_pine_samples_all$LON, lodgepole_pine_samples_all$LAT),
                    nrow= nrow(lodgepole_pine_samples_all),
                    ncol = 2)  # closed polygon

# Extract the convex hull around the samples...
ch <- chull(lp_sample_coords)

# Grab the coords of the convex hull
ch_coords <- lp_sample_coords[ch,]

# Convert the hull to a polygon...
poly1 <- sp::Polygon(ch_coords)
hull_Poly <- sp::Polygons(list(poly1), ID = "A")

# Convert hull polygon to a SpatialPolygon...
hull_SpatialPoly <- sp::SpatialPolygons(list(hull_Poly))

# Add in cordinate metadata...
proj4string(hull_SpatialPoly) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Crop the dggrid using the spatial polygon of the sampled populations...
croppedDat <- lp_sh[as(hull_SpatialPoly, "sf"),]

# Write the cropped grid as a new shape file...
st_write(croppedDat, "final.shp", SHPT = "polygon")
