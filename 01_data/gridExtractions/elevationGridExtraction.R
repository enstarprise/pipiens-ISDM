## script to extract elevation data

library(raster)
library(terra) 
library(sf)

# rm(list = ls())
# gc

scotland <- st_read("gridExtractions/scotlandBoundary/scotland_boundary.shp") %>%
  vect()
plot(scotland)
st_crs(scotland)

grid_clipped <- st_read("grids/grid_clipped_1km.gpkg")
grid_clipped <- st_transform(grid_clipped, st_crs(scotland))
st_crs(grid_clipped)

elevation <- rast("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/elevation_DEM/scotland_elevation.tif")
elevation <- project(elevation, crs(scotland))
st_crs(elevation) == st_crs(grid_clipped)

plot(elevation)

grid_terra <- vect(grid_clipped) # convrt to spatvector; terra format for extraction

# ----- extract values to 5km grid -----
grid_1km_elevation <- grid_clipped
grid_1km_elevation$elevation <- terra::extract(
  elevation, 
  grid_terra, 
  fun = mean,
  na.rm = TRUE
)[, 2] # will take a while 

library(tmap)
tm_shape(grid_1km_elevation) +
  tm_polygons("elevation",
              style = "quantile",
              # palette = "RdYlBu",
              title = "Elevation (m)") +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(legend.outside = TRUE)

write_csv(grid_1km_elevation, "csvs/elevation_1km.csv")








########### ISLE OF ARRAN 
elevation <- rast("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/elevation_DEM/scotland_elevation.tif")
elevation <- project(elevation, crs(scotland))
st_crs(elevation) == st_crs(grid_clipped)

plot(elevation)

mask <- st_read("~/OneDrive - University of Glasgow/PhD/modelling/dataIntegration/MARS/arran_mask.shp")

# Convert sf object to SpatVector for terra compatibility
mask_vect <- vect(mask)

# Crop elevation to the extent of the mask (rectangular bounding box)
elevation_cropped <- crop(elevation, mask_vect)

# Mask to the exact shape (sets values outside polygon to NA)
elevation_masked <- mask(elevation_cropped, mask_vect)

plot(elevation_masked)





