
library(geojsonsf)
library(ggplot2)
library(tidyverse)
library(sf)
library(terra)

# windSpeed <- geojson_sf("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/windSpeed/windspeeds/windspeeds.geojson")
# st_crs(windSpeed)

scotland <- st_read("gridExtractions/scotlandBoundary/scotland_boundary.shp") %>%
  vect()
plot(scotland)
st_crs(scotland)

windSpeed <-  st_read("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/windSpeed/windspeeds/windspeeds.geojson") |> 
  st_transform(st_crs(scotland)) |>
  vect()

grid_clipped <- st_read("grids/grid_clipped.gpkg")
grid_clipped <- st_transform(grid_clipped, st_crs(scotland))
st_crs(grid_clipped)
st_crs(windSpeed) == st_crs(grid_clipped)

# crop and mask and whatever
windSpeed_cropped <- crop(windSpeed, scotland)
windSpeed_masked <- mask(windSpeed_cropped, scotland_vect)

# 4. Convert back to sf if needed
windSpeed_final <- st_as_sf(windSpeed_masked)

# 5. Visual verification
plot(windSpeed_masked, "windspeed")  # Replace with actual attribute name
plot(scotland_vect, add=TRUE, border="red")
# ----- extract mean values to grid -----
grid_5km_windSpeed <- grid_clipped
grid_5km_livestock$livestock_density <- terra::extract(
  livestock_5km, #livestock_5km_filled
  grid_terra, 
  fun = mean,
  na.rm = TRUE
)[, 2]
