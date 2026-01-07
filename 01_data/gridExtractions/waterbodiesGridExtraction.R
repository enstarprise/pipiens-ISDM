#> script to extract waterbodies in grid cells

library(raster)
library(terra) 
library(sf)

# rm(list = ls())
# gc

scotland <- st_read("gridExtractions/scotlandBoundary/scotland_boundary.shp") %>%
  vect()
plot(scotland)
st_crs(scotland)

grid_clipped <- st_read("grids/grid_clipped_2point5km.gpkg")
# grid_clipped <- st_transform(grid_clipped, st_crs(scotland))
st_crs(grid_clipped)

watercourses <- st_read("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/waterbodies/oprvrs_essh_gb/data/WatercourseLink.shp")
hydronodes <- st_read("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/waterbodies/oprvrs_essh_gb/data/HydroNode.shp")
waterbodies <- st_read("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/waterbodies/GBR_wat/GBR_water_areas_dcw.shp")

