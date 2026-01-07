library(raster)
library(terra) 
library(sf)
library(tidyverse)
# rm(list = ls())
# gc

scotland <- st_read("gridExtractions/scotlandBoundary/scotland_boundary.shp") %>%
  vect()
plot(scotland)
st_crs(scotland)

grid_clipped <- st_read("grids/grid_clipped_area2km.gpkg")
grid_clipped <- st_transform(grid_clipped, st_crs(scotland))
st_crs(grid_clipped) == st_crs(scotland)

buildings <- vect("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/scotlandOSMBuildings_Others.shp/gis_osm_buildings_a_free_1.shp")
buildings <- project(buildings, crs(scotland))
st_crs(buildings) == st_crs(grid_clipped)

grid_terra <- vect(grid_clipped) # convrt to spatvector; terra format for extraction


# Perform spatial join to count buildings in each grid cell
buildings_sf <- st_as_sf(buildings)
buildings_counts <- st_join(buildings_sf, grid_clipped, join = st_within)

# Count POIs per grid cell
grid_buildings_counts <- buildings_counts %>%
  st_drop_geometry() %>%  # Remove geometry for counting
  group_by(grid_id) %>%  
  summarise(building_count = n())

# Join the counts back to the grid
grid_clipped <- grid_clipped %>%
  left_join(grid_buildings_counts, by = "grid_id")

# Replace NA values with 0 (for grid cells with no building structures)
grid_clipped$building_count[is.na(grid_clipped$building_count)] <- 0

# grid_clipped$area_km2 <- as.numeric(st_area(grid_clipped)) / 1e6 #in km instead of m
# 
# # Calculate building density (buildings per square kilometer)
# grid_clipped$building_density <- grid_clipped$building_count / grid_clipped$area_km2


library(ggplot2)
ggplot() +
  geom_sf(data = grid_clipped, aes(fill = building_count)) +
  scale_fill_viridis_c(name = "Building Count") +
  labs(title = "Building Count by Grid Cell") +
  theme_minimal()

ggplot() +
  geom_sf(data = grid_clipped, aes(fill = building_density)) +
  scale_fill_viridis_c(
    name = "Building Count",
    option = "plasma") +
  labs(
    # title = "Posterior Mean Intensity (λ)",
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )


clean_data <- grid_clipped %>%
  st_drop_geometry() # for saving into a csv

write.csv(clean_data, "01_data/csvs/buildings_area2km.csv", row.names = FALSE)




