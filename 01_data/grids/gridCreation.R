library(sf)
library(future.apply)
library(furrr)
plan(multisession)  # For parallel processing

scotland <- st_read("gridExtractions/scotlandBoundary/scotland_boundary.shp")
scotland <- st_transform(scotland,crs = 27700)

####### 2 km area
grid <- st_make_grid(scotland, cellsize = 1400, square = TRUE)
grid <- st_sf(grid_id = 1:length(grid), geometry = grid)
st_write(grid, "grids/grid_area2km.gpkg")

plan(multisession, workers = 10)  # use 10 cores # plan(multisession, workers = availableCores() - 1)
n_chunks <- 50  # 40/10 = 4 chunks of the dataset per core
grid_split <- split(grid, cut(1:nrow(grid), n_chunks, labels = FALSE))

# parallel intersection
grid_clipped_list <- future_map(
  grid_split,
  ~ st_intersection(.x, scotland),
  .options = furrr_options(seed = TRUE)
) 

grid_clipped_area2km <- do.call(rbind, grid_clipped_list) # combine results
plan(sequential) # reset to sequential
st_write(grid_clipped_area2km, "grids/grid_clipped_area2km.gpkg")

max(st_area(grid_clipped_area2km))
max(st_area(grid_clipped_area2km)/1e6) # 1.96 km area


# save as CSV 
grid_clipped_area2km <- st_read("grids/grid_clipped_area2km.gpkg")
grid_clipped_area2km$area <- st_area(grid_clipped_area2km)/1e6
grid_clipped_wkt <- grid_clipped_area2km |>
  mutate(geometry = st_as_text(geom))  # convert geometry to readable WKT

write_csv(grid_clipped_wkt, "grid_clipped_wkt.csv")
