rm(list = ls())

library(readxl)
library(dplyr)
library(lubridate)
library(svMisc);library(raster)
library(ncdf4) # package for netcdf manipulation
library(terra)
library(sf)
library(zoo)
library(tidyverse)

library(future)
library(furrr)
library(parallel)

# --- scotland boundary ---
scotland <- st_read("gridExtractions/scotlandBoundary/scotland_boundary.shp")
scotland <- st_transform(scotland,crs = 27700)

# grid_1km <- st_make_grid(scotland, cellsize = 1000, square = TRUE)
# grid_1km <- st_sf(grid_id = 1:length(grid_1km), geometry = grid_1km)
# st_write(grid_1km, "grids/grid_1km.gpkg")

# plan(multisession, workers = 10)  # use 10 cores # plan(multisession, workers = availableCores() - 1)
# n_chunks <- 40  # 40/10 = 4 chunks of the dataset per core
# grid_split <- split(grid_1km, cut(1:nrow(grid_1km), n_chunks, labels = FALSE))

# parallel intersection
# grid_clipped_list <- future_map(
#   grid_split,
#   ~ st_intersection(.x, scotland),
#   .options = furrr_options(seed = TRUE)
# ) # 5 minute run time

# grid_clipped_1km <- do.call(rbind, grid_clipped_list) # combine results
# plan(sequential) # reset to sequential
# st_write(grid_clipped_1km, "grids/grid_clipped_1km.gpkg")

grid_1km <- st_read("grids/grid_1km.gpkg")
grid_clipped_1km <- st_read("grids/grid_clipped_1km.gpkg")


#--------------------------
# 1. create Spatial version of grid
#--------------------------
grid_clipped_1km$ID <- 1:nrow(grid_clipped_1km)
grid_vect <- vect(grid_clipped_1km) # Convert sf grid to SpatVector

#--------------------------
# 2. function to process daily nc files into weekly
#--------------------------

process_weekly_temp_fast <- function(nc_file, varname = "tasmax", grid_vect) {
  message("Processing: ", nc_file)
  
  # Load raster stack
  r <- rast(nc_file, subds = varname)
  
  # Get time values and units
  nc <- nc_open(nc_file)
  time_vals <- ncvar_get(nc, "time")
  time_units <- ncatt_get(nc, "time", "units")$value
  nc_close(nc)
  
  # Safely extract the origin date
  origin_str <- sub(".*since ", "", time_units)
  origin_date <- tryCatch(as.Date(origin_str), error = function(e) NA)
  
  # Compute actual dates
  if (grepl("days", time_units)) {
    dates <- origin_date + time_vals
  } else if (grepl("hours", time_units)) {
    dates <- origin_date + lubridate::dhours(time_vals)
  } else {
    stop("Unrecognized time units: ", time_units)
  }
  
  # Name raster layers with date strings
  names(r) <- as.character(dates)
  
  # Group by weeks
  week_index <- floor_date(dates, unit = "week", week_start = 1)
  unique_weeks <- unique(week_index)
  
  # Average each week’s layers
  weekly_means <- lapply(unique_weeks, function(w) {
    idx <- which(week_index == w)
    if (length(idx) > 1) {
      r_week <- mean(r[[idx]], na.rm = TRUE)
    } else {
      r_week <- r[[idx]]  # if only one day in that week
    }
    names(r_week) <- as.character(w)
    return(r_week)
  })
  
  # Stack and extract to grid
  r_stack <- rast(weekly_means)
  grid_vect_proj <- project(grid_vect, crs(r_stack))
  extracted <- terra::extract(r_stack, grid_vect_proj, fun = mean, na.rm = TRUE, ID = FALSE)
  
  # Check and return in tidy format
  if (ncol(extracted) == 0) return(tibble())
  
  extracted$grid_id <- 1:nrow(extracted)
  
  extracted_long <- extracted %>%
    relocate(grid_id) %>%
    pivot_longer(
      cols = -grid_id,
      names_to = "week",
      values_to = "tasmax_avg"
    ) %>%
    mutate(week = as.Date(week))  # ensure week is date format
  
  return(extracted_long)
}

temp_files <- list.files("covariates/daily", 
                         pattern = "^tasmax_hadukgrid_uk_5km_day_.*\\.nc$", 
                         full.names = TRUE)

temp_weekly_all <- bind_rows(lapply(temp_files,
                                    process_weekly_temp_fast,
                                    grid_vect = grid_vect))


saveRDS(temp_weekly_all, "temp_weekly_all.RDS")


#--------------------------
# 2.1. save as raster brick
#--------------------------

# 1. Join weekly values to the 5x5 km grid geometry
temp_sf <- left_join(grid_clipped, temp_weekly_all, by = c("ID" = "grid_id"))

# 2. Create a blank template raster to match the grid resolution & extent
template_rast <- rast(ext(grid_clipped), resolution = 5000, crs = crs(grid_clipped))

# 3. Split by week, rasterize each week’s data
weekly_rasters <- temp_sf %>%
  group_by(week) %>%
  group_split() %>%
  map(~ {
    this_week <- unique(.x$week)
    r <- terra::rasterize(vect(.x), template_rast, field = "tasmax_avg")
    names(r) <- format(this_week, "%Y-%m-%d")
    r
  })

# 4. Combine into a single SpatRaster stack
temp_weekly_stack <- do.call(c, weekly_rasters)

# 5. Save to file
writeRaster(temp_weekly_stack, "temperature_weekly_stack.tif", overwrite = TRUE)
# Or NetCDF format (smaller and multidimensional)
writeCDF(temp_weekly_stack, "temperature_weekly_stack.nc", varname = "tasmax_avg", overwrite = TRUE)

names(temp_weekly_stack)
plot(temp_weekly_stack[[1]], main = names(temp_weekly_stack)[1])
hist(temp_stack[[1]], main = "Temp Distribution - Week 1")


#--------------------------
# 3. function to process daily rainfall files into week---
#--------------------------

process_weekly_rain_fast <- function(nc_file, varname = "rainfall", grid_vect) {
  message("Processing: ", nc_file)
  
  # Load raster stack
  r <- rast(nc_file, subds = varname)
  
  # Get time values and units
  nc <- nc_open(nc_file)
  time_vals <- ncvar_get(nc, "time")
  time_units <- ncatt_get(nc, "time", "units")$value
  nc_close(nc)
  
  # Safely extract the origin date
  origin_str <- sub(".*since ", "", time_units)
  origin_date <- tryCatch(as.Date(origin_str), error = function(e) NA)
  
  # Compute actual dates
  if (grepl("days", time_units)) {
    dates <- origin_date + time_vals
  } else if (grepl("hours", time_units)) {
    dates <- origin_date + lubridate::dhours(time_vals)
  } else {
    stop("Unrecognized time units: ", time_units)
  }
  
  # Name raster layers with date strings
  names(r) <- as.character(dates)
  
  # Group by weeks
  week_index <- floor_date(dates, unit = "week", week_start = 1)
  unique_weeks <- unique(week_index)
  
  # Average each week’s layers
  weekly_means <- lapply(unique_weeks, function(w) {
    idx <- which(week_index == w)
    if (length(idx) > 1) {
      r_week <- mean(r[[idx]], na.rm = TRUE)
    } else {
      r_week <- r[[idx]]  # if only one day in that week
    }
    names(r_week) <- as.character(w)
    return(r_week)
  })
  
  # Stack and extract to grid
  r_stack <- rast(weekly_means)
  grid_vect_proj <- project(grid_vect, crs(r_stack))
  extracted <- terra::extract(r_stack, grid_vect_proj, fun = mean, na.rm = TRUE, ID = FALSE)
  
  # Check and return in tidy format
  if (ncol(extracted) == 0) return(tibble())
  
  extracted$grid_id <- 1:nrow(extracted)
  
  extracted_long <- extracted %>%
    relocate(grid_id) %>%
    pivot_longer(
      cols = -grid_id,
      names_to = "week",
      values_to = "rainfall_avg"
    ) %>%
    mutate(week = as.Date(week))  # ensure week is date format
  
  return(extracted_long)
}


rain_files <- list.files("gridExtractions/covariates/daily", 
                         pattern = "^rainfall_hadukgrid_uk_5km_day_.*\\.nc$", 
                         full.names = TRUE)

rain_weekly_all <- bind_rows(lapply(rain_files, 
                                    process_weekly_rain_fast, 
                                    grid_vect = grid_vect))

saveRDS(rain_weekly_all, "rain_weekly_all.RDS")
#rain_weekly_all <- readRDS("rain_weekly_all.RDS")

#--------------------------
# 3.1. save as raster brick
#--------------------------

# 1. Join weekly values to the 5x5 km grid geometry
rain_sf <- left_join(grid_clipped, rain_weekly_all, by = c("ID" = "grid_id"))

# 2. Create a blank template raster to match the grid resolution & extent
#template_rast <- rast(ext(grid_clipped), resolution = 5000, crs = crs(grid_clipped))

# 3. Split by week, rasterize each week’s data
rain_weekly_rasters <- rain_sf %>%
  group_by(week) %>%
  group_split() %>%
  map(~ {
    this_week <- unique(.x$week)
    r <- terra::rasterize(vect(.x), template_rast, field = "rainfall_avg")
    names(r) <- format(this_week, "%Y-%m-%d")
    r
  })

# 4. Combine into a single SpatRaster stack
rain_weekly_stack <- do.call(c, rain_weekly_rasters)

# 5. Save to file
writeRaster(rain_weekly_stack, "rainfall_weekly_stack.tif", overwrite = TRUE)
# Or NetCDF format (smaller and multidimensional)
writeCDF(rain_weekly_stack, "rainfall_weekly_stack.nc", varname = "rainfall_avg", overwrite = TRUE)

#---------4wk cumulative precip-----------------
# 3.2. function to process daily rainfall files 
# into 4 week cumulative avg ---
#-----------------------------------------------

extract_4wk_cumulative_precip <- function(nc_files, target_dates, grid_vect, varname = "rainfall") {
  all_dates <- as.Date(target_dates)
  date_ranges <- lapply(all_dates, function(d) seq(d - 28, d - 1, by = "1 day"))
  unique_days <- sort(unique(do.call(c, date_ranges)))
  
  # 1. Extract dates from nc files
  load_nc_dates <- function(nc_file) {
    nc <- nc_open(nc_file)
    time_vals <- ncvar_get(nc, "time")
    time_units <- ncatt_get(nc, "time", "units")$value
    origin_str <- sub(".*since ", "", time_units)
    origin_date <- as.Date(origin_str)
    dates <- origin_date + time_vals
    nc_close(nc)
    tibble(file = nc_file, date = dates)
  }
  
  # 2. Build date-to-file index
  file_date_index <- bind_rows(lapply(nc_files, load_nc_dates))
  
  # 3. Filter for relevant files only
  needed_files <- file_date_index %>% filter(date %in% unique_days)
  
  # 4. Load rasters and extract for all dates needed
  all_extracts <- purrr::map_dfr(split(needed_files, needed_files$file), function(file_group) {
    file <- unique(file_group$file)
    dates <- file_group$date
    
    r <- rast(file, subds = varname)
    names(r) <- as.character(dates)
    
    r <- r[[which(names(r) %in% as.character(dates))]]
    r <- setNames(r, as.character(dates))
    
    # Project grid to raster CRS
    grid_vect_proj <- terra::project(grid_vect, crs(r))
    
    vals <- terra::extract(r, grid_vect_proj, ID = FALSE)
    vals$grid_id <- 1:nrow(vals)
    
    vals_long <- vals %>%
      pivot_longer(-grid_id, names_to = "date", values_to = "rain") %>%
      mutate(date = as.Date(date))
    
    return(vals_long)
  })
  
  # 5. Now for each target date, compute the 28-day cumulative per grid_id
  result <- purrr::map2_dfr(all_dates, 1:length(all_dates), function(target_date, idx) {
    window_dates <- seq(target_date - 28, target_date - 1, by = "1 day")
    df <- all_extracts %>%
      filter(date %in% window_dates, grid_id == idx) %>%
      summarise(rainfall_4wk = sum(rain, na.rm = TRUE)) %>%
      mutate(date_target = target_date, grid_id = idx)
    return(df)
  })
  
  return(result)
}

rain_files <- list.files("gridExtractions/covariates/daily", 
                         pattern = "^rainfall_hadukgrid_uk_5km_day_.*\\.nc$", 
                         full.names = TRUE)

# Assuming cs_df has a geometry column matching grid_vect order
cs_rainfall <- extract_4wk_cumulative_precip(
  nc_files = rain_files,
  target_dates = cs_sf$Date_found,
  grid_vect = grid_vect
)

# Same for adult_df
adult_rainfall <- extract_4wk_cumulative_precip(
  nc_files = rain_files,
  target_dates = adult_df$collection_date,
  grid_vect = grid_vect
)


#--------------------------
# 4. merge your weekly temperature and rainfall data 
#--------------------------
# duplicates for some combinations of grid_id and week, 
# there are multiple matching rows in both rain_weekly_all and temp_weekly_all
rain_weekly_all <- rain_weekly_all %>%
  group_by(grid_id, week) %>%
  summarise(rainfall_avg = mean(rainfall_avg, na.rm = TRUE), .groups = "drop")

temp_weekly_all <- temp_weekly_all %>%
  group_by(grid_id, week) %>%
  summarise(tasmax_avg = mean(tasmax_avg, na.rm = TRUE), .groups = "drop")

# 1. Merge rainfall and temperature by grid and week
climate_weekly <- full_join(rain_weekly_all, temp_weekly_all,
                            by = c("grid_id", "week"))

# 2. Merge with spatial grid
grid_all_weekly <- grid_clipped %>%
  left_join(climate_weekly, by = "grid_id")

# 3. Optional: Save as long-format CSV without geometry
climate_weekly_long <- grid_all_weekly %>%
  st_drop_geometry() %>%
  dplyr::select(grid_id, week, rainfall_avg, tasmax_avg)

write.csv(climate_weekly_long, "grid_weekly_climate_long.csv", row.names = FALSE)

# 4. Optional: Save with geometry as GeoPackage
st_write(grid_all_weekly, "grid_weekly_climate.gpkg", delete_dsn = TRUE)

