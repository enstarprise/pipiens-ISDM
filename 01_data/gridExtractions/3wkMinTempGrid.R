
################ 
#> MEAN TEMP CALCULATED BELOW
#> must first do 1wkMaxTempGrid.R first
#> the grid cells here are 2.5km^2 instead
#> of 1km, so the previous functions must
#> be adapted...
#################

library(sf)
library(ncdf4)
library(terra)   
library(future.apply)
library(progressr)
library(abind)
library(lubridate)
library(stringr)
library(readxl)
plan(multisession)  # For parallel processing

# --- scotland boundary ---
scotland <- st_read("01_data/gridExtractions/scotlandBoundary/scotland_boundary.shp")
scotland <- st_transform(scotland,crs = 27700)

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

grid_clipped <- st_read("01_data/grids/grid_clipped_1km.gpkg")
# grid_clipped <- st_read("01_data/grids/grid_clipped_2point5km.gpkg")

#-----------# 1. create Spatial version of grid---------------
grid_clipped$ID <- 1:nrow(grid_clipped)
grid_vect <- vect(grid_clipped) # Convert sf grid to SpatVector


# grid centers from sf
grid_coords <- st_centroid(grid_clipped) |> st_coordinates()
grid_df <- data.frame(
  grid_id = grid_clipped$grid_id,  
  x_bng = grid_coords[,1],
  y_bng = grid_coords[,2]
)


# ---- 0. Functions -----

# 1. Modified load_nc_data for temperature
load_nc_temp_data <- function(nc_dir, var_name = "tasmin") {
  lapply(nc_dir, function(file) {
    nc <- nc_open(file)
    temperature <- ncvar_get(nc, var_name)
    time <- ncvar_get(nc, "time")
    units <- ncatt_get(nc, "time", "units")$value
    origin <- as.POSIXct(sub("hours since ", "", units), tz = "UTC")
    dates <- as.Date(origin + hours(time))
    nc_close(nc)
    list(file = file, temperature = temperature, dates = dates)
  })
}

# 2. find_nearest_grid_bng (same as before)
find_nearest_grid_bng <- function(x, y, x_coords, y_coords) {
  x_idx <- which.min(abs(x_coords - x))
  y_idx <- which.min(abs(y_coords - y))
  return(c(x_idx, y_idx))
}

# 3. Modified to calculate 21-day mean min temp
get_mean_min_temp <- function(target_date, x_idx, y_idx, nc_data_list) {
  start_date <- target_date - 20  # 20-day window (inclusive)
  temp_values <- numeric(0)
  
  message("Target range: ", start_date, " to ", target_date)
  
  for (i in seq_along(nc_data_list)) {
    dates_in_file <- nc_data_list[[i]]$dates
    in_range <- dates_in_file >= start_date & dates_in_file <= target_date
    
    if (any(in_range)) {
      temp_vals <- nc_data_list[[i]]$temperature[x_idx, y_idx, in_range]
      temp_values <- c(temp_values, temp_vals)
    }
  }
  
  if (length(temp_values) == 0) {
    message("No temperature data found for this period")
    return(NA)
  }
  
  mean_temp <- mean(temp_values, na.rm = TRUE)
  message("\n7-DAY MEAN MIN TEMP: ", mean_temp)
  return(mean_temp)
}

# 4. Modified function to get correct NetCDF files for 20-day window
filter_nc_files_for_temp_range <- function(nc_dir, start_date, end_date) {
  nc_files <- list.files(nc_dir, pattern = "tasmin_hadukgrid_uk_1km_day_.*\\.nc$", full.names = TRUE)
  
  # Extract date ranges from filenames
  file_dates <- str_extract(nc_files, "\\d{8}-\\d{8}")
  
  matching_files <- sapply(seq_along(nc_files), function(i) {
    if (is.na(file_dates[i])) return(FALSE)
    
    file_date_range <- strsplit(file_dates[i], "-")[[1]]
    file_start <- as.Date(file_date_range[1], format = "%Y%m%d")
    file_end <- as.Date(file_date_range[2], format = "%Y%m%d")
    
    (file_start <= end_date) && (file_end >= start_date)
  })
  
  return(nc_files[matching_files])
}

# 5. Modified main extraction function for temperature
# extract_grid_temp_ts_optimized <- function(grid_sf, dates_to_extract, nc_dir) {
  # --- 1. Create a raster template from the first NetCDF file
#   nc_file <- list.files(nc_dir, pattern = "tasmin.*\\.nc$", full.names = TRUE)[1]
#   r_template <- rast(nc_file)  # <-- no varname
# 
#   # --- 2. Precompute cell indices for each 2.5 km polygon
#   message("Precomputing 1 km cell indices per polygon...")
#   cell_list <- lapply(1:nrow(grid_sf), function(i) {
#     cells(r_template, vect(grid_sf[i, ]))[[1]]
#   })
# 
# 
#   # Convert cell indices to x/y NetCDF indices
#   xy_idx_list <- lapply(cell_list, function(cells_idx) {
#     if (length(cells_idx) == 0) return(matrix(NA, ncol = 2))
#     rc <- terra::rowColFromCell(r_template, cells_idx)
#     cbind(rc[, 2], rc[, 1])
#   })
# 
#   # --- 3. Extract coordinate arrays for NetCDF indexing
#   nc_ref <- nc_open(nc_file)
#   x_coords <- ncvar_get(nc_ref, "projection_x_coordinate")
#   y_coords <- ncvar_get(nc_ref, "projection_y_coordinate")
#   nc_close(nc_ref)
# 
#   # --- 4. Process each target date in parallel
#   message("Extracting 21-day mean min temperatures...")
# 
#   results <- future_lapply(dates_to_extract, function(target_date) {
#     start_date <- target_date - 20
# 
#     # Find all relevant NetCDF files
#     nc_files <- filter_nc_files_for_temp_range(nc_dir, start_date, target_date)
#     if (length(nc_files) == 0) return(NULL)
# 
#     # Load and stack temperature data for the 21-day window
#     temp_data <- load_nc_temp_data(nc_files)
#     all_dates <- unlist(lapply(temp_data, `[[`, "dates"))
#     temp_stack <- abind::abind(lapply(temp_data, `[[`, "temperature"), along = 3)
#     date_mask <- (all_dates >= start_date) & (all_dates <= target_date)
#     if (!any(date_mask)) return(NULL)
# 
#     # --- 5. Extract for each polygon
#     temp_values <- vapply(seq_along(xy_idx_list), function(i) {
#       xy_idx <- xy_idx_list[[i]]
#       if (any(is.na(xy_idx))) return(NA_real_)
# 
#       cell_temps <- unlist(
#         mapply(function(xi, yi) {
#           temp_stack[xi, yi, date_mask]
#         }, xy_idx[,1], xy_idx[,2], SIMPLIFY = FALSE)
#       )
# 
#       mean(cell_temps, na.rm = TRUE)
#     }, numeric(1))
# 
# 
#     data.frame(
#       grid_id = grid_sf$grid_id,
#       date = target_date,
#       mean_min_temp_7d_celsius = round(temp_values, 1),
#       start_date = start_date,
#       end_date = target_date,
#       stringsAsFactors = FALSE
#     )
#   }, future.seed = TRUE)
# 
#   # --- 6. Combine and return
#   result_df <- do.call(rbind, results[!sapply(results, is.null)])
#   return(result_df)
# }

extract_grid_temp_ts_optimized <- function(grid_df, dates_to_extract, nc_dir, n_cores = NULL) {
  library(future.apply)
  library(parallel)
  library(ncdf4)
  library(abind)
  
  # Set up parallel backend
  if (is.null(n_cores)) {
    n_cores <- max(1, detectCores() - 1)
  }
  
  # Configure parallel strategy with explicit globals export
  plan(multisession, workers = n_cores)
  
  # Load first file to get spatial reference
  nc_ref <- nc_open(list.files(nc_dir, pattern = "tasmin.*\\.nc$", full.names = TRUE)[1])
  x_coords <- ncvar_get(nc_ref, "projection_x_coordinate")
  y_coords <- ncvar_get(nc_ref, "projection_y_coordinate")
  nc_close(nc_ref)
  
  # Pre-compute all grid indices (vectorized)
  grid_idx <- t(apply(grid_df[, c("x_bng", "y_bng")], 1, function(coord) {
    find_nearest_grid_bng(coord[1], coord[2], x_coords, y_coords)
  }))
  colnames(grid_idx) <- c("grid_x", "grid_y")
  grid_df <- cbind(grid_df, grid_idx)
  
  # Check if parallelization is worth it
  cat(sprintf("Processing %d dates with %d cores\n", length(dates_to_extract), n_cores))
  
  # Process dates in parallel with explicit globals
  results <- future_lapply(dates_to_extract, function(target_date) {
    # Explicitly load libraries in each worker
    library(ncdf4)
    library(abind)
    
    # Load required NetCDF files
    nc_files <- filter_nc_files_for_temp_range(nc_dir, target_date - 20, target_date)
    if (length(nc_files) == 0) return(NULL)
    
    # Load and stack temperature data
    temp_data <- load_nc_temp_data(nc_files)
    all_dates <- unlist(lapply(temp_data, `[[`, "dates"))
    temp_stack <- abind::abind(lapply(temp_data, `[[`, "temperature"), along = 3)
    
    # Create date mask for the 21-day window
    date_mask <- (all_dates >= (target_date - 20)) & (all_dates <= target_date)
    if (!any(date_mask)) return(NULL)
    
    # Vectorized extraction for all grids
    valid_mask <- !is.na(grid_df$grid_x) & !is.na(grid_df$grid_y)
    temp_values <- rep(NA_real_, nrow(grid_df))
    
    if (any(valid_mask)) {
      indices <- cbind(grid_df$grid_x[valid_mask], 
                       grid_df$grid_y[valid_mask])
      temp_values[valid_mask] <- apply(indices, 1, function(idx) {
        mean(temp_stack[idx[1], idx[2], date_mask], na.rm = TRUE)
      })
    }
    
    data.frame(
      grid_id = grid_df$grid_id,
      date = target_date,
      x_bng = grid_df$x_bng,
      y_bng = grid_df$y_bng,
      mean_min_temp_7d_celsius = round(temp_values, 1),
      start_date = target_date - 20,
      end_date = target_date,
      stringsAsFactors = FALSE
    )
  }, future.seed = TRUE,
  future.globals = structure(TRUE, add = c("grid_df", "nc_dir", 
                                           "filter_nc_files_for_temp_range",
                                           "load_nc_temp_data")))
  
  # Reset to sequential processing
  plan(sequential)
  
  # Combine and return results
  do.call(rbind, results[!sapply(results, is.null)])
}


# 6. get dates for extraction
get_unique_extraction_dates <- function(..., date_columns) {
  dfs <- list(...)
  all_dates <- unlist(mapply(function(df, col) df[[col]], dfs, date_columns, SIMPLIFY = FALSE))
  unique_dates <- sort(unique(as.Date(all_dates)))
  return(unique_dates)
}


# ---- 1. Setup ----
# Verify coordinate ranges with temperature files
nc_first <- nc_open(list.files("01_data/gridExtractions/covariates/daily/", 
                               pattern = "tasmin_hadukgrid_uk_1km_day_.*\\.nc$", 
                               full.names = TRUE)[1])
x_coords <- ncvar_get(nc_first, "projection_x_coordinate")
y_coords <- ncvar_get(nc_first, "projection_y_coordinate")
nc_close(nc_first)

print(paste("Grid X range:", min(grid_df$x_bng), "-", max(grid_df$x_bng)))
print(paste("NetCDF X range:", min(x_coords), "-", max(x_coords)))
print(paste("Grid Y range:", min(grid_df$y_bng), "-", max(grid_df$y_bng)))
print(paste("NetCDF Y range:", min(y_coords), "-", max(y_coords)))


# ---- 2. TEST FUNCTIONS -----
# Test 1: Basic file filtering
test_date <- as.Date("2023-03-015")
test_files <- filter_nc_files_for_temp_range(
  "01_data/gridExtractions/covariates/daily/", 
  test_date - 20, 
  test_date
)
print(test_files)

# Test file loading independently
test_files <- c(
  "01_data/gridExtractions/covariates/daily/tasmin_hadukgrid_uk_1km_day_20230601-20230630.nc",
  "01_data/gridExtractions/covariates/daily/tasmin_hadukgrid_uk_1km_day_20230701-20230731.nc"
)

test_data <- load_nc_temp_data(test_files)
# Check what dates are actually loaded
print("Dates in second file:");print(test_data[[1]]$dates)
print("Dates in second file:");print(test_data[[2]]$dates)

test_date <- as.Date("2023-07-01")
start_date <- test_date - 20
message("Checking date range: ", start_date, " to ", test_date)

# Check which dates fall in this range
for (i in seq_along(test_data)) {
  in_range <- test_data[[i]]$dates >= start_date & test_data[[i]]$dates <= test_date
  message("File ", basename(test_data[[i]]$file), " has ", sum(in_range), " matching date(s)")
  print(test_data[[i]]$dates[in_range])
}


# Test 3: Parallel extraction
test_dates <- as.Date(c("2023-06-15", "2023-07-15"))

# test_results <- extract_grid_temp_ts_optimized(
#   grid_clipped[1:5, ],  
#   test_dates, 
#   "gridExtractions/covariates/daily/"
# )

# print(test_results)

# Manually check values at known locations
test_idx <- find_nearest_grid_bng(213029.2, 531101.4, x_coords, y_coords)
print(test_idx)

# Extract values directly
nc <- nc_open(test_files[1])
raw_values <- ncvar_get(nc, "tasmin")[test_idx[1], test_idx[2], ]
print(raw_values)  # Should be 250-320 Kelvin (reasonable temperatures)
nc_close(nc)

# ---- 3. Full Extraction ----

survey_df_new <- read_csv("01_data/csvs/survey_df.csv")
survey_df_new$Setup_date <- as.Date(survey_df_new$Setup_date, 
                                    format = "%d/%m/%Y")


cs_df_new <- read_csv("01_data/csvs/cs_df.csv")
cs_df_new <- cs_df_new %>% filter(Verified_mosquito == "Yes" & 
                                    Species == "pipiens" & 
                                    Stage == "Adult")

cs_df_new$Date_found <- as.Date(cs_df_new$Date_found, 
                                format = "%d/%m/%Y")

dates_to_extract <- get_unique_extraction_dates(
  survey_df_new, cs_df_new, 
  date_columns = c("Setup_date", "Date_found")
)
dates_to_extract[1]

# dates_to_extract <- get_unique_extraction_dates(
#   surveillance_data, citizenScience_data,
#   date_columns = c("Collection_date", "Date_found")
# )

# Run the extraction
# temp_grid_ts <- extract_grid_temp_ts_optimized(
#   grid_df,
#   dates_to_extract,
#   "gridExtractions/covariates/daily/"  # Directory with temperature NetCDF files
# )
# temp_grid_ts <- extract_grid_temp_ts_optimized(grid_clipped,  # <- Use the sf object
#                                                dates_to_extract,
#                                                nc_dir = "gridExtractions/covariates/daily/")
# 
# result <- extract_grid_temp_ts_optimized(
#   grid_df = grid_df,
#   dates_to_extract = dates_to_extract,
#   nc_dir = "gridExtractions/covariates/daily/",
#   n_cores = 10 
# )


### #########jsut the new dates!!
# Filter for June–September
# citizenScience_data_newDates <- citizenScience_data %>%
#   filter(month(Date_found) %in% 6:9)
# 
# surveillance_data_newDates <- surveillance_data %>%
#   filter(month(Collection_date) %in% 2:4)
# 
# dates_to_extract <- get_unique_extraction_dates(
#   surveillance_data_newDates, citizenScience_data_newDates,
#   date_columns = c("Collection_date",  "Date_found")
# )

nc_dir = "01_data/gridExtractions/covariates/daily/" 

library(pbapply)
library(parallel)
library(sf)
library(stringr)      
library(lubridate)    
library(ncdf4)
library(abind)

# Extract coordinates from sf object
coords <- st_coordinates(st_centroid(grid_clipped))
grid_df <- data.frame(
  grid_id = grid_clipped$grid_id,
  x_bng = coords[, 1],
  y_bng = coords[, 2]
)

# Set up parallel cluster
cl <- makeCluster(detectCores() - 4)

# Export necessary objects and functions to cluster
clusterExport(cl, c("grid_df", "nc_dir", 
                    "extract_grid_temp_ts_optimized",
                    "filter_nc_files_for_temp_range",
                    "load_nc_temp_data",
                    "find_nearest_grid_bng"))

# Load required libraries on each worker (ADD stringr)
clusterEvalQ(cl, {
  library(ncdf4)
  library(abind)
  library(stringr) 
  library(lubridate)
})

# Run in parallel with progress bar
result <- pblapply(dates_to_extract, function(date) {
  extract_grid_temp_ts_optimized(
    grid_df = grid_df,
    dates_to_extract = date,
    nc_dir = nc_dir
  )
}, cl = cl)

# Stop cluster
stopCluster(cl)

# Combine all results
final_result <- do.call(rbind, result)

# View results
head(final_result)

# Check results
cat("Total rows:", nrow(final_result), "\n")
cat("Expected:", length(dates_to_extract) * nrow(grid_df), "\n")

tasmin_grid <- final_result

write.csv(tasmin_grid, 
          file = "01_data/csvs/tasmin_grid_21day_2point5km.csv",
          row.names = FALSE)


###### MEAN TEMP #####
# tasmin <- read_csv("csvs/tasmin_grid_2point5km.csv")
# tasmax <- read_csv("csvs/tasmax_grid.csv")
# 
# tasmean <- tasmax[, c("grid_id", "date", "mean_max_temp_7d_celsius")] %>%
#   inner_join(tasmin[, c("grid_id", "date", "mean_min_temp_7d_celsius")], by = c("grid_id", "date")) %>%
#   mutate(mean_temp_7d_celsius = (mean_max_temp_7d_celsius + mean_min_temp_7d_celsius) / 2)
# 
# write.csv(tasmean, 
#           file = "csvs/tasmean_grid.csv",
#           row.names = FALSE)
