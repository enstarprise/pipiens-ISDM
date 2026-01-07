library(sf)
library(ncdf4)
library(future.apply)
library(dplyr)
library(tidyr)
library(stringr)
library(abind)
library(lubridate)
library(terra)
library(ggplot2)
plan(multisession)
rm(list = ls())

# --- scotland boundary ---
scotland <- st_read("gridExtractions/scotlandBoundary/scotland_boundary.shp")
scotland <- st_transform(scotland,crs = 27700)

grid_clipped <- st_read("grids/grid_clipped_2point5km.gpkg")

#-----------# 1. create Spatial version of grid---------------
grid_clipped$ID <- 1:nrow(grid_clipped)
grid_vect <- vect(grid_clipped) # Convert sf grid to SpatVector

# Get CRS from any temperature raster
r0 <- rast(list.files("gridExtractions/covariates/daily", pattern = "tasmax.*\\.nc$", full.names = TRUE)[1])

# Reproject your grid to match raster CRS
grid_clipped <- st_transform(grid_clipped, crs(r0))

#-----------# 2. FUNCTIONS ---------------
### --- file loader --- ####
load_nc_temp_data_multi_fixed <- function(nc_files, var_names = c("tasmin","tasmax")) {
  lapply(nc_files, function(file) {
    nc <- nc_open(file)
    
    # Time variable with correct conversion for "hours since 1800-01-01"
    time <- ncvar_get(nc, "time")
    units <- ncatt_get(nc, "time", "units")$value
    
    # Convert hours since 1800-01-01 to dates
    origin <- as.POSIXct("1800-01-01 00:00:00", tz = "UTC")
    dates <- as.Date(origin + hours(time))
    
    # Only keep variables present in this file
    vars_to_read <- intersect(var_names, names(nc$var))
    if(length(vars_to_read) == 0) {
      nc_close(nc)
      stop("No requested variables in: ", file)
    }
    
    vars <- lapply(vars_to_read, function(v) ncvar_get(nc, v))
    names(vars) <- vars_to_read
    
    nc_close(nc)
    list(file = file, temperature = vars, dates = dates)
  })
}


#### ---- EXTRACTION OF 7 DAY MEAN TEMP (TASMIN & TASMAX) --- ####

extract_temp_simple_terra_fixed <- function(grid_sf, dates_to_extract, nc_dir) {
  
  # List temperature files
  tasmin_files <- list.files(nc_dir, pattern = "tasmin.*\\.nc$", full.names = TRUE)
  tasmax_files <- list.files(nc_dir, pattern = "tasmax.*\\.nc$", full.names = TRUE)
  
  if(length(tasmin_files) == 0 & length(tasmax_files) == 0) {
    stop("No temperature files found in: ", nc_dir)
  }
  
  message("Found ", length(tasmin_files), " tasmin files and ", length(tasmax_files), " tasmax files")
  
  results <- list()
  
  for(target_date in dates_to_extract) {
    start_date <- target_date - 6
    message("Processing 7-day window: ", format(as.Date(start_date, origin="1970-01-01")),
            " to ", format(as.Date(target_date, origin="1970-01-01")))
    
    tryCatch({
      # ---- tasmin stack ----
      tasmin_stack <- rast()
      for(f in tasmin_files) {
        file_dates <- as.Date(strsplit(sub(".*_(\\d{8}-\\d{8})\\.nc$", "\\1", f), "-")[[1]], "%Y%m%d")
        if(file_dates[1] <= target_date & file_dates[2] >= start_date) {
          r <- rast(f)
          time_vals <- time(r)
          if(!is.null(time_vals)) {
            date_mask <- as.Date(time_vals) >= start_date & as.Date(time_vals) <= target_date
            if(any(date_mask)) {
              tasmin_stack <- c(tasmin_stack, r[[which(date_mask)]])
            }
          }
        }
      }
      
      # ---- tasmax stack ----
      tasmax_stack <- rast()
      for(f in tasmax_files) {
        file_dates <- as.Date(strsplit(sub(".*_(\\d{8}-\\d{8})\\.nc$", "\\1", f), "-")[[1]], "%Y%m%d")
        if(file_dates[1] <= target_date & file_dates[2] >= start_date) {
          r <- rast(f)
          time_vals <- time(r)
          if(!is.null(time_vals)) {
            date_mask <- as.Date(time_vals) >= start_date & as.Date(time_vals) <= target_date
            if(any(date_mask)) {
              tasmax_stack <- c(tasmax_stack, r[[which(date_mask)]])
            }
          }
        }
      }
      
      # ---- Calculate combined 7-day mean ----
      combined_temp <- NULL
      if(nlyr(tasmin_stack) > 0 & nlyr(tasmax_stack) > 0) {
        mean_tasmin <- mean(tasmin_stack, na.rm = TRUE)
        mean_tasmax <- mean(tasmax_stack, na.rm = TRUE)
        combined_temp <- mean(c(mean_tasmin, mean_tasmax), na.rm = TRUE)
      } else if(nlyr(tasmin_stack) > 0) {
        combined_temp <- mean(tasmin_stack, na.rm = TRUE)
      } else if(nlyr(tasmax_stack) > 0) {
        combined_temp <- mean(tasmax_stack, na.rm = TRUE)
      } else {
        warning("No raster data for 7-day window: ", format(start_date), " to ", format(target_date))
        next
      }
      
      # ---- Extract values ----
      extracted <- extract(combined_temp, vect(grid_sf), fun = mean, na.rm = TRUE)
      
      results[[as.character(target_date)]] <- data.frame(
        grid_id = grid_sf$grid_id,
        start_date = as.Date(start_date, origin = "1970-01-01"),
        end_date   = as.Date(target_date, origin = "1970-01-01"),
        mean_temp_7d_celsius = round(extracted[,2], 1)
      )
      
      
    }, error = function(e) {
      warning("Error processing ", format(target_date), ": ", e$message)
    })
  }
  
  if(length(results) == 0) return(data.frame())
  do.call(rbind, results)
}

# -----------------------------
# Compute 4 lag windows per observed date
# -----------------------------
get_lagged_temps_terra_fixed <- function(grid_sf, obs_dates, nc_dir, max_lag_weeks = 4) {
  
  # 1. Build lag windows
  lag_df <- expand.grid(
    obs_date = as.Date(obs_dates),
    lag_week = 1:max_lag_weeks
  ) %>%
    mutate(
      end_date   = obs_date - 7 * (lag_week - 1),  # End of 7-day lag window
      start_date = end_date - 6                     # Start of 7-day window
    )
  
  unique_end_dates <- sort(unique(lag_df$end_date))
  
  # 2. Extract 7-day mean temps for all unique end_dates
  temp_df <- extract_temp_simple_terra_fixed(
    grid_sf = grid_sf,
    dates_to_extract = unique_end_dates,
    nc_dir = nc_dir
  )
  
  if(nrow(temp_df) == 0) {
    stop("No temperature data extracted. Check if lag windows fall within NetCDF date ranges.")
  }
  
  # 3. Join lag info with extracted temperatures
  lagged_df <- lag_df %>%
    left_join(temp_df, by = c("start_date", "end_date")) %>%
    select(grid_id, obs_date, lag_week, start_date, mean_temp_7d_celsius)
  
  # 4. Ensure uniqueness for pivot
  lagged_df_unique <- lagged_df %>%
    group_by(grid_id, obs_date, lag_week) %>%
    summarize(
      mean_temp_7d_celsius = mean(mean_temp_7d_celsius, na.rm = TRUE),
      start_date = min(start_date),
      .groups = "drop"
    )
  
  # 5. Pivot lag temperatures
  lagged_df_wide <- lagged_df_unique %>%
    pivot_wider(
      id_cols = c(grid_id, obs_date),
      names_from = lag_week,
      values_from = mean_temp_7d_celsius,
      names_glue = "lag{lag_week}"
    )
  
  # 6. Pivot lag start dates
  start_dates_wide <- lagged_df_unique %>%
    pivot_wider(
      id_cols = c(grid_id, obs_date),
      names_from = lag_week,
      values_from = start_date,
      names_glue = "start_lag{lag_week}"
    )
  
  # 7. Combine temperatures and start dates
  final_wide <- lagged_df_wide %>%
    left_join(start_dates_wide, by = c("grid_id", "obs_date"))
  
  return(final_wide)
}



# -----------------------------
# usage:
# -----------------------------
dates_to_extract <- union(adult_data$Collection_date, citizenScience_data$Date_found)
test_dates <- as.Date(c("2023-06-15", "2023-07-15"))

results <- get_lagged_temps_terra_fixed(
  grid_sf = grid_clipped,
  obs_dates = dates_to_extract,
  nc_dir = "gridExtractions/covariates/daily",
  max_lag_weeks = 4
)

write.csv(results, "csvs/results_lag.csv")


# Reshape results to long format
results_long <- results %>%
  pivot_longer(
    cols = starts_with("lag"),    # assumes columns are lag1, lag2, ...
    names_to = "lag",
    values_to = "temperature"
  )

# Convert lag names to more descriptive labels
results_long <- results_long %>%
  mutate(lag_label = case_when(
    lag == "lag1" ~ "1-week lag",
    lag == "lag2" ~ "2-week lag",
    lag == "lag3" ~ "3-week lag",
    lag == "lag4" ~ "4-week lag",
    TRUE ~ lag
  ))

# Compute mean temperature per lag
lag_means <- results_long %>%
  group_by(lag_label) %>%
  summarize(mean_temp = mean(temperature, na.rm = TRUE), .groups = "drop")

# Plot with facets for each lag
lag_hist_plot <- ggplot(results_long, aes(x = temperature)) +
  geom_histogram(bins = 20, fill = "skyblue", color = "white", alpha = 0.7) +
  geom_vline(data = lag_means, aes(xintercept = mean_temp), 
             linetype = "dashed", color = "red", size = 1) +
  facet_wrap(~lag_label, ncol = 2, scales = "free") +  # one facet per lag
  labs(
    title = "Temperature Distribution by Lag",
    x = "Temperature (°C)",
    y = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold")
  )
# Save the plot as PNG
ggsave(
  filename = "lag_histograms.png",  # file name
  plot = lag_hist_plot,             # the ggplot object
  width = 10,                       # width in inches
  height = 8,                       # height in inches
  dpi = 300                         # resolution
)



# -----------------------------
# test usage:
# -----------------------------
test_dates <- as.Date(c("2023-06-15", "2023-07-15"))

# First, let's verify our test dates are within file ranges
check_date_availability <- function(test_dates, nc_dir) {
  nc_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
  
  for(td in test_dates) {
    start_date <- td - 6
    cat("Checking date range:", as.character(start_date), "to", as.character(td), "\n")
    
    matching_files <- nc_files[sapply(nc_files, function(f) {
      date_match <- regmatches(f, regexpr("\\d{8}-\\d{8}", f))
      if(length(date_match) == 0) return(FALSE)
      dates <- as.Date(strsplit(date_match, "-")[[1]], "%Y%m%d")
      dates[1] <= td & dates[2] >= start_date
    })]
    
    cat("  Matching files:", length(matching_files), "\n")
    if(length(matching_files) > 0) {
      print(basename(matching_files))
    }
    cat("\n")
  }
}

check_date_availability(test_dates, "gridExtractions/covariates/daily")

# Test the simple version
test_dates <- as.Date(c("2023-06-15", "2023-07-15"))

simple_result <- extract_temp_simple_terra(
  grid_sf = grid_clipped,
  dates_to_extract = test_dates,
  nc_dir = "gridExtractions/covariates/daily"
)
head(simple_result)

lagged_result <- get_lagged_temps_terra_fixed(
  grid_sf = grid_clipped,
  obs_dates = test_dates,
  nc_dir = "gridExtractions/covariates/daily",
  max_lag_weeks = 4
)
head(lagged_result)



#### run the extraction ####

dates_to_extract <- union(adult_data$Collection_date, citizenScience_data$Date_found)

results <- get_lagged_temps_terra_fixed(
  grid_sf = grid_clipped,
  obs_dates = dates_to_extract,
  nc_dir = "gridExtractions/covariates/daily",
  max_lag_weeks = 4
)






