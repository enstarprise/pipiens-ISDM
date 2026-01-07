# ---- 0. data prep ----
library(ncdf4)
library(lubridate)
library(dplyr)
library(future.apply)
library(zoo)
library(sf)
library(ggplot2)
library(readxl)
library(stringr)

scotland <- st_read("gridExtractions/scotlandBoundary/scotland_boundary.shp")
scotland <- st_transform(scotland,crs = 27700)
# grid_1km <- st_read("grids/grid_1km.gpkg")
# grid_clipped <- st_read("grids/grid_clipped.gpkg")
grid_clipped <- st_read("grids/grid_clipped_1km.gpkg")


# grid centers from sf
grid_coords <- st_centroid(grid_clipped) |> st_coordinates()
grid_df <- data.frame(
  grid_id = grid_clipped$grid_id,  
  x_bng = grid_coords[,1],
  y_bng = grid_coords[,2]
)

nc_files <- list.files("gridExtractions/covariates/daily/",
                       pattern = "rainfall.*\\.nc$", full.names = TRUE)

nc_first <- nc_open(nc_files[1])
x_coords <- ncvar_get(nc_first, "projection_x_coordinate")
y_coords <- ncvar_get(nc_first, "projection_y_coordinate")
nc_close(nc_first)

# are points are within the ncdf grid bounds
print(paste("Your x range:", min(grid_df$x_bng), "-", max(grid_df$x_bng)))
print(paste("NetCDF x range:", min(x_coords), "-", max(x_coords)))

print(paste("Your y range:", min(grid_df$y_bng), "-", max(grid_df$y_bng)))
print(paste("NetCDF y range:", min(y_coords), "-", max(y_coords)))

# Keep only grids within NetCDF bounds
valid_grids <- grid_df %>%
  filter(
    x_bng >= min(x_coords),
    x_bng <= max(x_coords),
    y_bng >= min(y_coords),
    y_bng <= max(y_coords)
  )

# Visualize valid vs invalid points
ggplot() +
  geom_rect(aes(xmin=min(x_coords), xmax=max(x_coords),
                ymin=min(y_coords), ymax=max(y_coords)),
            fill="lightblue", alpha=0.3) +
  geom_point(data=grid_df, aes(x_bng, y_bng), color="red", size=1) +
  geom_point(data=valid_grids, aes(x_bng, y_bng), color="green", size=1) +
  labs(title="Grid Coverage (Green=Valid, Red=Outside NetCDF)") +
  theme_bw()

# ---- 0. the Functions ----
# 1. load_nc_data 
load_nc_data <- function(nc_dir) {
  lapply(nc_dir, function(file) {
    nc <- nc_open(file)
    rainfall <- ncvar_get(nc, "rainfall")
    time <- ncvar_get(nc, "time")
    units <- ncatt_get(nc, "time", "units")$value
    origin <- as.POSIXct(sub("hours since ", "", units), tz = "UTC")
    dates <- as.Date(origin + hours(time))
    nc_close(nc)
    list(file = file, rainfall = rainfall, dates = dates)
  })
}

# 2. find_nearest_grid for BNG coordinates
find_nearest_grid_bng <- function(x, y, x_coords, y_coords) {
  x_idx <- which.min(abs(x_coords - x))
  y_idx <- which.min(abs(y_coords - y))
  return(c(x_idx, y_idx))
}

# 3. 
get_28d_cumulative_preloaded <- function(target_date, x_idx, y_idx, nc_data_list) {
  start_date <- target_date - 27
  total <- 0
  
  # message("\n==== DEBUGGING 28-DAY SUM ====")
  message("Target range: ", start_date, " to ", target_date)
  
  for (i in seq_along(nc_data_list)) {
    dates_in_file <- nc_data_list[[i]]$dates
    in_range <- dates_in_file >= start_date & dates_in_file <= target_date
    
    if (any(in_range)) {
      rainfall_vals <- nc_data_list[[i]]$rainfall[x_idx, y_idx, in_range]
      dates_vals <- dates_in_file[in_range]
      
      # message("\nFile: ", basename(nc_data_list[[i]]$file))
      # print(data.frame(
      #   Date = dates_vals,
      #   Rainfall = rainfall_vals,
      #   Cumulative = cumsum(rainfall_vals)
      # ))
      
      total <- total + sum(rainfall_vals, na.rm = TRUE)
    }
  }
  
  message("\nTOTAL 28-DAY RAINFALL: ", total)
  return(total)
}

# 4. Function to get the correct NetCDF files for the 28-day window
filter_nc_files_for_date_range <- function(nc_dir, start_date, end_date) {
  # Get all rainfall files
  nc_files <- list.files(nc_dir, pattern = "rainfall.*\\.nc$", full.names = TRUE)
  
  # Extract date ranges from filenames
  file_dates <- str_extract(nc_files, "\\d{8}-\\d{8}")
  
  # Find files that overlap with our target date range
  matching_files <- sapply(seq_along(nc_files), function(i) {
    if (is.na(file_dates[i])) return(FALSE)
    
    # Parse file's start and end dates
    file_date_range <- strsplit(file_dates[i], "-")[[1]]
    file_start <- as.Date(file_date_range[1], format = "%Y%m%d")
    file_end <- as.Date(file_date_range[2], format = "%Y%m%d")
    
    # Check if file's date range overlaps with our target range
    (file_start <= end_date) && (file_end >= start_date)
  })
  
  return(nc_files[matching_files])
}

# 5. extract_grid_rainfall_ts_parallel
extract_grid_rainfall_ts_optimized <- function(grid_df, dates_to_extract, nc_dir) {
  # Load first file to get spatial reference
  nc_ref <- nc_open(list.files(nc_dir, pattern = "rainfall.*\\.nc$", full.names = TRUE)[1])
  x_coords <- ncvar_get(nc_ref, "projection_x_coordinate")
  y_coords <- ncvar_get(nc_ref, "projection_y_coordinate")
  nc_close(nc_ref)
  
  # Pre-compute all grid indices (vectorized)
  grid_idx <- t(apply(grid_df[, c("x_bng", "y_bng")], 1, function(coord) {
    find_nearest_grid_bng(coord[1], coord[2], x_coords, y_coords)
  }))
  colnames(grid_idx) <- c("grid_x", "grid_y")
  grid_df <- cbind(grid_df, grid_idx)
  
  # Initialize NA counter
  na_counter <- list(
    total_cells = 0,
    na_cells = 0,
    na_per_date = setNames(rep(0, length(dates_to_extract)), as.character(dates_to_extract))
  )
  
  # Process dates in parallel
  results <- future_lapply(dates_to_extract, function(target_date) {
    # Load required NetCDF files
    nc_files <- filter_nc_files_for_date_range(nc_dir, target_date - 27, target_date)
    if (length(nc_files) == 0) return(NULL)
    
    # Load and stack rainfall data
    nc_data <- load_nc_data(nc_files)
    all_dates <- unlist(lapply(nc_data, `[[`, "dates"))
    rainfall_stack <- abind::abind(lapply(nc_data, `[[`, "rainfall"), along = 3)
    
    # Create date mask for the 28-day window
    date_mask <- (all_dates >= (target_date - 27)) & (all_dates <= target_date)
    if (!any(date_mask)) return(NULL)
    
    # Vectorized extraction for all grids
    rainfall_values <- sapply(1:nrow(grid_df), function(i) {
      if (any(is.na(grid_df[i, c("grid_x", "grid_y")]))) {
        return(NA)
      }
      sum(rainfall_stack[grid_df$grid_x[i], grid_df$grid_y[i], date_mask], na.rm = TRUE)
    })
    
    # Update NA counter
    date_na_count <- sum(is.na(rainfall_values))
    na_counter$na_per_date[as.character(target_date)] <<- date_na_count
    na_counter$total_cells <<- na_counter$total_cells + nrow(grid_df)
    na_counter$na_cells <<- na_counter$na_cells + date_na_count
    
    data.frame(
      grid_id = grid_df$grid_id,
      date = target_date,
      x_bng = grid_df$x_bng,
      y_bng = grid_df$y_bng,
      rainfall_28d = rainfall_values,
      start_date = target_date - 27,
      end_date = target_date,
      stringsAsFactors = FALSE
    )
  }, future.seed = TRUE)
  
  # Combine and return results
  list(
    rainfall_data = do.call(rbind, results[!sapply(results, is.null)]),
    na_counter = na_counter
  )
}

# 6. get dates for extraction
get_unique_extraction_dates <- function(..., date_columns) {
  dfs <- list(...)
  all_dates <- unlist(mapply(function(df, col) df[[col]], dfs, date_columns, SIMPLIFY = FALSE))
  unique_dates <- sort(unique(as.Date(all_dates)))
  return(unique_dates)
}

# ---- 1. Setup ----
library(sf)
library(ncdf4)
library(future.apply)
plan(multisession, workers = 8) 


# Verify coordinate ranges
nc_first <- nc_open(list.files("gridExtractions/covariates/daily/", 
                               pattern = "rainfall.*\\.nc$", 
                               full.names = TRUE)[1])
x_coords <- ncvar_get(nc_first, "projection_x_coordinate")
y_coords <- ncvar_get(nc_first, "projection_y_coordinate")
nc_close(nc_first)

print(paste("Grid X range:", min(grid_df$x_bng), "-", max(grid_df$x_bng)))
print(paste("NetCDF X range:", min(x_coords), "-", max(x_coords)))
print(paste("Grid Y range:", min(grid_df$y_bng), "-", max(grid_df$y_bng)))
print(paste("NetCDF Y range:", min(y_coords), "-", max(y_coords)))

# ---- 2. Run Tests ----
# Test 1: Basic file filtering
test_date <- as.Date("2023-07-01")
test_files <- filter_nc_files_for_date_range(
  "gridExtractions/covariates/daily/", 
  test_date - 27, 
  test_date
)
print(test_files)

# Test file loading independently
test_files <- c(
  "gridExtractions/covariates/daily/rainfall_hadukgrid_uk_1km_day_20230601-20230630.nc",
  "gridExtractions/covariates/daily/rainfall_hadukgrid_uk_1km_day_20230701-20230731.nc"
)

test_data <- load_nc_data(test_files)
# Check what dates are actually loaded
print("Dates in second file:");print(test_data[[1]]$dates)
print("Dates in second file:");print(test_data[[2]]$dates)

test_date <- as.Date("2023-07-01")
start_date <- test_date - 27
message("Checking date range: ", start_date, " to ", test_date)

# Check which dates fall in this range
for (i in seq_along(test_data)) {
  in_range <- test_data[[i]]$dates >= start_date & test_data[[i]]$dates <= test_date
  message("File ", basename(test_data[[i]]$file), " has ", sum(in_range), " matching date(s)")
  print(test_data[[i]]$dates[in_range])
}


# Test 2: Single point extraction
# Use known good indices (try center of grid)
x_idx <- 90  # Midpoint of 180 grid cells
y_idx <- 135 # Midpoint of 290 grid cells

# Check raw values
for (i in seq_along(test_data)) {
  date_indices <- which(test_data[[i]]$dates >= start_date & test_data[[i]]$dates <= test_date)
  if (length(date_indices) > 0) {
    vals <- test_data[[i]]$rainfall[x_idx, y_idx, date_indices]
    message("\nFile: ", basename(test_data[[i]]$file))
    print(data.frame(
      Date = test_data[[i]]$dates[date_indices],
      Rainfall = vals
    ))
  }
}

# Test 3: Parallel extraction
test_dates <- as.Date(c("2023-06-15", "2023-07-15"))
test_results <- extract_grid_rainfall_ts_optimized(
  grid_df[1:5,], 
  test_dates, 
  "gridExtractions/covariates/daily/"
)
print(test_results)

# Test 4: highland test
# Approximate Fort William area in BNG
highland_point <- data.frame(
  x_bng = 210000,
  y_bng = 770000
)

# Find grid indices
nc_first <- nc_open(nc_files[1])
x_coords <- ncvar_get(nc_first, "projection_x_coordinate")
y_coords <- ncvar_get(nc_first, "projection_y_coordinate")
nc_close(nc_first)

highland_idx <- find_nearest_grid_bng(
  highland_point$x_bng,
  highland_point$y_bng,
  x_coords,
  y_coords
)

# Test extraction
test_rain <- get_28d_cumulative_preloaded(
  as.Date("2023-07-01"),
  highland_idx[1],
  highland_idx[2],
  test_data
)
print(test_rain)

# ---- 4. Full Extraction ----
survey_df_new <- read_csv("01_data/csvs/survey_df.csv")
survey_df_new$Setup_date <- as.Date(survey_df_new$Setup_date, 
        format = "%d/%m/%Y")


cs_df_new <- read_csv("01_data/csvs/cs_df.csv")
cs_df_new <- cs_df_new %>% filter(Verified_mosquito == "Yes" & Species == "pipiens" & Stage == "Adult")

cs_df_new$Date_found <- as.Date(cs_df_new$Date_found, 
                                          format = "%d/%m/%Y")

dates_to_extract <- get_unique_extraction_dates(
  survey_df_new, cs_df_new, 
  date_columns = c("Setup_date", "Date_found")
)
dates_to_extract[1]


# rainfall_grid_ts <- extract_grid_rainfall_ts_optimized(
#   grid_df,
#   dates_to_extract,
#   "gridExtractions/covariates/daily/"
# )

# Run with progress monitoring
library(progressr)
with_progress({
  rainfall_grid_ts <- extract_grid_rainfall_ts_optimized(grid_df, 
                                                         dates_to_extract,
                                                         "gridExtractions/covariates/daily/")
})
head(rainfall_grid_ts$rainfall_data)
rainfall_grid_ts$rainfall_data <- rainfall_grid_ts$rainfall_data %>%
  mutate(
    start_date = date - 27,
    end_date = date,
    .after = date  # places these right after the date column
  )

# see all date ranges for a specific grid_id
head(rainfall_grid_ts$rainfall_data %>% 
  filter(grid_id == 3972) %>%
  select(grid_id, date, rainfall_28d))

head(rainfall_grid_ts[["rainfall_data"]])

rainfall_grid <- rainfall_grid_ts$rainfall_data

#---- 5. grid_id check ------
table(unique(rainfall_grid_ts$rainfall_data$grid_id) %in% unique(grid_clipped$grid_id))

length(unique(rainfall_grid_ts[["rainfall_data"]][["grid_id"]]))

missing_grids <- setdiff(grid_clipped$grid_id, rainfall_grid_ts$rainfall_data$grid_id) # ...

original_ids <- unique(grid_clipped$grid_id)
result_ids <- unique(rainfall_grid_ts$rainfall_data$grid_id)

missing_ids <- setdiff(original_ids, result_ids)
if (length(missing_ids) > 0) {
  warning("Some grid IDs are missing from results: ", paste(missing_ids, collapse=", "))
}

# ---- 6. coordinate check ----
# Get original coordinates from grid_clipped
original_coords <- st_coordinates(st_centroid(grid_clipped)) |> 
  as.data.frame() |>
  cbind(grid_id = grid_clipped$grid_id)

# Merge with results
coord_check <- merge(
  original_coords,
  unique(rainfall_grid_ts$rainfall_data[, c("grid_id", "x_bng", "y_bng")]),
  by = "grid_id"
)

# Check for mismatches
mismatches <- coord_check[abs(coord_check$X - coord_check$x_bng) > 1 | 
                            abs(coord_check$Y - coord_check$y_bng) > 1, ]

# ---- 7. plot the grids -----

ggplot() +
  geom_rect(aes(xmin=min(x_coords), xmax=max(x_coords),
                ymin=min(y_coords), ymax=max(y_coords)),
            fill="lightblue", alpha=0.3) +
  geom_point(data=rainfall_grid_ts$rainfall_data$grid_id, aes(x_bng, y_bng), color="red", size=1) +
  geom_point(data=grid_df, aes(x_bng, y_bng), color="green", size=1) +
  labs(title="Grid Coverage (Green=Valid, Red=Outside NetCDF)") +
  theme_bw()


# ----- 8. save it -----
write.csv(rainfall_grid, 
          file = "01_data/csvs/rainfall_grid_new.csv",
          row.names = FALSE)


# ------ 9. plot it -----

library(readr)
library(ggplot2)
rainfall_grid <- read_csv("csvs/rainfall_grid.csv")

#> average cumulative rainfall across the study period
grid_avg_rain <- rainfall_grid %>%
  group_by(grid_id) %>%
  summarize(avg_rain = mean(rainfall_28d, na.rm = TRUE)) %>%
  ungroup()

# grid_clipped <- st_read("grids/grid_clipped_1km.shp")

grid_avg_sf <- grid_clipped %>%
  left_join(grid_avg_rain, by = "grid_id")

plot = ggplot(grid_avg_sf) +
  geom_sf(aes(fill = avg_rain), color = NA) +
  scale_fill_viridis_c(option = "plasma", name = "Avg 4wk\nCumulative Rainfall") +
  labs(title = "Average 4 week Cumulative Rainfall by Grid Cell",
       subtitle = "Across all dates during study period") +
  theme_void()

ggsave(filename = "maps/new_mean_rain_study_period.png",
       plot = ggplot(grid_avg_sf) +
         geom_sf(aes(fill = avg_rain), color = NA) +
         scale_fill_viridis_c(option = "plasma", name = "Avg 4wk\nCumulative Rainfall") +
         labs(title = "Average 4 week Cumulative Rainfall by Grid Cell",
              subtitle = "Across all dates during study period") +
         theme_void(),
       device = "png",
       dpi =  300,
       units = "in",
       width = 7,
       height = 8
)


map.temp<-grid_avg_sf

map.temp<-st_make_valid(map.temp)

map.temp<-st_transform(map.temp, 27700)

res <- 1000  # Grid cell size: 10km

map.2 <- st_make_grid(map.temp, cellsize = res, what = "centers")  # Generate grid
map.2 <- st_sf(geometry = map.2)  # Convert to sf object
map.2<-map.2[st_transform(grid_avg_sf, 27700),]


map.2$avg_rain <- st_join(map.2, map.temp)$avg_rain

map.r<-rast(data.frame(x=st_coordinates(map.2)[,1],
                       y=st_coordinates(map.2)[,2],
                       z=map.2$avg_rain),
            crs = st_crs(map.2)$wkt)

library(ggplot2)
library(cowplot)
library(mapview)
library(tidyterra)
library(terra)
avg_rain_map <- ggplot() +
  geom_spatraster(data = map.r, 
                  mask_projection = TRUE) +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid = element_blank()
  )   +
  scale_fill_viridis_c(
    option = "plasma",
    name = "",  
    na.value = "transparent"
  ) +
  theme_void()

ggsave(filename = "maps/new_mean_rain_study_period.png",
       plot = avg_rain_map,
       device = "png",
       dpi =  300,
       units = "in",
       width = 7,
       height = 8
)


########### ISLE OF ARRAN 


mask<-st_read("~/OneDrive - University of Glasgow/PhD/modelling/dataIntegration/MARS/arran_mask.shp")
plot(mask)
# mask<-st_transform(mask, st_crs(avg_rain_map))
# ARRAN<-avg_rain_map[mask,]

# Transform mask to match the raster CRS
mask <- st_transform(mask, st_crs(map.r))

# Now create the plot with the mask
avg_rain_map <- ggplot() +
  geom_spatraster(data = map.r) +
  geom_sf(data = mask, fill = NA, color = "black", linewidth = 1) +  # Add mask outline
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid = element_blank()
  ) +
  scale_fill_viridis_c(
    option = "plasma",
    name = "",  
    na.value = "transparent"
  ) +
  theme_void()

# If you want to mask the raster data to only show within the polygon:
# Create a masked version of the raster
# Crop to the mask extent first, then mask
map.r_cropped <- crop(map.r, mask)
map.r_masked <- mask(map.r_cropped, mask)


# Then plot the masked version
avg_rain_map_arran <- ggplot() +
  geom_spatraster(data = map.r_masked) +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid = element_blank()
  ) +
  scale_fill_viridis_c(
    option = "plasma",
    name = "",  
    na.value = "transparent"
  ) +
  theme_void()


ggsave(filename = "maps/ARRAN_mean_rain_study_period.png",
       plot = avg_rain_map_arran,
       device = "png",
       dpi =  300,
       units = "in",
       width = 7,
       height = 8
)


# Simple histogram using terra
hist(map.r_masked, 
     main = "Distribution of Rainfall Values\nin Isle of Arran",
     xlab = "Rainfall")

# Transform mask to match the CRS of map.2
mask <- st_transform(mask, st_crs(map.2))
map.2_masked <- map.2[mask, ]

# Create histogram from masked points
histogram_masked_arran <- ggplot(map.2_masked, aes(x = avg_rain)) +
  geom_histogram(bins = 30, fill = "gray", colour = "white") +
  labs(title = "Distribution of Rainfall Values\nin Isle of Arran",
       x = "Rainfall", 
       y = "Frequency") +
  theme_minimal()

print(histogram_masked_arran)

ggsave(filename = "maps/ARRAN_mean_rain_histogram_masked_arran.png",
       plot = histogram_masked_arran,
       device = "png",
       dpi =  300,
       units = "in",
       width = 7,
       height = 8
)

hist(map.r, 
     main = "Distribution of Rainfall Values\nacross Scotland",
     xlab = "Rainfall")
histogram_from_points <- ggplot(map.2, aes(x = avg_rain)) +
  geom_histogram(bins = 30, fill = "gray", colour = "white") +
  labs(title = "Distribution of Rainfall Values\nacross Scotland",
       x = "Rainfall", 
       y = "Frequency") +
  theme_minimal()

print(histogram_from_points)

ggsave(filename = "maps/SCOTLAND_mean_rain_histogram.png",
       plot = histogram_from_points,
       device = "png",
       dpi =  300,
       units = "in",
       width = 7,
       height = 8
)





