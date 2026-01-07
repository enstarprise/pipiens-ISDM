library(terra)
library(sf)
library(dplyr)
library(lubridate)



grid_clipped <- st_read("grids/grid_clipped_1km.gpkg")
# grid_clipped <- st_transform(grid_clipped, st_crs(scotland))
st_crs(grid_clipped)
library(terra)
library(sf)
library(dplyr)
library(lubridate)

library(terra)
library(sf)
library(dplyr)
library(lubridate)

# Function to process a single NDVI tiff file
process_single_ndvi <- function(tiff_path, grid_sf) {
  
  # Read NDVI raster
  ndvi_rast <- rast(tiff_path)
  
  # Rescale DN values
  ndvi_rast <- ndvi_rast * 0.004 - 0.08
  
  # Reproject grid to raster CRS
  grid_reproj <- st_transform(grid_sf, crs(ndvi_rast))
  
  # Extract mean NDVI for each polygon
  # Note: fun parameter now must be a function that handles na.rm internally
  ndvi_values <- terra::extract(ndvi_rast, vect(grid_reproj), 
                                fun = function(x) mean(x, na.rm = TRUE))
  
  # Parse date
  filename <- basename(tiff_path)
  date_str <- substr(filename, 15, 22)
  ndvi_date <- ymd(date_str)
  
  # Build output
  result <- data.frame(
    grid_id = grid_sf$grid_id,
    ndvi_date = ndvi_date,
    ndvi_mean = ndvi_values[, 2]   # Extracted value
  )
  
  return(result)
}

# Batch process multiple NDVI tiff files
batch_process_ndvi <- function(tiff_dir, grid_sf, pattern = "NDVI300.*\\.tif$", 
                               output_file = NULL) {
  
  # Get all tiff files matching the pattern
  tiff_files <- list.files(tiff_dir, pattern = pattern, full.names = TRUE)
  
  if (length(tiff_files) == 0) {
    stop("No NDVI tiff files found in the specified directory.")
  }
  
  cat(sprintf("Found %d NDVI files to process.\n", length(tiff_files)))
  
  # Process each file and combine results
  all_results <- lapply(seq_along(tiff_files), function(i) {
    cat(sprintf("Processing file %d of %d: %s\n", 
                i, length(tiff_files), basename(tiff_files[i])))
    
    tryCatch({
      process_single_ndvi(tiff_files[i], grid_sf)
    }, error = function(e) {
      warning(sprintf("Error processing %s: %s", 
                      basename(tiff_files[i]), e$message))
      return(NULL)
    })
  })
  
  # Remove NULL entries (failed processing)
  all_results <- Filter(Negate(is.null), all_results)
  
  # Combine all results into one data frame
  final_results <- dplyr::bind_rows(all_results)
  
  # Sort by date and grid_id - use .data pronoun to avoid NSE issues
  final_results <- final_results %>%
    dplyr::arrange(.data$ndvi_date, .data$grid_id)
  
  # Save to file if specified
  if (!is.null(output_file)) {
    write.csv(final_results, output_file, row.names = FALSE)
    cat(sprintf("\nResults saved to: %s\n", output_file))
  }
  
  cat(sprintf("\nProcessing complete. Total records: %d\n", nrow(final_results)))
  
  return(final_results)
}

# Usage:
results <- batch_process_ndvi(
  tiff_dir = "~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/vegetation/copernicusNDVI/NDVItiff",
  grid_sf = grid_clipped,
  output_file = "ndvi_by_grid_date.csv"
)
# View summary
head(results)
summary(results)

# Check for any missing values
sum(is.na(results$ndvi_mean))




##### ndvi for each date and grid_id
#### must first load lambda_grid_climate

dates_to_keep <- union(unique(survey_df$date), unique(cs_df$date)) 
dates_to_keep <- as.Date(dates_to_keep)
all_dates_sorted <- sort(unique(dates_to_keep))
date_lookup <- setNames(seq_along(all_dates_sorted), as.character(all_dates_sorted))

# rain and temp: varies by time
lambda_grid_climate <- expand.grid(grid_id = grid_clipped$grid_id, date = dates_to_keep)
lambda_grid_climate$date_index <- date_lookup[as.character(lambda_grid_climate$date)]


join_ndvi_to_climate <- function(climate_df, ndvi_results) {
  
  # Ensure date columns are Date type
  climate_df$date <- as.Date(climate_df$date)
  ndvi_results$date <- as.Date(ndvi_results$date)
  
  # Store original column names
  original_cols <- names(climate_df)
  
  # For each row in climate_df, find the closest NDVI date for that grid_id
  climate_with_ndvi <- climate_df %>%
    left_join(
      ndvi_results %>% select(grid_id, date, ndvi_mean),
      by = "grid_id",
      suffix = c("", "_ndvi"),
      relationship = "many-to-many"
    ) %>%
    mutate(
      date_diff = abs(as.numeric(date - date_ndvi))
    ) %>%
    group_by(across(all_of(original_cols))) %>%
    slice_min(date_diff, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    rename(ndvi_date = date_ndvi) %>%
    select(-date_diff)
  
  # Report matching statistics
  cat(sprintf("Total rows: %d\n", nrow(climate_with_ndvi)))
  cat(sprintf("Rows with NDVI data: %d (%.1f%%)\n", 
              sum(!is.na(climate_with_ndvi$ndvi_mean)),
              100 * sum(!is.na(climate_with_ndvi$ndvi_mean)) / nrow(climate_with_ndvi)))
  cat(sprintf("Rows without NDVI data: %d (%.1f%%)\n", 
              sum(is.na(climate_with_ndvi$ndvi_mean)),
              100 * sum(is.na(climate_with_ndvi$ndvi_mean)) / nrow(climate_with_ndvi)))
  
  # Show date difference summary for matched rows
  if (sum(!is.na(climate_with_ndvi$ndvi_mean)) > 0) {
    matched_rows <- climate_with_ndvi %>% filter(!is.na(ndvi_mean))
    date_diffs <- abs(as.numeric(matched_rows$date - matched_rows$ndvi_date))
    cat(sprintf("\nDate difference statistics (days):\n"))
    cat(sprintf("  Min: %d, Median: %d, Mean: %.1f, Max: %d\n",
                min(date_diffs), median(date_diffs), 
                mean(date_diffs), max(date_diffs)))
  }
  
  return(climate_with_ndvi)
}

test <- lambda_grid_climate
# Step 3: Join NDVI data to your climate dataframe
test <- join_ndvi_to_climate(
  climate_df = lambda_grid_climate,
  ndvi_results = results
)


write_csv(test, "csvs/ndvi_by_observed_date_2point5km.csv")

########

#### mapping

scotland <- st_read("gridExtractions/scotlandBoundary/scotland_boundary.shp") %>%
  vect()
plot(scotland)
st_crs(scotland)

grid_clipped <- st_read("grids/grid_clipped_2point5km.gpkg")
grid_clipped <- st_transform(grid_clipped, st_crs(scotland))
st_crs(grid_clipped) == st_crs(scotland)

ndvi <- read_csv("csvs/ndvi_by_observed_date_2point5km.csv")

grid_clipped <- grid_clipped %>%
  left_join(ndvi, by = "grid_id")


library(ggplot2)
ndvi_map <- ggplot() +
  geom_sf(data = grid_clipped, aes(fill = ndvi_mean)) +
  scale_fill_viridis_c(name = "Normalised Difference\nVegetation Index") +
  labs(title = "NDVI value") +
  theme_minimal()

ggsave(filename = "maps/ndvi_values_study_period.png",
       plot = ndvi_map,
       device = "png",
       dpi =  300,
       units = "in",
       width = 7,
       height = 8
)









