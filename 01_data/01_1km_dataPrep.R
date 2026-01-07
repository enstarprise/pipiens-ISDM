################################################################################
# Script: 01_dataPrep.R
# Description: Install and load required R packages, and prepare data for #
#             mosquito distribution using and integrated species distribution
#             models (ISDMs)
#             
# Author: shennice knight
# Date: 8 Dec 2025
# R Version: 4.5.1 (2025-06-13)
################################################################################

# Clear workspace
# rm(list = ls())

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

#===============================================================================
# PACKAGE MANAGEMENT
#===============================================================================

# Define required packages
required_packages <- c(
  # Data manipulation and I/O
  "readr",        # Fast reading of delimited files
  "writexl",      # Write Excel files
  "readxl",       # Read Excel files
  "tidyverse",    # Suite of data science packages
  "dplyr",        # Data manipulation
  "tidyr",        # Data tidying
  "lubridate",    # Date/time manipulation
  
  # Spatial analysis
  "sp",           # Classes for spatial data
  "sf",           # Simple features for spatial data
  "spdep",        # Spatial dependence (for neighborhood structures)
  
  # Modeling
  "nimble",       # BUGS-language flexible modeling
  "Matrix",       # Sparse and dense matrix classes
  
  # Utilities
  "pbapply"       # Progress bars for apply functions
)

#' Install missing packages
#'
#' @param packages Character vector of package names
#' @return NULL (installs packages as side effect)
install_missing <- function(packages) {
  installed <- packages %in% installed.packages()[, "Package"]
  
  if (any(!installed)) {
    cat("Installing missing packages:", 
        paste(packages[!installed], collapse = ", "), "\n")
    install.packages(packages[!installed])
  } else {
    cat("All required packages already installed.\n")
  }
}

# Install any missing packages
install_missing(required_packages)

# Suppress startup messages for cleaner output
suppressPackageStartupMessages({
  library(readr)
  library(writexl)
  library(readxl)
  library(tidyverse)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(nimble)
  library(Matrix)
  library(sp)
  library(sf)
  library(spdep)
  library(pbapply)
})


# ================================================================================
# HELPER FUNCTIONS
# ================================================================================

#' Handles dates (numeric serial dates), Unix timestamps,
#' POSIXct objects, and character strings.
#'
#' @param df Data frame containing date column
#' @param col_name Character string naming the date column (default: "date")
#' @return Data frame with converted date column
#' @examples
#' df <- data.frame(date = c(44562, 44563, 44564))
#' df <- fix_date_column(df, "date")
fix_date_column <- function(df, col_name = "date") {
  
  # Check if column exists
  if (!(col_name %in% names(df))) {
    warning(paste("Column", col_name, "not found in data frame"))
    return(df)
  }
  
  # Already a proper Date - no conversion needed
  if (inherits(df[[col_name]], "Date")) {
    message(paste("Column", col_name, "is already Date format"))
    return(df)
  }
  
  # Convert from POSIXct/POSIXt
  if (inherits(df[[col_name]], "POSIXct") || 
      inherits(df[[col_name]], "POSIXt")) {
    df[[col_name]] <- as.Date(df[[col_name]])
    message(paste("Converted", col_name, "from POSIXct to Date"))
    
    # Convert from numeric (Excel serial or Unix timestamp)
  } else if (is.numeric(df[[col_name]])) {
    max_val <- max(df[[col_name]], na.rm = TRUE)
    
    # Unix timestamp (seconds since 1970-01-01)
    if (max_val > 1000000) {
      df[[col_name]] <- as.Date(as.POSIXct(df[[col_name]], 
                                           origin = "1970-01-01"))
      message(paste("Converted", col_name, "from Unix timestamp to Date"))
      
      # Excel serial date (days since 1899-12-30)
    } else {
      df[[col_name]] <- as.Date(df[[col_name]], origin = "1899-12-30")
      message(paste("Converted", col_name, "from Excel serial to Date"))
    }
    
    # Convert from character
  } else {
    df[[col_name]] <- as.Date(as.character(df[[col_name]]))
    message(paste("Converted", col_name, "from character to Date"))
  }
  
  return(df)
}

#' Standardize (z-score) a numeric vector
#'
#' Centers and scales a vector to mean 0 and standard deviation 1
#'
#' @param x Numeric vector
#' @param na.rm Logical; remove NA values (default: TRUE)
#' @return Standardized numeric vector
standardize <- function(x, na.rm = TRUE) {
  (x - mean(x, na.rm = na.rm)) / sd(x, na.rm = na.rm)
}


#' Impute missing values using normal distribution
#'
#' Replaces NA values with random draws from a normal distribution
#' based on observed mean and standard deviation.
#'
#' @param x Numeric vector with missing values
#' @param min_sd Minimum standard deviation for imputation (default: 0.1)
#' @param verbose Logical; print imputation details (default: FALSE)
#' @return Numeric vector with imputed values
#' @examples
#' x <- c(1, 2, NA, 4, 5, NA)
#' x_imputed <- impute_missing_normally(x, verbose = TRUE)
impute_missing_normally <- function(x, min_sd = 0.1, verbose = FALSE) {
  
  # Input validation
  if (!is.numeric(x)) {
    stop("Input must be numeric vector")
  }
  
  na_count <- sum(is.na(x))
  
  # No missing values
  if (na_count == 0) {
    if (verbose) message("No missing values found")
    return(x)
  }
  
  # All values missing
  if (na_count == length(x)) {
    stop("All values are missing - cannot impute")
  }
  
  # Calculate imputation parameters
  x_mean <- mean(x, na.rm = TRUE)
  x_sd <- max(sd(x, na.rm = TRUE), min_sd)
  
  # Impute missing values
  missing_indices <- which(is.na(x))
  x[missing_indices] <- rnorm(na_count, mean = x_mean, sd = x_sd)
  
  # Report if verbose
  if (verbose) {
    message(sprintf(
      "Imputed %d missing values (%.1f%%) with mean = %.3f, sd = %.3f",
      na_count, 
      100 * na_count / length(x),
      x_mean, 
      x_sd
    ))
  }
  
  return(x)
}

#' Check for required input files
#'
#' Validates that all necessary data files exist before processing
#'
#' @param file_list Character vector of file paths
#' @return Logical; TRUE if all files exist, FALSE otherwise (with warnings)
check_required_files <- function(file_list) {
  
  missing_files <- !file.exists(file_list)
  
  if (any(missing_files)) {
    cat("\nERROR: Missing required files:\n")
    cat(paste("  -", file_list[missing_files], collapse = "\n"), "\n")
    return(FALSE)
  } else {
    cat("\nAll required files found:\n")
    cat(paste("  ✓", file_list, collapse = "\n"), "\n")
    return(TRUE)
  }
}

#' Create lookup table for mapping IDs to indices
#'
#' @param ids Vector of unique identifiers
#' @return Named integer vector (IDs as names, indices as values)
create_id_lookup <- function(ids) {
  unique_ids <- sort(unique(ids))
  setNames(seq_along(unique_ids), as.character(unique_ids))
}

#' Print data frame summary with NA counts
#'
#' @param df Data frame to summarize
#' @param name Optional name for the data frame
#' @return NULL (prints summary as side effect)
print_df_summary <- function(df, name = NULL) {
  
  if (!is.null(name)) {
    cat("\n", rep("=", 60), "\n", sep = "")
    cat("Summary:", name, "\n")
    cat(rep("=", 60), "\n", sep = "")
  }
  
  cat("Dimensions:", nrow(df), "rows x", ncol(df), "columns\n\n")
  
  cat("Column summary:\n")
  for (col in names(df)) {
    na_count <- sum(is.na(df[[col]]))
    na_pct <- 100 * na_count / nrow(df)
    cat(sprintf("  %-25s: %s, %d NAs (%.1f%%)\n", 
                col, class(df[[col]])[1], na_count, na_pct))
  }
  cat("\n")
}

# ============================================================================
# DATA LOADING SCRIPT - Updated for new covariate structure
# ============================================================================
data_dir <- "."

# Define all required input files
input_files <- list(
  survey = file.path(data_dir, "01_data/csvs/survey_df.csv"),
  cs = file.path(data_dir, "01_data/csvs/cs_df.csv"),
  grid = file.path(data_dir, "01_data/csvs/grid_clipped_1km_wkt.csv"),
  z_land = file.path(data_dir, "01_data/covariates/z_land_1km.csv"),
  z_climate = file.path(data_dir, "01_data/covariates/z_climate_1km.csv")
)

# Check that all files exist
if (!check_required_files(unlist(input_files))) {
  stop("Missing required input files. Please check file paths.")
}


# ============================================================================
# LOAD GRID AND OBSERVED DATA 
# ============================================================================

# SURVEY DATA
survey_df <- read_csv(input_files$survey)

# Code the trap types into integer categories; missing trap types are BGP
which(is.na(survey_df$Trap_type)) # BGP 
survey_df$Trap_type[is.na(survey_df$Trap_type)] <- "BGP"
survey_df <- survey_df %>%
  mutate(trap_type_idx = as.integer(factor(Trap_type)))
which(is.na(survey_df$trap_type_idx))

# Create DATE column 
survey_df$date <- survey_df$Setup_date
survey_df <- fix_date_column(survey_df, "date")

# Validate survey data
stopifnot("Culex" %in% names(survey_df))
stopifnot("date" %in% names(survey_df))

print_df_summary(survey_df, "Survey Data")
which(is.na(survey_df$date)) 
survey_df <- survey_df[!is.na(survey_df$date), ] # 2 obs with unknown set up dates etc

# Standardize microclimate variables
survey_df <- survey_df %>%
  mutate(
    z_RH = standardize(TT_mean_RH1),           # Relative humidity
    z_WS_rain = standardize(WS_total_rain)  #  rainfall
  )


# Extract for modeling
# z_RH <- as.numeric(survey_df$z_RH)
# z_WS_rain <- as.numeric(survey_df$z_WS_rain)

# CS DATA 
cs_df <- read_csv(input_files$cs)
cs_df$date <- cs_df$Date_found
cs_df <- fix_date_column(cs_df, "date")
cs_df$presence <- 1

cs_df <- cs_df %>% filter(Verified_mosquito == "Yes" & Species == "pipiens" & Stage == "Adult")

# Validate CS data
stopifnot("presence" %in% names(cs_df))
stopifnot("date" %in% names(cs_df))

cat(sprintf("Citizen science data loaded: %d observations from %s to %s\n",
            nrow(cs_df),
            min(cs_df$date),
            max(cs_df$date)))

print_df_summary(cs_df, "Citizen Science Data")



# ============================================================================
# LOAD GRID (WKT format)
# grid_sf <- st_read("01_data/grids/grid_clipped_area2km.gpkg")
# ============================================================================

grid <- read_csv(input_files$grid) %>%
  st_as_sf(wkt = "geometry", crs = 27700) 

# Ensure grid_id is numeric and rename to match covariate files
grid <- grid %>%
  mutate(
    grid_id = as.numeric(grid_id),
    # grid_id_area2km = grid_id  # Create matching column name
  )

# Validate grid data
stopifnot("grid_id" %in% names(grid))
stopifnot("area" %in% names(grid))

cat(sprintf("Grid data loaded: %d grid cells\n", nrow(grid)))
cat(sprintf("Total area: %.1f km²\n", sum(grid$area)))

print_df_summary(grid, "Grid Data")

# ============================================================================
# ADD GRID ID COLUMN TO THE DATAFRAMES
# ============================================================================

survey_df <- survey_df %>%
  drop_na(c("Trap_Lat", "Trap_Long")) # remove coordinates with missing value (1 obs.)

cs_df <- cs_df %>%
  drop_na(c("Latitude", "Longitude")) # remove coordinates with missing value (~6 obs.)

survey_sf <- st_as_sf(survey_df, coords = c("Trap_Long", "Trap_Lat"), crs = 27700)  # WGS84
cs_sf <- st_as_sf(cs_df, coords = c("Longitude", "Latitude"), crs = 27700)  # WGS84

survey_sf <- st_transform(survey_sf, crs = st_crs(grid))
cs_sf <- st_transform(cs_sf, crs = st_crs(grid))

survey_with_grid <- st_join(survey_sf, grid[c("grid_id", "area")])
cs_with_grid <- st_join(cs_sf, grid[c("grid_id", "area")])

survey_df <- survey_with_grid %>%
  st_drop_geometry()

cs_df <- cs_with_grid %>%
  st_drop_geometry()

# ============================================================================
# LOAD TIME-VARYING CLIMATE COVARIATES (z-scored)
# ============================================================================

z_climate_data <- read_csv(input_files$z_climate, show_col_types = FALSE)
z_climate_data$date <- as.Date(z_climate_data$date)

cat(sprintf("\nTime-varying climate covariates loaded: %d records\n", nrow(z_climate_data)))
cat(sprintf("Date range: %s to %s\n", 
            min(z_climate_data$date), 
            max(z_climate_data$date)))



# ============================================================================
# CONVERT CLIMATE DATA TO MATRICES (grids × dates)
# ============================================================================

cat("\nConverting climate data to matrix format...\n")

# Get sorted unique dates for column ordering
climate_dates_sorted <- sort(unique(z_climate_data$date))

# Create temperature matrix (rows = grids, columns = dates)
z_temp_min <- z_climate_data %>%
  dplyr::select(grid_id, date, z_temp_min) %>%
  pivot_wider(
    names_from = date,
    values_from = z_temp_min
  ) %>%
  arrange(grid_id) %>%
  dplyr::select(-grid_id) %>%
  as.matrix()

# Create rainfall matrix (rows = grids, columns = dates)
z_rain <- z_climate_data %>%
  dplyr::select(grid_id, date, z_rain) %>%
  pivot_wider(
    names_from = date,
    values_from = z_rain
  ) %>%
  arrange(grid_id) %>%
  dplyr::select(-grid_id) %>%
  as.matrix()

z_climate_cols <- grep("^z_", names(z_climate_data), value = TRUE)
cat(paste("  -", z_climate_cols, collapse = "\n"), "\n")

cat(sprintf("  Temperature matrix: %d grids × %d dates\n", nrow(z_temp_min), ncol(z_temp_min)))
cat(sprintf("  Rainfall matrix: %d grids × %d dates\n", nrow(z_rain), ncol(z_rain)))
cat(sprintf("  Temperature completeness: %.1f%%\n", 100 * mean(!is.na(z_temp_min))))
cat(sprintf("  Rainfall completeness: %.1f%%\n", 100 * mean(!is.na(z_rain))))

# ============================================================================
# LOAD STATIC LAND COVARIATES (z-scored)
# ============================================================================

z_land_data <- read_csv(input_files$z_land, show_col_types = FALSE)
# z_land_data <- subset(z_land_data, select = -c(z_saltwater_km2, z_urban_km2, z_suburban_km2))

cat(sprintf("Static land covariates loaded: %d grid cells\n", nrow(z_land_data)))
cat("Available z-scored land covariates:\n")
z_land_cols <- grep("^z_", names(z_land_data), value = TRUE)
cat(paste("  -", z_land_cols, collapse = "\n"), "\n")

# Summary statistics
cat("\nLand covariate ranges (z-scores):\n")
for (col in z_land_cols) {
  cat(sprintf("  %s: [%.2f, %.2f]\n", 
              col, 
              min(z_land_data[[col]], na.rm = TRUE),
              max(z_land_data[[col]], na.rm = TRUE)))
}

# ============================================================================
# MERGE COVARIATES WITH OBSERVATION DATA
# ============================================================================
# 
# # Merge static land covariates with survey data
# survey_df <- survey_df %>%
#   left_join(z_land_data, by = c("grid_id" = "grid_id_area2km"))
# 
# # Merge static land covariates with CS data
# cs_df <- cs_df %>%
#   left_join(z_land_data, by = c("grid_id" = "grid_id_area2km"))
# 
# # Merge time-varying climate covariates with survey data
# survey_df <- survey_df %>%
#   left_join(z_climate_data, by = c("grid_id" = "grid_id_area2km", "date" = "date"))
# 
# # Merge time-varying climate covariates with CS data
# cs_df <- cs_df %>%
#   left_join(z_climate_data, by = c("grid_id" = "grid_id_area2km", "date" = "date"))
# 
# # Check for missing covariates
# cat("\n=== Missing Covariate Check ===\n")
# cat(sprintf("Survey data - rows with missing land covariates: %d\n", 
#             sum(!complete.cases(survey_df[, z_land_cols]))))
# cat(sprintf("Survey data - rows with missing climate covariates: %d\n", 
#             sum(!complete.cases(survey_df[, z_climate_cols]))))
# cat(sprintf("CS data - rows with missing land covariates: %d\n", 
#             sum(!complete.cases(cs_df[, z_land_cols]))))
# cat(sprintf("CS data - rows with missing climate covariates: %d\n", 
#             sum(!complete.cases(cs_df[, z_climate_cols]))))
# 

# ============================================================================
# DATA VALIDATION SUMMARY
# ============================================================================

data_summary <- data.frame(
  Dataset = c("Survey", "Citizen Science", "Grid Cells",
              "Static Land Covariates", "Climate Covariates"),
  Records = c(nrow(survey_df), nrow(cs_df), nrow(grid),
              nrow(z_land_data), nrow(z_climate_data)),
  Variables = c(
    ncol(survey_df),
    ncol(cs_df),
    ncol(grid),
    length(z_land_cols),
    length(z_climate_cols)
  ),
  Status = rep("✓ Loaded", 5)
)

print(data_summary, row.names = FALSE)

# ==============================================================================
# 1. GRID AND DATA LOADING
# ==============================================================================
cat("Loading grid and climate data...\n")

# --- Grid data ---
lambda_grid <- grid[,c("grid_id", "geom", "area")]

# --- Date indexing ---
dates_to_keep <- union(unique(survey_df$date), unique(cs_df$date))
dates_to_keep <- as.Date(dates_to_keep)
all_dates_sorted <- sort(unique(dates_to_keep))
date_lookup <- setNames(seq_along(all_dates_sorted), as.character(all_dates_sorted))

cat("Data types fixed:\n")
cat("  grid grid_id:", class(grid$grid_id), "\n")
cat("  dates_to_keep:", class(dates_to_keep), "\n")

# ==============================================================================
# 3. CONSTANTS AND INDEXING
# ==============================================================================
cat("Setting up constants and indexing...\n")

n_dates <- length(dates_to_keep)
grid_levels <- unique(lambda_grid$grid_id)
grid_id_lookup <- setNames(seq_along(grid_levels), grid_levels)
n_grids_total <- length(grid_levels)

obs_grid_idx <- sort(unique(union(survey_df$grid_id, cs_df$grid_id)))
n_grids_obs <- length(obs_grid_idx)

area_km2 <- lambda_grid %>%
  distinct(grid_id, area) %>%
  arrange(grid_id) %>%
  pull(area)

obs_area <- full_join(survey_df[, c("grid_id", "area")],
                      cs_df[, c("grid_id", "area")],
                      by = c("grid_id", "area")) %>%
  distinct()

area_grid_mat <- matrix(area_km2, nrow = n_grids_total, ncol = n_dates)


obs_grid_ids <- unique(c(survey_df$grid_id, cs_df$grid_id))
n_grids_obs <- length(obs_grid_ids)

obs_grid_idx <- match(obs_grid_ids, grid$grid_id)
po_grid_idx <- match(cs_df$grid_id, grid$grid_id)
site_to_grid <- match(survey_df$grid_id, obs_grid_ids)

obs_area_values <- grid$area[match(obs_grid_ids, grid$grid_id)]
area_grid <- grid$area

# --- Map grid and date indices ---
survey_df$grid_id_1 <- grid_id_lookup[as.character(survey_df$grid_id)]
cs_df$grid_id_1 <- grid_id_lookup[as.character(cs_df$grid_id)]
survey_df$date_1 <- date_lookup[as.character(survey_df$date)]
cs_df$date_1 <- date_lookup[as.character(cs_df$date)]


# --- Create presence-only data with spatio-temporal indexing ---
# Following Fidino et al. approach:
# - Calculate lambda and thin_prob for ALL grids (for background)
# - Only use actual presence points in likelihood (not 0s and 1s)
# - Use ones trick with normalizing denominator (divided by n_obs_c)

cat("Preparing presence-only data for Bernoulli ones trick...\n")

# --- Extract only presence records ---
cs_presences <- cs_df %>%
  dplyr::select(grid_id, date, grid_id_1, date_1) %>%
  distinct() %>%
  arrange(grid_id_1, date_1)

# Bernoulli ones trick: observe 1 at each presence location
cs_presences$ones <- 1

# ==============================================================================
# 4. PREPARING FINAL COVARIATE MATRICES 
# ==============================================================================

cat("Preparing covariate matrices to be saved...\n")

## ensure grid and date ordering
z_climate_data <- z_climate_data %>%
  mutate(
    grid_id_1 = grid_id_lookup[as.character(grid_id)],
    date_1 = date_lookup[as.character(date)]
  )
# Temperature matrix
z_temp <- z_climate_data %>%
  dplyr::select(grid_id_1, date_1, z_temp_min) %>%
  pivot_wider(
    names_from = date_1,
    values_from = z_temp_min
  ) %>%
  arrange(grid_id_1) %>%
  dplyr::select(all_of(as.character(seq_len(n_dates)))) %>%  # ensures columns = date_1
  as.matrix()

# Rainfall matrix
z_rain <- z_climate_data %>%
  dplyr::select(grid_id_1, date_1, z_rain) %>%
  pivot_wider(
    names_from = date_1,
    values_from = z_rain
  ) %>%
  arrange(grid_id_1) %>%
  dplyr::select(all_of(as.character(seq_len(n_dates)))) %>%
  as.matrix()

# LAnd matrices
z_land_data <- z_land_data %>%
  arrange(grid_id) %>%
  # dplyr::select(starts_with("z_")) %>%
  # dplyr::select(-c(z_saltwater_km2, z_urban_km2, z_suburban_km2, z_poi_count, z_reports)) %>%
  as.matrix()

n_land_covs <- ncol(z_land_data[, setdiff(grep("^z_", colnames(z_land_data), value = TRUE), 
                                          c("z_saltwater_km2", "z_urban_km2", "z_suburban_km2", "z_poi_count", "z_reports"))])


# ==============================================================================
# SAVE DATA IN AN RDATA FILE
# ==============================================================================

# Save loaded data as RData for quick loading in subsequent scripts
save(
  survey_df, cs_df, cs_presences, grid,
  n_dates, n_grids_total, n_grids_obs, n_land_covs,
  obs_grid_idx, site_to_grid,
  
  area_grid,# area_grid_mat,
  
  dates_to_keep, grid_levels,
  date_lookup, grid_id_lookup,
  file = "01_data/loaded_data_1km.RData"
)

save(
  z_land_data,
  file = "01_data/env_land_data_1km.RData"
)

save(
  z_climate_data, 
  file = "01_data/env_climate_data_1km.RData"
)

save(
  z_temp, 
  file = "01_data/env_temp_data_1km.RData"
)
saveRDS(z_temp, "01_data/env_temp_data_1km.rds")

save(
  z_rain, 
  file = "01_data/env_rain_data_1km.RData"
)
saveRDS(z_rain, "01_data/env_rain_data_1km.rds")

# Save with maximum compression
saveRDS(z_temp, "01_data/env_temp_data_1km.rds", compress = "xz")
save(z_temp, file = "01_data/env_temp_data_1km.RData", compress = "xz")

# Split by rows (time periods or spatial points)
mid_point <- floor(nrow(z_rain) / 2)

z_rain_part1 <- z_rain[1:mid_point, ]
z_rain_part2 <- z_rain[(mid_point + 1):nrow(z_rain), ]


save(z_rain_part1, file = "01_data/env_rain_data_1km_part1.RData", compress = "xz")
save(z_rain_part2, file = "01_data/env_rain_data_1km_part2.RData", compress = "xz")
saveRDS(z_rain_part1, "01_data/env_rain_data_1km_part1.rds", compress = "xz")
saveRDS(z_rain_part2, "01_data/env_rain_data_1km_part2.rds", compress = "xz")

# On cluster, load and recombine:
# z_rain_part1 <- readRDS("env_rain_data_1km_part1.rds")
# z_rain_part2 <- readRDS("env_rain_data_1km_part2.rds")
# z_rain <- rbind(z_rain_part1, z_rain_part2)

cat("✓ All data loaded and saved to 'loaded_data_1km.RData'\n")

# HOW TO LOAD THE RDATA FILE
# loaded_data <- load("01_data/loaded_data_1km.RData")
