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
  "cmdstanr",       #cmdstanr
  "Matrix",       # Sparse and dense matrix classes
  
  # Utilities
  "pbapply"       # Progress bars for apply functions
)

# # Install any missing packages
# install_missing(required_packages)

# Suppress startup messages for cleaner output
suppressPackageStartupMessages({
  library(readr)
  library(writexl)
  library(readxl)
  library(tidyverse)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(cmdstanr)
  library(Matrix)
  library(sp)
  library(sf)
  library(spdep)
  library(pbapply)
})

# ============================================================================
# 1. LOAD AND ORGANIZE COVARIATES
# ============================================================================
cat("Loading data...\n")
load("01_data/processedCovariates/1km/loaded_data_1km.RData")
load("01_data/processedCovariates/1km/env_land_data_1km.RData")
load("01_data/processedCovariates/1km/env_temp_data_1km.RData")
load("01_data/processedCovariates/1km/env_rain_data_1km_part1.RData")
load("01_data/processedCovariates/1km/env_rain_data_1km_part2.RData")

z_rain <- rbind(z_rain_part1, z_rain_part2)
remove(z_rain_part1, z_rain_part2)
gc()

# Organize land covariates
z_land <- as.data.frame(z_land_data) %>%
  dplyr::arrange(grid_id) %>%
  dplyr::select(starts_with("z_")) %>%
  dplyr::select(-c(z_saltwater_km2, z_urban_km2, z_suburban_km2, z_poi_count, z_reports)) %>%
  as.matrix()

z_poi <- as.data.frame(z_land_data) %>%
  dplyr::arrange(grid_id) %>%
  pull(z_poi_count)

z_reports <- as.data.frame(z_land_data) %>%
  dplyr::arrange(grid_id) %>%
  pull(z_reports)

# ============================================================================
# 2. SUBSET TO OBSERVED GRIDS ONLY 
# ============================================================================
# Get unique grid IDs that have observations
observed_grid_ids <- sort(unique(c(survey_df$grid_id, cs_df$grid_id)))
n_grids_observed <- length(observed_grid_ids)

# Get all grid IDs from the full data
all_grid_ids <- as.data.frame(z_land_data) %>%
  dplyr::arrange(grid_id) %>%
  pull(grid_id)

cat(sprintf("  Reducing from %d to %d grids (%.1f%% reduction)\n", 
            length(all_grid_ids), n_grids_observed,
            100 * (1 - n_grids_observed/length(all_grid_ids))))

# Find positions of observed grids in the full data
obs_grid_positions <- match(observed_grid_ids, all_grid_ids)

# Subset all covariate matrices/vectors
z_temp_obs <- z_temp[obs_grid_positions, , drop = FALSE]
z_rain_obs <- z_rain[obs_grid_positions, , drop = FALSE]
z_land_obs <- z_land[obs_grid_positions, , drop = FALSE]
z_poi_obs <- z_poi[obs_grid_positions]
z_reports_obs <- z_reports[obs_grid_positions]
area_grid_obs <- area_grid[obs_grid_positions]

cat(sprintf("  Climate matrices: %d × %d\n", 
            nrow(z_temp_obs), ncol(z_temp_obs)))

# ============================================================================
# 3. PREPARE MISSING DATA FOR STAN IMPUTATION (MAR approach)
# ============================================================================

# Function to prepare missing data indicators for Stan
prepare_missing_data <- function(x, var_name = "variable") {
  na_count <- sum(is.na(x))
  na_pct <- 100 * na_count / length(x)
  
  cat(sprintf("%s:\n", var_name))
  cat(sprintf("  Missing: %d / %d (%.1f%%)\n", na_count, length(x), na_pct))
  
  # Create missing indicator (1 = missing, 0 = observed)
  missing_indicator <- is.na(x) * 1L
  
  # Create observed-only matrix (fill NA with placeholder 0)
  # These placeholder values won't be used where missing_indicator == 1
  x_obs <- x
  x_obs[is.na(x_obs)] <- 0
  
  list(
    obs = x_obs,
    missing = missing_indicator,
    n_missing = na_count
  )
}

# Prepare temperature with missing data tracking
temp_data <- prepare_missing_data(z_temp_obs, "Temperature")
z_temp_stan <- temp_data$obs
temp_missing <- temp_data$missing
n_temp_missing <- temp_data$n_missing

# Prepare rainfall with missing data tracking
rain_data <- prepare_missing_data(z_rain_obs, "Rainfall")
z_rain_stan <- rain_data$obs
rain_missing <- rain_data$missing
n_rain_missing <- rain_data$n_missing

# Land covariates impute regularly no MAR missingness based on other values
cat("\nLand covariates:\n")
if (any(is.na(z_land_obs))) {
  cat("  WARNING: Missing values found in land covariates\n")
  cat("  Imputing with column means...\n")
  for (j in 1:ncol(z_land_obs)) {
    col_na <- is.na(z_land_obs[, j])
    if (any(col_na)) {
      z_land_obs[col_na, j] <- mean(z_land_obs[, j], na.rm = TRUE)
      cat(sprintf("    Column %d: imputed %d values\n", j, sum(col_na)))
    }
  }
} else {
  cat("  ✓ No missing values (complete predictors for imputation)\n")
}
z_land_clean <- z_land_obs

## elevation would also have missing values in the full matrix, but this is a subset and no values are missing 

# Observation-level covariates
cat("\nObservation-level covariates:\n")
z_RH_clean <- survey_df$z_RH
z_WS_rain_clean <- survey_df$z_WS_rain

if (any(is.na(z_RH_clean))) {
  na_count <- sum(is.na(z_RH_clean))
  z_RH_clean[is.na(z_RH_clean)] <- mean(z_RH_clean, na.rm = TRUE)
  cat(sprintf("  z_RH: imputed %d values\n", na_count))
}

if (any(is.na(z_WS_rain_clean))) {
  na_count <- sum(is.na(z_WS_rain_clean))
  z_WS_rain_clean[is.na(z_WS_rain_clean)] <- mean(z_WS_rain_clean, na.rm = TRUE)
  cat(sprintf("  z_WS_rain: imputed %d values\n", na_count))
}

# Verify critical data structures
stopifnot(!any(is.na(z_land_clean)))  # Must be complete for imputation
stopifnot(!any(is.na(z_RH_clean)))
stopifnot(!any(is.na(z_WS_rain_clean)))

cat("\n✓ Missing data prepared for Stan imputation\n")
cat(sprintf("  Temperature: %d values to impute\n", n_temp_missing))
cat(sprintf("  Rainfall: %d values to impute\n", n_rain_missing))
cat(sprintf("  Total parameters to estimate: %d\n", 
            n_temp_missing + n_rain_missing))

# ============================================================================
# 4. CREATE INDEX MAPPINGS
# ============================================================================
cat("\n=== Creating index mappings ===\n")

# Map observed grid IDs to sequential indices 1:n_grids_observed
grid_id_to_idx <- setNames(seq_along(observed_grid_ids), observed_grid_ids)

# Update survey data with new grid indices
survey_df_stan <- survey_df %>%
  mutate(grid_idx = grid_id_to_idx[as.character(grid_id)]) %>%
  filter(!is.na(Culex))

# Update citizen science data with new grid indices
cs_presences_stan <- cs_presences %>%
  mutate(grid_idx = grid_id_to_idx[as.character(grid_id)])

cat(sprintf("  Survey observations: %d\n", nrow(survey_df_stan)))
cat(sprintf("  Presence-only records: %d\n", nrow(cs_presences_stan)))

# ============================================================================
# 5. CHOOSING N MAX FOR MARGINALISATION
# ============================================================================
cat("\n=== Setting N_max for marginalization ===\n")

epsilon <- 1e-8
phi <- 2.0  # Initial guess for overdispersion
lambda_max <- max(survey_df_stan$Culex) * 2  # Conservative estimate

# Method 1: Based on quantile
N_max_quantile <- qnbinom(1 - epsilon, size = phi, mu = lambda_max)

# Method 2: Based on mean + sd rule
N_max_rule <- ceiling(lambda_max + 6 * sqrt(lambda_max + lambda_max^2 / phi))

# Use the more conservative (larger) value
N_max <- max(N_max_quantile, N_max_rule, max(survey_df_stan$Culex) * 50 + 100)

cat(sprintf("  Method 1 (quantile): %d\n", N_max_quantile))
cat(sprintf("  Method 2 (mean + 6*SD): %d\n", N_max_rule))
cat(sprintf("  Using N_max = %d\n", N_max))

N_multiplier <- ceiling((N_max - max(survey_df_stan$Culex)) / max(survey_df_stan$Culex))
cat(sprintf("  N_multiplier = %d\n", N_multiplier))

N_multiplier <- 5

# ============================================================================
# 6. PREPARE STAN DATA LIST WITH IMPUTATION STRUCTURE
# ============================================================================
cat("\n=== Preparing Stan data list ===\n")

stan_data <- list(
  # Dimensions
  n_grids_total = n_grids_observed,
  n_grids_obs = n_grids_observed,
  n_dates = ncol(z_temp_stan),
  n_obs_y = nrow(survey_df_stan),
  n_obs_po = nrow(cs_presences_stan),
  n_land_covs = ncol(z_land_clean),
  
  # Indices
  obs_grid_idx = 1:n_grids_observed,
  site_to_grid = survey_df_stan$grid_idx,
  date_y = as.integer(survey_df_stan$date_1),
  trap_type = as.integer(survey_df_stan$trap_type_idx),
  po_grid_idx = as.integer(cs_presences_stan$grid_idx),
  date_po = as.integer(cs_presences_stan$date_1),
  
  # Observations
  y = as.integer(survey_df_stan$Culex),
  ones = rep(1L, nrow(cs_presences_stan)),
  
  # Climate covariates WITH MISSING DATA STRUCTURE
  z_temp_obs = z_temp_stan,
  z_rain_obs = z_rain_stan,
  temp_missing = temp_missing,
  rain_missing = rain_missing,
  n_temp_missing = n_temp_missing,
  n_rain_missing = n_rain_missing,
  
  # Complete covariates (predictors for imputation)
  z_land = z_land_clean,
  z_poi = z_poi_obs,
  z_reports = z_reports_obs,
  z_RH = z_RH_clean,
  z_WS_rain = z_WS_rain_clean,
  
  # Other
  area_grid = area_grid_obs,
  CONSTANT = 10000,
  N_multiplier = N_multiplier
)

temp_mean <- mean(z_temp_obs[!temp_missing])
temp_sd   <- sd(z_temp_obs[!temp_missing])

rain_mean <- mean(z_rain_obs[!rain_missing])
rain_sd   <- sd(z_rain_obs[!rain_missing])

init_fun <- function() {
  list(
    # --------------------
    # Abundance model
    # --------------------
    beta0      = rnorm(1, 0, 0.2),
    beta_temp  = rnorm(1, 0, 0.2),
    beta_rain  = rnorm(1, 0, 0.2),
    beta_land  = rnorm(n_land_covs, 0, 0.2),
    
    # --------------------
    # Detection model
    # --------------------
    alpha0         = qlogis(0.1),
    alpha_RH       = rnorm(1, 0, 0.2),
    alpha_WS_rain  = rnorm(1, 0, 0.2),
    alpha_trap_raw = rnorm(4, 0, 0.1),
    log_sigma_trap = log(0.5),
    
    # --------------------
    # Thinning model
    # --------------------
    delta0        = rnorm(1, -1, 0.2),
    delta_poi     = rnorm(1, 0, 0.2),
    delta_reports = rnorm(1, 0, 0.2),
    
    # --------------------
    # Dispersion
    # --------------------
    phi = abs(rnorm(1, 2, 0.5)) + 0.5,
    
    # --------------------
    # Imputation: temperature
    # --------------------
    gamma_temp0         = rnorm(1, 0, 0.2),
    gamma_temp_land     = rnorm(n_land_covs, 0, 0.2),
    gamma_temp_time_raw = rnorm(n_dates, 0, 0.1),
    sigma_temp          = 0.5, #  sigma_temp          = runif(1, 0.3, 1.2),
    sigma_temp_time     = 0.5, #  sigma_temp_time     = runif(1, 0.1, 0.8),
   
   
    
    # --------------------
    # Imputation: rain
    # --------------------
    gamma_rain0         = rnorm(1, 0, 0.2),
    gamma_rain_land     = rnorm(n_land_covs, 0, 0.2),
    gamma_rain_time_raw = rnorm(n_dates, 0, 0.1),
    sigma_rain          = 0.5, # sigma_rain          = runif(1, 0.4, 1.5),
    sigma_rain_time     = 0.5,# sigma_rain_time     = runif(1, 0.1, 1.0),
    
    nu_rain = 4,
    
    # --------------------
    # Imputed values
    # --------------------
    z_temp_imputed_raw = rnorm(n_temp_missing, temp_mean, temp_sd),
    z_rain_imputed_raw = rnorm(n_rain_missing, rain_mean, rain_sd)
  )
}



# ============================================================================
# 7. FINAL VERIFICATION
# ============================================================================

# Check for NAs in required complete data
check_nas <- function(data_list, name = "stan_data") {
  na_vars <- character(0)
  
  for (var_name in names(data_list)) {
    x <- data_list[[var_name]]
    if (is.list(x)) next  # Skip nested lists
    
    # Skip the _obs matrices (they have placeholders, not real NAs)
    if (grepl("_obs$", var_name)) next
    
    # Skip missing indicators (they're supposed to indicate missingness)
    if (grepl("missing$", var_name)) next
    
    if (any(is.na(x))) {
      na_vars <- c(na_vars, var_name)
    }
  }
  
  if (length(na_vars) > 0) {
    stop(sprintf("NAs found in %s: %s", name, paste(na_vars, collapse = ", ")))
  }
  
  cat(sprintf("✓ No unexpected NAs in %s\n", name))
}

check_nas(stan_data)

# Verify dimensions match
stopifnot(nrow(stan_data$z_temp_obs) == stan_data$n_grids_total)
stopifnot(ncol(stan_data$z_temp_obs) == stan_data$n_dates)
stopifnot(nrow(stan_data$z_rain_obs) == stan_data$n_grids_total)
stopifnot(ncol(stan_data$z_rain_obs) == stan_data$n_dates)
stopifnot(nrow(stan_data$z_land) == stan_data$n_grids_total)
stopifnot(length(stan_data$y) == stan_data$n_obs_y)


# ============================================================================
# 8. PRINT SUMMARY STATISTICS
# ============================================================================

cat(sprintf("  Grids: %d (%.1f%% of total)\n", 
            stan_data$n_grids_total,
            100 * stan_data$n_grids_total / length(all_grid_ids)))
cat(sprintf("  Dates: %d\n", stan_data$n_dates))
cat(sprintf("  Survey observations: %d\n", stan_data$n_obs_y))
cat(sprintf("  Presence-only records: %d\n", stan_data$n_obs_po))

cat(sprintf("\nMissing data to impute:\n"))
cat(sprintf("  Temperature: %d / %d (%.1f%%)\n", 
            stan_data$n_temp_missing,
            stan_data$n_grids_total * stan_data$n_dates,
            100 * stan_data$n_temp_missing / (stan_data$n_grids_total * stan_data$n_dates)))
cat(sprintf("  Rainfall: %d / %d (%.1f%%)\n", 
            stan_data$n_rain_missing,
            stan_data$n_grids_total * stan_data$n_dates,
            100 * stan_data$n_rain_missing / (stan_data$n_grids_total * stan_data$n_dates)))
cat(sprintf("  Total imputation parameters: %d\n",
            stan_data$n_temp_missing + stan_data$n_rain_missing))

cat("\nResponse variable (y):\n")
cat(sprintf("  Range: %d to %d\n", min(stan_data$y), max(stan_data$y)))
cat(sprintf("  Mean: %.2f\n", mean(stan_data$y)))
cat(sprintf("  Median: %.0f\n", median(stan_data$y)))
cat(sprintf("  Zeros: %d (%.1f%%)\n", 
            sum(stan_data$y == 0),
            100 * mean(stan_data$y == 0)))

cat("\nTrap types:\n")
print(table(stan_data$trap_type))

cat("\nCovariates:\n")
cat(sprintf("  Land covariates: %d (complete, used for imputation)\n", 
            stan_data$n_land_covs))
cat(sprintf("  Temperature range (observed): [%.2f, %.2f]\n", 
            min(stan_data$z_temp_obs[stan_data$temp_missing == 0]), 
            max(stan_data$z_temp_obs[stan_data$temp_missing == 0])))
cat(sprintf("  Rainfall range (observed): [%.2f, %.2f]\n", 
            min(stan_data$z_rain_obs[stan_data$rain_missing == 0]), 
            max(stan_data$z_rain_obs[stan_data$rain_missing == 0])))

cat(sprintf("\nMarginalization:\n"))
cat(sprintf("  N_max = %d\n", 
            max(stan_data$y) * stan_data$N_multiplier + 100))

# ============================================================================
# 9. SAVE PREPARED DATA
# ======= MUST UNCOMMENT ===========
# ============================================================================

save(stan_data, init_fun, observed_grid_ids, n_land_covs, n_dates,
     n_temp_missing, temp_mean, temp_sd,
     n_rain_missing, rain_mean, rain_sd,
     file = "01_data/processedCovariates/1km/stan_data_init_obsgrid_imputation_1km.RData")


# saveRDS(stan_data, "01_data/processedCovariates/1km/stan_data_imputation_1km.rds")
# saveRDS(init_fun, "01_data/processedCovariates/1km/stan_init_imputation_1km.rds")
# saveRDS(observed_grid_ids, "01_data/processedCovariates/1km/observed_grid_ids_1km.rds")

# # Save original matrices for comparison/validation
# saveRDS(list(
#   z_temp_original = z_temp_obs,
#   z_rain_original = z_rain_obs,
#   temp_missing_pattern = temp_missing,
#   rain_missing_pattern = rain_missing
# ), "01_data/processedCovariates/1km/original_climate_with_missingness_1km.rds")

cat("✓ Data saved successfully\n")
cat("\nFiles saved:\n")
cat("  - stan_data_imputation_1km.rds (Stan data list)\n")
cat("  - observed_grid_ids_1km.rds (grid ID mapping)\n")
cat("  - original_climate_with_missingness_1km.rds (for validation)\n")

cat("\n=== DATA PREPARATION COMPLETE ===\n")
cat("\nNext steps:\n")
cat("1. Fit the Stan model with imputation\n")
cat("2. Check convergence of imputation parameters (gamma_*, sigma_*)\n")
cat("3. Validate imputed values match observed distribution\n")
cat("4. Run posterior predictive checks\n")
cat("5. Compare results to complete-case analysis (sensitivity check)\n")