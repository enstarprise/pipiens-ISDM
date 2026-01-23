# Clear workspace
# rm(list = ls())
gc()

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
  "nimble",      
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
  library(nimble)
  library(Matrix)
  library(sp)
  library(sf)
  library(spdep)
  library(pbapply)
})

# ============================================================================
# LOAD AND ORGANIZE COVARIATES
# ============================================================================
cat("Loading data...\n")
load("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/AIM 1/pipiens-ISDM/01_data/loaded_data_cov.RData")


# Organize land covariates
# FIX: USING THE ONE WHERE ALL NA'S ARE GIVEN THEIR MEAN
z_land <- as.data.frame(z_land_data_clean) %>%
  dplyr::arrange(grid_id) %>%
  dplyr::select(starts_with("z_")) %>%
  dplyr::select(-c(z_saltwater_km2, z_urban_km2, z_suburban_km2, z_poi_count_grouped, z_reports)) %>%
  as.matrix()

z_poi <- as.data.frame(z_land_data) %>%
  dplyr::arrange(grid_id) %>%
  pull(z_poi_count_grouped)

z_reports <- as.data.frame(z_land_data) %>%
  dplyr::arrange(grid_id) %>%
  pull(z_reports)

CONSTANT = 10000


# ==============================================================================
#  MODEL DATA AND CONSTANTS
# ==============================================================================
cat("Preparing model data and constants...\n")

data <- list(
  z_rain = z_rain_clean,
  z_rain2 = z_rain2,
  z_temp = z_temp_clean,
  z_land = z_land,
  z_poi = z_poi,
  # z_ndvi = z_land_data$z_ndvi,
  z_RH = as.numeric(survey_df$z_RH),
  z_WS_rain = as.numeric(survey_df$z_WS_rain),
  z_reports = z_reports,
  y = as.numeric(survey_df$Culex),
  ones = cs_presences$ones,  # vector of 1s for each presence
  # area_grid_mat = area_grid_mat,
  area_grid = area_grid
)

constants <- list(
  n_land_covs = n_land_covs,
  n_grids_total = n_grids_total,
  n_dates = n_dates,
  n_grids_obs = n_grids_obs,
  obs_grid_idx = obs_grid_idx,
  po_grid_idx = cs_presences$grid_id_1,  # grid index for each presence
  n_obs_y = nrow(survey_df),
  n_obs_po = nrow(cs_presences),  # number of presence-only records
  date_y = as.integer(survey_df$date_1),
  date_po = as.integer(cs_presences$date_1),  # date index for each presence
  site_to_grid = site_to_grid,
  trap_type = as.numeric(survey_df$trap_type_idx),
  CONSTANT = CONSTANT
)

cat("Data and constants prepared.\n")
cat("Note: Model will calculate lambda and thin_prob for ALL", 
    n_grids_total, "grids to compute background.\n\n")
cat("Data and constants prepared successfully.\n\n")

# ==============================================================================
# 8. INITIAL VALUES
# ==============================================================================
# --- Function to impute missing values ---
impute_missing_normally <- function(x, min_sd = 0.1, verbose = FALSE) {
  if (!is.numeric(x)) stop("Input must be numeric")
  na_count <- sum(is.na(x))
  if (na_count == 0) {
    if (verbose) message("No missing values found")
    return(x)
  }
  if (na_count == length(x)) stop("All values are missing - cannot impute")
  
  x_mean <- mean(x, na.rm = TRUE)
  x_sd <- max(sd(x, na.rm = TRUE), min_sd)
  missing_indices <- which(is.na(x))
  x[missing_indices] <- rnorm(na_count, mean = x_mean, sd = x_sd)
  
  if (verbose) {
    message(sprintf("Imputed %d missing values with mean = %.2f, sd = %.2f",
                    na_count, x_mean, x_sd))
  }
  return(x)
}


# --- Impute missing values ---
# z_ndvi_inits <- impute_missing_normally(z_ndvi, verbose = TRUE)
z_RH_inits <- impute_missing_normally(survey_df$z_RH, verbose = TRUE)
z_WS_rain_inits <- impute_missing_normally(survey_df$z_WS_rain, verbose = TRUE)
z_land_inits <- impute_missing_normally(z_land, verbose = TRUE)

# --- Initialize N matrix ---
N_init <- matrix(NA_integer_, nrow = n_grids_obs, ncol = n_dates)
N_init[,] <- pmax(1L, ceiling(mean(data$y, na.rm = TRUE) + 5L))

for(k in seq_len(nrow(survey_df))) {
  y_obs <- as.integer(survey_df$Culex[k])
  if(is.na(y_obs)) next
  s_local <- site_to_grid[k]
  d <- constants$date_y[k]
  N_init[s_local, d] <- max(N_init[s_local, d], y_obs + 5L)
}

# --- Initial values function ---
inits <- function() {
  
  log_lambda_init <- matrix(
    rnorm(n_grids_total * n_dates, mean = -2, sd = 0.5),
    nrow = n_grids_total, 
    ncol = n_dates
  )
  lambda_grid_init <- exp(log_lambda_init)
  
  # Initialize lambda_thinned to small positive values
  lambda_thinned_init <- matrix(
    exp(rnorm(n_grids_total * n_dates, -3, 0.5)),
    nrow = n_grids_total, 
    ncol = n_dates
  )
  
  # Initialize thinning probabilities
  logit_p_init <- rnorm(n_grids_total, mean = -1, sd = 0.5)
  p_init <- plogis(logit_p_init)
  
  # Calculate lambda_thinned consistently
  lambda_thinned_init <- matrix(NA_real_, nrow = n_grids_total, ncol = n_dates)
  for(g in 1:n_grids_total) {
    for(t in 1:n_dates) {
      lambda_thinned_init[g, t] <- lambda_grid_init[g, t] * 
        data$area_grid[g] * 
        p_init[g]
    }
  }
  
  # Calculate background to match model code (sum, not mean!)
  background_init <- numeric(n_dates)
  for(t in 1:n_dates) {
    background_init[t] <- sum(lambda_thinned_init[, t]) / constants$n_obs_po
  }
  
  list(
    # z_ndvi = z_ndvi_inits,
    z_RH = z_RH_inits,
    z_WS_rain = z_WS_rain_inits,
    z_land = z_land_inits,
    
    # Intensity model parameters
    beta0 = rnorm(1, -2, 0.3),
    beta_rain = rnorm(1, 0, 0.2),
    beta_temp = rnorm(1, 0, 0.2),
    beta_ndvi = rnorm(1, 0, 0.2),
    beta_buildings = rnorm(1, 0, 0.2),
    beta_land = rnorm(n_land_covs, 0, 0.2),
    
    phi = runif(1, 1, 3),
    
    # Detection parameters (only the ones with priors!)
    alpha0 = rnorm(1, qlogis(0.10), 0.2),
    alpha_RH = rnorm(1, 0, 0.1),
    alpha_WS_rain = rnorm(1, 0, 0.1),
    alpha_trap_raw = rnorm(4, 0, 0.2),  # Only 4 values
    sigma_trap = runif(1, 0.3, 0.8),
    # DON'T initialize alpha_trap (derived)
    # DON'T initialize p_trap (derived)
    
    # Thinning parameters
    delta0 = rnorm(1, -1, 0.3),
    delta_poi = rnorm(1, 0, 0.2),
    delta_reports = rnorm(1, 0, 0.2),
    
    # State variables (initialize in order of dependency)
    log_lambda = log_lambda_init,
    lambda_grid = lambda_grid_init,
    logit_p = logit_p_init,
    p = p_init,
    lambda_thinned = lambda_thinned_init,
    background = background_init,
    N = N_init
  )
}

initList <- inits()

