# ============================================================================
# PREPARE DATA FOR STAN (Marginalized Model)
# ============================================================================

library(tidyverse)
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
# load("01_data/processedCovariates/1km/loaded_data_1km.RData")
# load("01_data/processedCovariates/1km/env_land_data_1km.RData")
# load("01_data/processedCovariates/1km/env_temp_data_1km.RData")
# load("01_data/processedCovariates/1km/env_rain_data_1km_part1.RData")
# load("01_data/processedCovariates/1km/env_rain_data_1km_part2.RData")

# z_rain <- rbind(z_rain_part1, z_rain_part2)
# remove(z_rain_part1, z_rain_part2)
# gc()

load("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/AIM 1/pipiens-ISDM/01_data/loaded_data_cov_2point5km.RData")

# Organize land covariates
z_land <- as.data.frame(z_land_data_clean) %>%
  dplyr::arrange(grid_id) %>% # grid_id_area2km
  dplyr::select(starts_with("z_")) %>%
  dplyr::select(-c(z_saltwater_km2, z_urban_km2, z_suburban_km2, 
                   z_poi_count, z_reports, z_buildings_count, 
                   z_poi_log, z_livestock_density, z_freshwater_km2_log)) %>% # z_poi_count_grouped
  as.matrix()

z_poi <- as.data.frame(z_land_data) %>%
  dplyr::arrange(grid_id) %>% # grid_id_area2km
  pull(z_poi_log) # z_poi_count

z_reports <- as.data.frame(z_land_data) %>%
  dplyr::arrange(grid_id) %>% # grid_id_area2km
  pull(z_reports)

cat(sprintf("Full data: %d grids × %d dates\n", nrow(z_temp_clean), ncol(z_temp_clean)))


# ============================================================================
# COORDS FOR GP APPROX
# ============================================================================

# grid_coords_km <- coords_2point5km / 1000
# grid_coords_centered <- scale(grid_coords_km, center = TRUE, scale = FALSE)

# Domain boundaries (for basis function construction)
# L_x <- diff(range(grid_coords_centered[, 1]))
# L_y <- diff(range(grid_coords_centered[, 2]))

# ============================================================================
# SUBSET TO OBSERVED GRIDS ONLY (for computational efficiency)
# ============================================================================
cat("\nSubsetting to observed grids only...\n")

# Get unique grid IDs that have observations
observed_grid_ids <- sort(unique(c(survey_df$grid_id, cs_df$grid_id)))
n_grids_observed <- length(observed_grid_ids)

# Get all grid IDs from the full data
all_grid_ids <- as.data.frame(z_land_data) %>%
  dplyr::arrange(grid_id) %>% # grid_id_area2km
  pull(grid_id) # grid_id_area2km

cat(sprintf("  Reducing from %d to %d grids (%.1f%% reduction)\n", 
            length(all_grid_ids), n_grids_observed,
            100 * (1 - n_grids_observed/length(all_grid_ids))))

# Find positions of observed grids in the full data
obs_grid_positions <- match(observed_grid_ids, all_grid_ids)

# Subset all covariate matrices/vectors
z_temp_obs <- z_temp_clean[obs_grid_positions, , drop = FALSE]
z_rain_obs <- z_rain_clean[obs_grid_positions, , drop = FALSE]
z_land_obs <- z_land[obs_grid_positions, , drop = FALSE]
z_poi_obs <- z_poi[obs_grid_positions]
z_reports_obs <- z_reports[obs_grid_positions]
area_grid_obs <- area_grid[obs_grid_positions]

cat(sprintf("  Climate matrices: %d × %d\n", 
            nrow(z_temp_obs), ncol(z_temp_obs)))

# Convert BNG coordinates to km and center
stopifnot(!any(is.na(obs_grid_positions)))
coords_obs <- coords_2point5km[obs_grid_positions, , drop = FALSE]

#=============================================================================
# HSGP DATA PREPARATION - CORRECT IMPLEMENTATION
#=============================================================================

# 1. Convert to kilometers
grid_coords_km <- coords_obs / 1000

# 2. Center coordinates (this is what the paper requires!)
grid_coords_centered <- scale(grid_coords_km, center = TRUE, scale = FALSE)

# 3. Compute half-ranges
S_x <- max(abs(grid_coords_centered[, 1]))
S_y <- max(abs(grid_coords_centered[, 2]))

# 4. Set boundary factor
c <- 1.5  # Paper recommends 1.2-1.5

# 5. Compute boundaries
L_x <- c * S_x
L_y <- c * S_y

# 6. Choose number of basis functions
# Rule of thumb: M_sqrt should scale with (c * S) / rho 
# posterior rho is around 7: 1.5*837/7 = ~180
# With large domain and small rho, need more basis functions
M_sqrt <- 20  # Conservative choice (400 total basis functions)

# 7. Diagnostic output
cat("HSGP Configuration:\n")
cat("==================\n")
cat("Domain extent:\n")
cat("  X range:", range(grid_coords_centered[,1]), "km\n")
cat("  Y range:", range(grid_coords_centered[,2]), "km\n")
cat("\nHalf-ranges (S):\n")
cat("  S_x =", round(S_x, 2), "km\n")
cat("  S_y =", round(S_y, 2), "km\n")
cat("\nBoundary factor: c =", c, "\n")
cat("\nBoundaries (L = c × S):\n")
cat("  L_x =", round(L_x, 2), "km\n")
cat("  L_y =", round(L_y, 2), "km\n")
cat("\nBasis functions:\n")
cat("  M_sqrt =", M_sqrt, "\n")
cat("  Total M =", M_sqrt^2, "\n")
cat("\nExpected lengthscale range:\n")
cat("  Suggested: ρ > c*S/M_sqrt ≈", round(c*mean(c(S_x,S_y))/M_sqrt, 2), "km\n")

# 8. Verify coverage
if (!all(abs(grid_coords_centered[, 1]) <= L_x) || 
    !all(abs(grid_coords_centered[, 2]) <= L_y)) {
  stop("ERROR: Some coordinates fall outside [-L, L] domain!")
} else {
  cat("\n✓ All coordinates within domain boundaries\n")
}


# ============================================================================
# 3. HANDLE MISSING VALUES (Stan doesn't accept NAs)
# ============================================================================
cat("\nHandling missing values...\n")

impute_missing_values <- function(x, method = "mean", verbose = TRUE) {
  if (!any(is.na(x))) {
    if (verbose) cat("  No missing values\n")
    return(x)
  }
  
  na_count <- sum(is.na(x))
  na_pct <- 100 * na_count / length(x)
  
  if (method == "mean") {
    if (is.matrix(x)) {
      # For matrices, impute row-wise first, then column-wise, then overall
      for (i in 1:nrow(x)) {
        row_na <- is.na(x[i, ])
        if (any(row_na) && !all(row_na)) {
          x[i, row_na] <- mean(x[i, ], na.rm = TRUE)
        }
      }
      for (j in 1:ncol(x)) {
        col_na <- is.na(x[, j])
        if (any(col_na) && !all(col_na)) {
          x[col_na, j] <- mean(x[, j], na.rm = TRUE)
        }
      }
      # Final pass: any remaining NAs get overall mean
      if (any(is.na(x))) {
        x[is.na(x)] <- mean(x, na.rm = TRUE)
      }
    } else {
      # For vectors
      x[is.na(x)] <- mean(x, na.rm = TRUE)
    }
  } else if (method == "zero") {
    x[is.na(x)] <- 0
  }
  
  if (verbose) {
    cat(sprintf("  Imputed %d values (%.1f%%) using %s\n", 
                na_count, na_pct, method))
  }
  
  return(x)
}

# Impute missing values in all covariates
cat("Temperature:\n")
z_temp_clean <- impute_missing_values(z_temp_obs)

cat("Rainfall:\n")
z_rain_clean <- impute_missing_values(z_rain_obs)

cat("Land covariates:\n")
z_land_clean <- impute_missing_values(z_land_obs)

cat("Observation-level covariates:\n")
z_RH_clean <- impute_missing_values(survey_df$z_RH)
z_WS_rain_clean <- impute_missing_values(survey_df$z_WS_rain)

# Verify no NAs remain
stopifnot(!any(is.na(z_temp_clean)))
stopifnot(!any(is.na(z_rain_clean)))
stopifnot(!any(is.na(z_land_clean)))
stopifnot(!any(is.na(z_RH_clean)))
stopifnot(!any(is.na(z_WS_rain_clean)))

cat("✓ All missing values handled\n\n")

# ============================================================================
# 4. CREATE INDEX MAPPINGS
# ============================================================================
cat("Creating index mappings...\n")

# Map observed grid IDs to sequential indices 1:n_grids_observed
grid_id_to_idx <- setNames(seq_along(observed_grid_ids), observed_grid_ids)

# Update survey data with new grid indices
survey_df_stan <- survey_df %>%
  mutate(grid_idx = grid_id_to_idx[as.character(grid_id)]) %>%
  filter(!is.na(Culex))  # Remove any observations with missing response

# Update citizen science data with new grid indices
cs_presences_stan <- cs_presences %>%
  mutate(grid_idx = grid_id_to_idx[as.character(grid_id)])

cat(sprintf("  Survey observations: %d\n", nrow(survey_df_stan)))
cat(sprintf("  Presence-only records: %d\n", nrow(cs_presences_stan)))


# ============================================================================
# CHOOSING N MAX FOR MARGINALISATION
# ============================================================================

# epsilon <- 1e-8
# phi = 0.2
# lambda = 250
# N_max <- qnbinom(1 - epsilon, size = phi, mu = lambda)
# N_max
# 
# (50*74) + 100
# 
# N_max = ceiling(lambda + 6 * sqrt(lambda + lambda^2 / phi))
# N_max

# ============================================================================
# 5. PREPARE STAN DATA LIST
# ============================================================================
cat("\nPreparing Stan data list...\n")

# Note: In Stan model, we only need observed grids
# n_grids_total = n_grids_obs in the marginalized Stan model

stan_data <- list(
  # Dimensions
  n_grids_total = n_grids_observed,  # Only observed grids
  n_grids_obs = n_grids_observed,    # Same as total (all are observed)
  n_dates = ncol(z_temp_clean),
  n_obs_y = nrow(survey_df_stan),
  n_obs_po = nrow(cs_presences_stan),
  n_land_covs = ncol(z_land_clean),
  
  # Indices
  obs_grid_idx = 1:n_grids_observed,  # Identity mapping (all grids observed)
  site_to_grid = survey_df_stan$grid_idx,
  date_y = as.integer(survey_df_stan$date_1),
  trap_type = as.integer(survey_df_stan$trap_type_idx),
  po_grid_idx = as.integer(cs_presences_stan$grid_idx),
  date_po = as.integer(cs_presences_stan$date_1),
  
  # Observations
  y = as.integer(survey_df_stan$Culex),
  ones = rep(1L, nrow(cs_presences_stan)),
  
  # Covariates (all cleaned, no NAs)
  z_temp = z_temp_clean,
  z_rain = z_rain_clean,
  # z_ndvi = matrix(0, nrow = n_grids_observed, ncol = ncol(z_temp_clean)),  # Placeholder
  z_land = z_land_clean,
  z_poi = z_poi_obs,
  z_reports = z_reports_obs,
  z_RH = z_RH_clean,
  z_WS_rain = z_WS_rain_clean,
  
  # GP approximation
  M_sqrt = 20,  # Since M = 49, M_sqrt = 7
  grid_coords = grid_coords_centered,
  L_x = L_factor * L_x,  # Increase boundary slightly
  L_y = L_factor * L_y,
  
  # Other
  area_grid = area_grid_obs,
  CONSTANT = 10000,
  N_multiplier = 5,  # sum over N from y to y + 5*max(y) + 100; was 50
  grainsize = 10
)

# Final verification: no NAs in any data
check_nas <- function(data_list, name = "stan_data") {
  na_check <- sapply(data_list, function(x) {
    if (is.list(x)) return(NA)  # Skip nested lists
    any(is.na(x))
  })
  
  if (any(na_check, na.rm = TRUE)) {
    problem_vars <- names(data_list)[which(na_check)]
    stop(sprintf("NAs found in %s: %s", name, paste(problem_vars, collapse = ", ")))
  }
  
  cat(sprintf("✓ No NAs found in %s\n", name))
}

check_nas(stan_data)

cat("\nStan data prepared successfully!\n")
cat(sprintf("  Grids: %d\n", stan_data$n_grids_total))
cat(sprintf("  Dates: %d\n", stan_data$n_dates))
cat(sprintf("  Survey obs: %d\n", stan_data$n_obs_y))
cat(sprintf("  PO records: %d\n", stan_data$n_obs_po))
cat(sprintf("  Max N for marginalization: %d\n", 
            max(stan_data$y) * stan_data$N_multiplier + 100))


init_fun <- function() {
  list(
    # --------------------
    # Abundance model
    # --------------------
    beta0      = rnorm(1, 2, 0.1),
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
    delta0        = rnorm(1, 0, 0.2),
    delta_poi     = rnorm(1, 0, 0.2),
    delta_reports = rnorm(1, 0, 0.2),
    
    # --------------------
    # Dispersion
    # --------------------
    #phi = abs(rnorm(1, 2, 0.5)) + 0.5, # Ensure positive, add small offset
    phi = 1,
    
    # --------------------
    # GP; HSGP parameters
    # --------------------
    # test.stan
    # sigma_spatial = 0.01,  # Start TINY
    # rho_spatial = 50,      # km (safe for 2.5 km grid)
    # beta_gp_raw   = matrix(0, stan_data$M, n_dates)
    #test2.stan
    # rho_gp = abs(rnorm(1, 30, 5)) + 1.0,    # Add positive offset # Reasonable length-scale (km)
    # alpha_gp = abs(rnorm(1, 0.5, 0.1)) + 0.1, # Add positive offset alpha_gp = 0.5,  # Start smal
    alpha_gp = 0.5,
    rho_gp = 5,
    beta_gp_raw = rnorm(stan_data$M_sqrt^2, 0, 0.1)
  )
}



# ============================================================================
# 6. SAVE PREPARED DATA
# ============================================================================

saveRDS(stan_data, "01_data/processedCovariates/2.5km/stan_data_obsgrid_2.5km.rds")
saveRDS(observed_grid_ids, "01_data/processedCovariates/2.5km/observed_grid_ids_2.5km.rds")


save(stan_data, init_fun, observed_grid_ids, n_land_covs, n_dates, n_grids_observed, n_grids_obs,
     file = "01_data/processedCovariates/2.5km/stan_data_init_obsgrid_2.5km.RData")


# saveRDS(stan_data, "01_data/processedCovariates/1km/stan_data_obsgrid_1km.rds")
# saveRDS(observed_grid_ids, "01_data/processedCovariates/1km/observed_grid_ids_1km.rds")

# save(stan_data, init_fun, observed_grid_ids, n_land_covs, n_dates,
#      file = "01_data/processedCovariates/1km/stan_data_init_obsgrid_1km.RData")

