# ============================================================================
# PREPARE DATA FOR STAN (Marginalized Model)
# ============================================================================

library(tidyverse)
# ============================================================================
# 1. LOAD AND ORGANIZE COVARIATES
# ============================================================================
cat("Loading data...\n")
# load("01_data/processedCovariates/1km/loaded_data_1km.RData")
# load("01_data/processedCovariates/1km/env_land_data_1km.RData")
# load("01_data/processedCovariates/1km/env_temp_data_1km.RData")
# load("01_data/processedCovariates/1km/env_rain_data_1km_part1.RData")
# load("01_data/processedCovariates/1km/env_rain_data_1km_part2.RData")
# 
# z_rain <- rbind(z_rain_part1, z_rain_part2)
# remove(z_rain_part1, z_rain_part2)
# gc()

load("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/AIM 1/pipiens-ISDM/01_data/loaded_data_cov_2point5km.RData")

# Organize land covariates
z_land <- as.data.frame(z_land_data_clean) %>%
  dplyr::arrange(grid_id) %>% # grid_id_area2km
  dplyr::select(starts_with("z_")) %>%
  dplyr::select(-c(z_saltwater_km2, z_urban_km2, z_suburban_km2, z_poi_count, z_reports, z_buildings_count, z_poi_log)) %>% # z_poi_count_grouped
  as.matrix()

z_poi <- as.data.frame(z_land_data) %>%
  dplyr::arrange(grid_id) %>% # grid_id_area2km
  pull(z_poi_count) # z_poi_count_grouped

z_reports <- as.data.frame(z_land_data) %>%
  dplyr::arrange(grid_id) %>% # grid_id_area2km
  pull(z_reports)

cat(sprintf("Full data: %d grids × %d dates\n", nrow(z_temp_clean), ncol(z_temp_clean)))


# ============================================================================
# 2. HANDLE MISSING VALUES (Stan doesn't accept NAs)
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
# cat("Temperature:\n")
# z_temp_clean <- impute_missing_values(z_temp)
# 
# cat("Rainfall:\n")
# z_rain_clean <- impute_missing_values(z_rain)
# 
# cat("Land covariates:\n")
# z_land_clean <- impute_missing_values(z_land)

cat("Observation-level covariates:\n")
z_RH_clean <- impute_missing_values(survey_df$z_RH)
z_WS_rain_clean <- impute_missing_values(survey_df$z_WS_rain)

# Verify no NAs remain
stopifnot(!any(is.na(z_temp_clean)))
stopifnot(!any(is.na(z_rain_clean)))
stopifnot(!any(is.na(z_land)))
stopifnot(!any(is.na(z_RH_clean)))
stopifnot(!any(is.na(z_WS_rain_clean)))

cat("✓ All missing values handled\n\n")


# ============================================================================
# COORDS FOR GP APPROX
# ============================================================================

grid_coords_km <- coords_2point5km / 1000

# Domain boundaries (for basis function construction)
L_x <- diff(range(grid_coords_km[, 1]))
L_y <- diff(range(grid_coords_km[, 2]))

scale_coords_to_L <- function(coords, L) {
  # Center coordinates
  coords_centered <- scale(coords, center = TRUE, scale = FALSE)
  
  # Scale to [-L, L] range
  max_abs <- max(abs(coords_centered))
  if (max_abs > 0) {
    coords_scaled <- coords_centered * (L / max_abs)
  } else {
    coords_scaled <- coords_centered
  }
  return(coords_scaled)
}

# Use L_x and L_y from your calculation, can be adjusted
# Typically L is about 1.5 times the maximum coordinate value
L_factor <- 1.5  # Adjust as needed

# Scale coordinates
grid_coords_scaled <- scale_coords_to_L(grid_coords_km, L_factor * max(L_x, L_y))


# ============================================================================
# 4. CREATE INDEX MAPPINGS
# ============================================================================
cat("Creating index mappings...\n")

# Get unique grid IDs that have observations
# already loaded in
# observed_grid_ids <- sort(unique(c(survey_df$grid_id, cs_df$grid_id)))
# n_grids_obs <- length(observed_grid_ids)
# 
# # Get all grid IDs from the full data
# all_grid_ids <- as.data.frame(z_land_data) %>%
#   dplyr::arrange(grid_id) %>% # grid_id_area2km
#   pull(grid_id) # grid_id_area2km
# 
# n_grids_total <- length(all_grid_ids)
# 
# # Map observed grid IDs to sequential indices 1:n_grids_total
# ## already loaded in as the grid_id_lookup dict
# grid_id_to_idx <- setNames(seq_along(all_grid_ids), all_grid_ids)
# 
# # # Update the dfs with new grid indices
# ## already loaded in and created
# survey_df$grid_idx <- grid_id_lookup[as.character(survey_df$grid_id)]
# cs_presences$grid_idx <- grid_id_lookup[as.character(cs_presences$grid_id)]

cat(sprintf("  Survey observations: %d\n", nrow(survey_df)))
cat(sprintf("  Presence-only records: %d\n", nrow(cs_presences)))

# ============================================================================
# CHOOSING N MAX FOR MARGINALISATION
# ============================================================================

# epsilon <- 1e-8
# phi = 0.2
# lambda = 250
# N_max <- qnbinom(1 - epsilon, size = phi, mu = lambda)
# N_max
# (50*74) + 100
# N_max = ceiling(lambda + 6 * sqrt(lambda + lambda^2 / phi))
# N_max

# ============================================================================
# 5. PREPARE STAN DATA LIST
# ============================================================================
cat("\nPreparing Stan data list...\n")

stan_data <- list(
  # Dimensions
  n_grids_total = n_grids_total,  
  n_grids_obs = n_grids_obs,   # Only observed grids
  n_dates = ncol(z_temp_clean),
  n_obs_y = nrow(survey_df),
  n_obs_po = nrow(cs_presences),
  n_land_covs = ncol(z_land),
  
  # Indices
  obs_grid_idx = 1:n_grids_obs,  # Identity mapping (all grids)
  site_to_grid = survey_df$grid_id_1,
  date_y = as.integer(survey_df$date_1),
  trap_type = as.integer(survey_df$trap_type_idx),
  po_grid_idx = as.integer(cs_presences$grid_id_1),
  survey_grid_idx = as.integer(survey_df$grid_id_1),
  date_po = as.integer(cs_presences$date_1),

  # Observations
  y = as.integer(survey_df$Culex),
  ones = rep(1L, nrow(cs_presences)),
  
  # Covariates (all cleaned, no NAs)
  z_temp = z_temp_clean,
  z_rain = z_rain_clean,
  # z_ndvi = matrix(0, nrow = n_grids_observed, ncol = ncol(z_temp_clean)),  # Placeholder
  z_land = z_land,
  z_poi = z_poi,
  z_reports = z_reports,
  z_RH = z_RH_clean,
  z_WS_rain = z_WS_rain_clean,
  
  # GP approximation
  M_sqrt = 7,  # Since M = 49, M_sqrt = 7
  grid_coords = grid_coords_scaled,
  L_x = L_factor * L_x,  # Increase boundary slightly
  L_y = L_factor * L_y,
  
  # Other
  area_grid = area_grid,
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
    phi = abs(rnorm(1, 2, 0.5)) + 0.5, # Ensure positive, add small offset
    
    # --------------------
    # GP; HSGP parameters
    #     #marginalised_approxGP.stan // test2.stan
    # --------------------
    rho_gp = abs(rnorm(1, 30, 5)) + 1.0,    # Add positive offset # Reasonable length-scale (km)
    alpha_gp = abs(rnorm(1, 0.5, 0.1)) + 0.1, # Add positive offset alpha_gp = 0.5,  # Start small
    beta_gp_raw = rnorm(stan_data$M_sqrt^2, 0, 0.1)
  )
}




# ============================================================================
# 6. SAVE PREPARED DATA
# ============================================================================

saveRDS(stan_data, "01_data/processedCovariates/2.5km/stan_data_2.5km_full.rds")
saveRDS(obs_grid_idx, "01_data/processedCovariates/2.5km/observed_grid_ids_2.5km_full.rds")


save(stan_data, init_fun, obs_grid_idx, n_land_covs, n_dates,
     file = "01_data/processedCovariates/2.5km/stan_data_init_2.5km_full.RData")

