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
  dplyr::select(-c(z_saltwater_km2, z_urban_km2, z_suburban_km2, 
                   z_poi_count, z_reports, z_buildings_count, 
                   z_poi_log, z_livestock_density, z_freshwater_km2_log)) %>% 
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


#=============================================================================
# HSGP DATA PREPARATION - CORRECT IMPLEMENTATION
#=============================================================================

# 1. Convert to kilometers
grid_coords_km <- coords_2point5km / 1000

# 2. Center coordinates (this is what the paper requires!)
grid_coords_centered <- scale(grid_coords_km, center = TRUE, scale = FALSE)

# 3. Compute half-ranges
S_x <- max(abs(grid_coords_centered[, 1]))
S_y <- max(abs(grid_coords_centered[, 2]))

range(grid_coords_centered[,1])
range(grid_coords_centered[,2])

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

# =============================================================================
# NNGP PREPROCESSING: BUILD NEAREST-NEIGHBOUR STRUCTURE VIA spNNGP
# =============================================================================
# This script uses the official NNMatrix() wrapper (from the Stan NNGP case
# study by Lu Zhang) which internally calls spNNGP::spConjNNGP() to build the
# neighbour index efficiently.
# https://mc-stan.org/learn-stan/case-studies/nngp.html#latent-nngp-model-for-simulation-study
# https://github.com/LuZhangstat/NNGP_STAN/tree/master

library(spNNGP)
source("NNMatrix.R")   # Load the NNMatrix wrapper and helpers


# -----------------------------------------------------------------------------
# PARAMETERS
# -----------------------------------------------------------------------------
M <- 10   # Number of nearest neighbours. 10 is a reasonable default;
# increase to 15 for better approximation at the cost of speed.

# -----------------------------------------------------------------------------
# STEP 1: DEFINE COORDINATES
# -----------------------------------------------------------------------------
# coords_2point5km should be an [n_grids_total x 2] matrix of grid centroids
# in km (or whatever unit your distances are in — must match rho_gp prior).
# This is the same coordinate object you previously passed to Stan as
# grid_coords[n_grids_total, 2].

# Replace with your actual coordinates, e.g.:
# coords_2point5km <- st_coordinates(your_grid_sf) / 1000   # convert m -> km

N <- nrow(coords_2point5km)
cat("Number of grid cells (N):", N, "\n")
cat("Number of nearest neighbours (M):", M, "\n")

# -----------------------------------------------------------------------------
# STEP 2: BUILD NEAREST-NEIGHBOUR STRUCTURE
# -----------------------------------------------------------------------------
# NNMatrix() handles:
#   - Ordering of locations (default: x-axis ordering)
#   - kNN search (fast "cb" algorithm; use "brute" if on a regular grid
#     with many identical x-coordinates)
#   - Computing NN_ind, NN_dist, NN_distM

NN.matrix <- NNMatrix(
  coords         = coords_2point5km,
  n.neighbors    = M,
  n.omp.threads  = 2,
  search.type    = "cb"    # use "brute" if coords fall on a regular grid
)

str(NN.matrix)

# Quick sanity check: visualise neighbours of one location
# Check_Neighbors(NN.matrix$coords.ord, M, NN.matrix, ind = 500)

# -----------------------------------------------------------------------------
# STEP 3: REORDER ALL GRID-INDEXED DATA
# -----------------------------------------------------------------------------
# NN.matrix$ord is the permutation vector. Every array indexed by grid cell
# (1:n_grids_total) must be reordered using this before passing to Stan.
#
# inv_ord maps Stan's internal ordered indices back to original grid indices
# (useful for plotting posterior spatial fields).

ord     <- NN.matrix$ord
inv_ord <- order(ord)   # inverse permutation: ordered -> original

# Reorder all grid-indexed Stan data arrays:
z_land_ord      <- z_land[ord, ]
z_poi_ord       <- z_poi[ord]
z_reports_ord   <- z_reports[ord]
area_grid_ord   <- area_grid[ord]
z_temp_ord      <- z_temp_clean[ord, ]
z_rain_ord      <- z_rain_clean[ord, ]

# Remap observation-level grid indices.
# If site_to_grid[k] = g means "observation k is at original grid g",
# after reordering, that grid's new position is inv_ord[g]:
site_to_grid_ord <- inv_ord[site_to_grid]
po_grid_idx_ord  <- inv_ord[as.integer(cs_presences$grid_id_1)]
obs_grid_idx_ord <- inv_ord[obs_grid_idx]

# -----------------------------------------------------------------------------
# STEP 4: BUILD THE STAN DATA LIST
# -----------------------------------------------------------------------------
# Add the NNGP fields to existing stan_data.
# These REPLACE the old HSGP fields: M_sqrt, L_x, L_y, grid_coords.

stan_data <- list(
  # --- Unchanged dimensions ---
  n_grids_total = n_grids_total,  
  n_grids_obs = n_grids_obs,   # Only observed grids
  n_dates = ncol(z_temp_clean),
  n_obs_y = nrow(survey_df),
  n_obs_po = nrow(cs_presences),
  n_land_covs = ncol(z_land),

  # --- Reordered grid arrays ---
  z_temp        = z_temp_clean[ord, ],
  z_rain        = z_rain_clean[ord, ],
  z_land        = z_land[ord, ],
  z_poi         = z_poi[ord],
  z_reports     = z_reports[ord],
  area_grid     = area_grid[ord],

  # --- Remapped observation indices ---
  site_to_grid  = inv_ord[site_to_grid],
  po_grid_idx   = inv_ord[po_grid_idx],
  obs_grid_idx  = inv_ord[obs_grid_idx],
  date_y        = as.integer(survey_df$date_1),          # unchanged (not grid-indexed)
  date_po       = as.integer(cs_presences$date_1),
  trap_type     = as.integer(survey_df$trap_type_idx),
  z_RH = z_RH_clean,
  z_WS_rain = z_WS_rain_clean,
  y             = as.integer(survey_df$Culex),
  ones          = rep(1L, nrow(cs_presences)),
  survey_grid_idx = as.integer(survey_df$grid_id_1),
 

  # --- Other unchanged scalars ---
  CONSTANT      = 10000,
  N_multiplier  = 40,
  grainsize     = 10,

  # --- NEW: NNGP fields (replaces M_sqrt, L_x, L_y, grid_coords) ---
  M        = M,
  NN_ind   = NN.matrix$NN_ind,
  NN_dist  = NN.matrix$NN_dist,
  NN_distM = NN.matrix$NN_distM
)



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
  
  # HSGP parameters
  grid_coords = grid_coords_centered,
  L_x = result$L_x,
  L_y = result$L_y,
  M_sqrt = 50,
  grid_coords = grid_coords_centered,  # already centered, no more scaling!
  
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
    # GPparameters
    #     #marginalised_approxGP.stan // test2.stan
    # --------------------
    w_gp     = rep(0, N),
    #   alpha_gp = 1,
    #   rho_gp   = 5
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

