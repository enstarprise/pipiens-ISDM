# ============================================================================
# PREPARE STAN DATA — SWITCHABLE MODEL TYPE AND GRID SCOPE
# ============================================================================
# Three flags at the top control everything downstream.
#
# MODEL_TYPE:
#   "HSGP"          → Hilbert Space GP (original)
#   "NNGP"          → latent NNGP (standard)
#   "NNGP_centered" → latent NNGP centred at intercept (better mixing)
#   "nonspatial"    → no spatial random effect (baseline for LOO comparison)
#
# GRID_SCOPE:
#   "all"      → all n_grids_total grid cells (full spatial coverage)
#   "observed" → only grids with at least one survey or PO observation
#                (faster; good for model testing/development)
#
# QUADRATIC:
#   "none"  → only z_temp and z_rain passed  (beta_temp, beta_rain)
#   "rain"  → z_rain and z_rain_sq passed    (adds beta_rain2)
#   "both"  → all four matrices passed       (adds beta_temp2, beta_rain2)
#
#   z_temp_sq = z_temp^2, z_rain_sq = z_rain^2 (squared z-scores).
#   Stan always receives z_temp and z_rain. Quadratic matrices are
#   ADDITIONAL fields — both linear and quadratic terms are in the model
#   simultaneously when active. Inactive quadratic matrices are passed
#   as zeros so beta_temp2/beta_rain2 sample from their prior only.
#
# NOTE: NNGP with GRID_SCOPE = "all" and >10k grids will be slow to sample.
#       Use "observed" for development, "all" for final runs.
# ============================================================================

MODEL_TYPE <- "NNGP_centered"   # "HSGP" | "NNGP" | "NNGP_centered"
GRID_SCOPE <- "observed"        # "all"  | "observed"
QUADRATIC  <- "none"            # "none" | "rain" | "both"

stopifnot(MODEL_TYPE %in% c("HSGP", "NNGP", "NNGP_centered", "nonspatial"))
stopifnot(GRID_SCOPE %in% c("all", "observed"))
stopifnot(QUADRATIC  %in% c("none", "rain", "both"))
cat(sprintf("\n=== Spatial model: %s | Grid scope: %s | Quadratic: %s ===\n\n",
            MODEL_TYPE, GRID_SCOPE, QUADRATIC))

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
# 1. LOAD DATA
# ============================================================================
cat("Loading data...\n")

load("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/AIM 1/pipiens-ISDM/01_data/loaded_data_cov_2point5km.RData")

z_land <- as.data.frame(z_land_data_clean) %>%
  dplyr::arrange(grid_id) %>%
  dplyr::select(starts_with("z_")) %>%
  dplyr::select(-c(z_saltwater_km2, z_urban_km2, z_suburban_km2,
                   z_poi_count, z_reports, z_buildings_count,
                   z_poi_log, z_livestock_density, z_freshwater_km2_log)) %>%
  as.matrix()

z_poi <- as.data.frame(z_land_data) %>%
  dplyr::arrange(grid_id) %>%
  pull(z_poi_log)

z_reports <- as.data.frame(z_land_data) %>%
  dplyr::arrange(grid_id) %>%
  pull(z_reports)

all_grid_ids <- as.data.frame(z_land_data) %>%
  dplyr::arrange(grid_id) %>%
  pull(grid_id)

cat(sprintf("Full data: %d grids x %d dates\n", nrow(z_temp_clean), ncol(z_temp_clean)))

# ============================================================================
# 2. GRID SCOPE: SELECT ROWS TO USE
# ============================================================================

if (GRID_SCOPE == "observed") {

  cat("\nSubsetting to observed grids only...\n")

  observed_grid_ids <- sort(unique(c(survey_df$grid_id, cs_df$grid_id)))
  n_grids_total     <- length(observed_grid_ids)
  grid_positions    <- match(observed_grid_ids, all_grid_ids)

  stopifnot(!any(is.na(grid_positions)))

  cat(sprintf("  %d -> %d grids (%.1f%% reduction)\n",
              length(all_grid_ids), n_grids_total,
              100 * (1 - n_grids_total / length(all_grid_ids))))

  # Grid ID -> sequential index mapping for this scope
  grid_id_to_idx <- setNames(seq_along(observed_grid_ids), observed_grid_ids)

} else {  # "all"

  cat("\nUsing all grids...\n")

  observed_grid_ids <- all_grid_ids
  n_grids_total     <- length(all_grid_ids)
  grid_positions    <- seq_along(all_grid_ids)   # identity

  grid_id_to_idx <- setNames(seq_along(all_grid_ids), all_grid_ids)

  cat(sprintf("  %d grids\n", n_grids_total))
}

# Subset all grid-indexed arrays to the selected scope
coords_scope    <- coords_2point5km[grid_positions, , drop = FALSE]
z_temp_scope    <- z_temp_clean[grid_positions, , drop = FALSE]
z_rain_scope    <- z_rain_clean[grid_positions, , drop = FALSE]
z_land_scope    <- z_land[grid_positions, , drop = FALSE]
z_poi_scope     <- z_poi[grid_positions]
z_reports_scope <- z_reports[grid_positions]
area_grid_scope <- area_grid[grid_positions]

# Map observation-level grid IDs to sequential indices within selected scope
survey_df_stan <- survey_df %>%
  mutate(grid_idx = grid_id_to_idx[as.character(grid_id)]) %>%
  filter(!is.na(Culex), !is.na(grid_idx))

cs_presences_stan <- cs_presences %>%
  mutate(grid_idx = grid_id_to_idx[as.character(grid_id)]) %>%
  filter(!is.na(grid_idx))

cat(sprintf("  Survey obs: %d | PO records: %d\n",
            nrow(survey_df_stan), nrow(cs_presences_stan)))

# ============================================================================
# 3. MISSING VALUE IMPUTATION
# ============================================================================
cat("\nHandling missing values...\n")

impute_missing_values <- function(x, method = "mean", verbose = TRUE) {
  if (!any(is.na(x))) {
    if (verbose) cat("  No missing values\n")
    return(x)
  }
  na_count <- sum(is.na(x))
  na_pct   <- 100 * na_count / length(x)
  if (method == "mean") {
    if (is.matrix(x)) {
      for (i in 1:nrow(x)) {
        row_na <- is.na(x[i, ])
        if (any(row_na) && !all(row_na)) x[i, row_na] <- mean(x[i, ], na.rm = TRUE)
      }
      for (j in 1:ncol(x)) {
        col_na <- is.na(x[, j])
        if (any(col_na) && !all(col_na)) x[col_na, j] <- mean(x[, j], na.rm = TRUE)
      }
      if (any(is.na(x))) x[is.na(x)] <- mean(x, na.rm = TRUE)
    } else {
      x[is.na(x)] <- mean(x, na.rm = TRUE)
    }
  } else if (method == "zero") {
    x[is.na(x)] <- 0
  }
  if (verbose) cat(sprintf("  Imputed %d values (%.1f%%) using %s\n", na_count, na_pct, method))
  return(x)
}

cat("Temperature:\n");    z_temp_scope    <- impute_missing_values(z_temp_scope)
cat("Rainfall:\n");       z_rain_scope    <- impute_missing_values(z_rain_scope)
cat("Land covariates:\n"); z_land_scope   <- impute_missing_values(z_land_scope)
cat("Obs-level RH:\n");   z_RH_clean      <- impute_missing_values(survey_df$z_RH)
cat("Obs-level rain:\n"); z_WS_rain_clean <- impute_missing_values(survey_df$z_WS_rain)

stopifnot(!any(is.na(z_temp_scope)), !any(is.na(z_rain_scope)),
          !any(is.na(z_land_scope)), !any(is.na(z_RH_clean)),
          !any(is.na(z_WS_rain_clean)))
cat("All missing values handled\n\n")

# ============================================================================
# 3B. QUADRATIC TERMS
# ============================================================================
# Quadratic matrices are built by squaring the already-standardised z-scores.
# They are SEPARATE from the linear matrices — both terms enter the model
# simultaneously. z_temp_sq and z_rain_sq are ALWAYS built and passed to Stan
# as zero matrices when inactive, so beta_temp2/beta_rain2 are always declared
# in Stan but contribute nothing to the likelihood when their matrix is zero.
#
#   QUADRATIC = "none" → z_temp_sq = 0, z_rain_sq = 0
#   QUADRATIC = "rain" → z_temp_sq = 0, z_rain_sq = z_rain^2
#   QUADRATIC = "both" → z_temp_sq = z_temp^2, z_rain_sq = z_rain^2
#
# z_temp and z_rain are ALWAYS passed unchanged as the linear terms.

zeros_scope <- matrix(0, nrow = nrow(z_temp_scope), ncol = ncol(z_temp_scope))

z_temp_sq_scope <- if (QUADRATIC == "both") z_temp_scope^2 else zeros_scope
z_rain_sq_scope <- if (QUADRATIC %in% c("rain", "both")) z_rain_scope^2 else zeros_scope

cat(sprintf("Quadratic terms: %s\n\n", switch(QUADRATIC,
  "none" = "none  (z_temp_sq and z_rain_sq = zero matrices)",
  "rain" = "rain  (z_rain_sq = z_rain^2; z_temp_sq = 0)",
  "both" = "both  (z_temp_sq = z_temp^2, z_rain_sq = z_rain^2)"
)))

# ============================================================================
# 4. SHARED CONSTANTS
# ============================================================================

n_land_covs  <- ncol(z_land_scope)
n_dates      <- ncol(z_temp_scope)
n_grids_obs  <- n_grids_total   # for "observed": obs == total; for "all": same
                                 # (Stan model uses n_grids_total throughout)
CONSTANT     <- 10000
N_MULTIPLIER <- 40
GRAINSIZE    <- 10

# ============================================================================
# 5A. HSGP-SPECIFIC SETUP
# ============================================================================

if (MODEL_TYPE == "nonspatial") {

  cat("--- Non-spatial setup (no GP) ---

")
  ord          <- seq_len(n_grids_total)
  inv_ord      <- seq_len(n_grids_total)
  spatial_data <- list()   # no spatial fields passed to Stan

} else if (MODEL_TYPE == "HSGP") {

  cat("--- HSGP setup ---\n")

  grid_coords_km       <- coords_scope / 1000
  grid_coords_centered <- scale(grid_coords_km, center = TRUE, scale = FALSE)

  S_x    <- max(abs(grid_coords_centered[, 1]))
  S_y    <- max(abs(grid_coords_centered[, 2]))
  c_fac  <- 1.5
  L_x    <- c_fac * S_x
  L_y    <- c_fac * S_y
  M_sqrt <- 20

  cat(sprintf("  X [%.1f, %.1f] km | Y [%.1f, %.1f] km\n",
              min(grid_coords_centered[,1]), max(grid_coords_centered[,1]),
              min(grid_coords_centered[,2]), max(grid_coords_centered[,2])))
  cat(sprintf("  L_x = %.1f | L_y = %.1f | M_sqrt = %d (%d basis fns)\n",
              L_x, L_y, M_sqrt, M_sqrt^2))
  cat(sprintf("  Suggested min rho: %.1f km\n", c_fac * mean(c(S_x, S_y)) / M_sqrt))

  if (!all(abs(grid_coords_centered[, 1]) <= L_x) ||
      !all(abs(grid_coords_centered[, 2]) <= L_y))
    stop("Some coordinates fall outside [-L, L] domain!")
  cat("  All coordinates within domain\n\n")

  ord     <- seq_len(n_grids_total)   # no reordering for HSGP
  inv_ord <- seq_len(n_grids_total)

  spatial_data <- list(
    grid_coords = grid_coords_centered,
    L_x         = L_x,
    L_y         = L_y,
    M_sqrt      = M_sqrt
  )
}

# ============================================================================
# 5B. NNGP-SPECIFIC SETUP (shared for NNGP and NNGP_centered)
# ============================================================================

if (MODEL_TYPE %in% c("NNGP", "NNGP_centered")) {

  cat(sprintf("--- %s setup ---\n", MODEL_TYPE))

  library(spNNGP)
  source("NNMatrix.R")

  M_nngp <- 10

  cat(sprintf("  N = %d grids | M = %d neighbours\n", nrow(coords_scope), M_nngp))

  # "brute" is safer for regular grids with repeated x-coordinates
  search_type <- ifelse(GRID_SCOPE == "observed", "cb", "brute")
  cat(sprintf("  search.type = '%s'\n", search_type))

  NN.matrix <- NNMatrix(
    coords        = coords_scope,
    n.neighbors   = M_nngp,
    n.omp.threads = 2,
    search.type   = search_type
  )

  ord     <- NN.matrix$ord
  inv_ord <- order(ord)

  cat("  Neighbour structure built.\n")

  # ---------------------------------------------------------------------------
  # FIX: Replace zero-padding in NN_ind with valid index padding
  # ---------------------------------------------------------------------------
  # NNMatrix() pads rows for early locations (i < M) with 0 in NN_ind because
  # they have fewer than M predecessors. Stan validates ALL array values before
  # running any functions, so it rejects zeros in NN_ind[lower=1] immediately,
  # even though the Stan function's `dim` clipping means those padded slots are
  # never actually accessed in computation.
  #
  # Fix: replace zeros in each row with the last valid (non-zero) index in
  # that row. The padded values are never used (Stan slices 1:dim), so any
  # valid index >= 1 works — repeating the last real neighbour is conventional.

  NN_ind_fixed <- NN.matrix$NN_ind
  for (i in seq_len(nrow(NN_ind_fixed))) {
    row       <- NN_ind_fixed[i, ]
    valid_idx <- which(row > 0)
    if (length(valid_idx) < M_nngp) {
      last_valid        <- row[max(valid_idx)]
      row[row == 0]     <- last_valid
      NN_ind_fixed[i, ] <- row
    }
  }
  stopifnot(all(NN_ind_fixed >= 1))
  cat(sprintf("  NN_ind zero-padding fixed (%d early-location rows affected).\n\n",
              sum(rowSums(NN.matrix$NN_ind == 0) > 0)))

  spatial_data <- list(
    M        = M_nngp,
    NN_ind   = NN_ind_fixed,          # zeros replaced with last valid neighbour
    NN_dist  = NN.matrix$NN_dist,
    NN_distM = NN.matrix$NN_distM
  )
}

# ============================================================================
# 6. APPLY NNGP ORDERING TO ALL GRID-INDEXED ARRAYS
# ============================================================================
# For HSGP, ord = 1:N so this is a no-op.
# For NNGP, this reorders every grid-indexed array to match NN.matrix$ord,
# and remaps observation-level grid indices to the new positions.

z_temp_ord       <- z_temp_scope[ord, ]
z_rain_ord       <- z_rain_scope[ord, ]
z_temp_sq_ord    <- z_temp_sq_scope[ord, ]   # zero matrix when QUADRATIC != "both"
z_rain_sq_ord    <- z_rain_sq_scope[ord, ]   # zero matrix when QUADRATIC == "none"
z_land_ord       <- z_land_scope[ord, ]
z_poi_ord        <- z_poi_scope[ord]
z_reports_ord    <- z_reports_scope[ord]
area_grid_ord    <- area_grid_scope[ord]

# inv_ord[g] gives the new position of original grid g in the ordered system
site_to_grid_ord <- inv_ord[survey_df_stan$grid_idx]
po_grid_idx_ord  <- inv_ord[cs_presences_stan$grid_idx]
obs_grid_idx_ord <- inv_ord[seq_len(n_grids_total)]   # identity for obs-only scope

# ============================================================================
# 7. BUILD STAN DATA LIST
# ============================================================================
cat("Building stan_data...\n")

stan_data <- c(
  list(
    # Dimensions
    n_grids_total = n_grids_total,
    n_grids_obs   = n_grids_obs,
    n_dates       = n_dates,
    n_obs_y       = nrow(survey_df_stan),
    n_obs_po      = nrow(cs_presences_stan),
    n_land_covs   = n_land_covs,

    # Indices (reordered)
    obs_grid_idx    = obs_grid_idx_ord,
    site_to_grid    = site_to_grid_ord,
    date_y          = as.integer(survey_df_stan$date_1),
    trap_type       = as.integer(survey_df_stan$trap_type_idx),
    po_grid_idx     = po_grid_idx_ord,
    survey_grid_idx = site_to_grid_ord,
    date_po         = as.integer(cs_presences_stan$date_1),

    # Observations
    y    = as.integer(survey_df_stan$Culex),
    ones = rep(1L, nrow(cs_presences_stan)),

    # Covariates (grid-indexed, reordered)
    z_temp    = z_temp_ord,
    z_rain    = z_rain_ord,
    z_temp_sq = z_temp_sq_ord,   # zero matrix when QUADRATIC != "both"
    z_rain_sq = z_rain_sq_ord,   # zero matrix when QUADRATIC == "none"
    z_land    = z_land_ord,
    z_poi     = z_poi_ord,
    z_reports = z_reports_ord,
    area_grid = area_grid_ord,

    # Observation-level covariates (not grid-indexed, no reordering needed)
    z_RH      = z_RH_clean,
    z_WS_rain = z_WS_rain_clean,

    # Scalars
    CONSTANT     = CONSTANT,
    N_multiplier = N_MULTIPLIER,
    grainsize    = GRAINSIZE
  ),
  spatial_data   # appends HSGP or NNGP fields
)

# ============================================================================
# 8. INITIAL VALUES
# ============================================================================

init_fun <- function() {

  base_inits <- list(
    # Abundance
    beta0      = rnorm(1, 2, 0.1),
    beta_temp  = rnorm(1, 0, 0.2),
    beta_rain  = rnorm(1, 0, 0.2),
    beta_temp2 = rnorm(1, 0, 0.2),   # prior-only when z_temp_sq = 0
    beta_rain2 = rnorm(1, 0, 0.2),   # prior-only when z_rain_sq = 0
    beta_land  = rnorm(n_land_covs, 0, 0.2),

    # Detection
    alpha0         = qlogis(0.1),
    alpha_RH       = rnorm(1, 0, 0.2),
    alpha_WS_rain  = rnorm(1, 0, 0.2),
    alpha_trap_raw = rnorm(4, 0, 0.1),
    log_sigma_trap = log(0.5),

    # Thinning
    delta0        = rnorm(1, 0, 0.2),
    delta_poi     = rnorm(1, 0, 0.2),
    delta_reports = rnorm(1, 0, 0.2),

    # Dispersion
    phi = 1,

    # GP hyperparameters (same names for HSGP and NNGP)
    alpha_gp = 0.5,
    rho_gp   = 5
  )

  spatial_inits <- switch(MODEL_TYPE,
    "HSGP" = list(
      beta_gp_raw = rnorm(stan_data$M_sqrt^2, 0, 0.1)
    ),
    "NNGP" = list(
      w_gp = rep(0, n_grids_total)
    ),
    "NNGP_centered" = list(
      # Initialise w_b1 at beta0 so the zero-mean residual starts near zero,
      # avoiding a large initial gradient in the NNGP log-density evaluation
      w_b1 = rep(base_inits$beta0, n_grids_total)
    ),
    "nonspatial" = list()   # no spatial parameters to initialise
  )

  c(base_inits, spatial_inits)
}

# ============================================================================
# 9. VALIDATION
# ============================================================================

check_nas <- function(data_list, name = "stan_data") {
  na_check <- sapply(data_list, function(x) {
    if (is.list(x)) return(NA)
    any(is.na(x))
  })
  if (any(na_check, na.rm = TRUE))
    stop(sprintf("NAs in %s: %s", name,
                 paste(names(data_list)[which(na_check)], collapse = ", ")))
  cat(sprintf("No NAs in %s\n", name))
}

check_nas(stan_data)

cat("\nStan data ready:\n")
cat(sprintf("  Model type  : %s\n", MODEL_TYPE))
cat(sprintf("  Grid scope  : %s (%d grids)\n", GRID_SCOPE, n_grids_total))
cat(sprintf("  Quadratic   : %s\n", QUADRATIC))
cat(sprintf("  Dates       : %d\n", n_dates))
cat(sprintf("  Survey obs  : %d\n", nrow(survey_df_stan)))
cat(sprintf("  PO records  : %d\n", nrow(cs_presences_stan)))
cat(sprintf("  Land covs   : %d\n", n_land_covs))
if (MODEL_TYPE == "HSGP") {
  cat(sprintf("  Basis fns    : %d (M_sqrt = %d)\n", stan_data$M_sqrt^2, stan_data$M_sqrt))
} else if (MODEL_TYPE %in% c("NNGP", "NNGP_centered")) {
  cat(sprintf("  NN neighbours: %d\n", stan_data$M))
} else {
  cat("  Spatial      : none\n")
}

# ============================================================================
# 10. SAVE
# ============================================================================

out_stem <- sprintf(
  "01_data/processedCovariates/2.5km/stan_data_2.5km_%s_%s_%s",
  MODEL_TYPE, GRID_SCOPE, QUADRATIC
)

saveRDS(stan_data, paste0(out_stem, ".rds"))

save(stan_data, init_fun,
     observed_grid_ids, ord, inv_ord,
     n_land_covs, n_dates, n_grids_total,
     MODEL_TYPE, GRID_SCOPE, QUADRATIC,
     file = paste0(out_stem, ".RData"))

cat(sprintf("\nSaved: %s(.rds / .RData)\n", out_stem))
