# ============================================================================
# PREPARE STAN DATA — HSGP / NNGP SWITCHABLE
# ============================================================================
# Set MODEL_TYPE at the top to switch between spatial GP implementations.
# Everything else is shared.
#
# MODEL_TYPE = "HSGP"           → original Hilbert Space GP
# MODEL_TYPE = "NNGP"           → standard latent NNGP
# MODEL_TYPE = "NNGP_centered"  → latent NNGP centred at intercept (better mixing)
# ============================================================================

MODEL_TYPE <- "NNGP_centered"   # <-- change this line only

stopifnot(MODEL_TYPE %in% c("HSGP", "NNGP", "NNGP_centered"))
cat(sprintf("\n=== Spatial model: %s ===\n\n", MODEL_TYPE))

library(tidyverse)

# ============================================================================
# 1. LOAD DATA
# ============================================================================
cat("Loading data...\n")

load("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/AIM 1/pipiens-ISDM/01_data/loaded_data_cov_2point5km.RData")

# Land covariates
z_land <- as.data.frame(z_land_data_clean) %>%
  dplyr::arrange(grid_id) %>%
  dplyr::select(starts_with("z_")) %>%
  dplyr::select(-c(z_saltwater_km2, z_urban_km2, z_suburban_km2,
                   z_poi_count, z_reports, z_buildings_count,
                   z_poi_log, z_livestock_density, z_freshwater_km2_log)) %>%
  as.matrix()

z_poi <- as.data.frame(z_land_data) %>%
  dplyr::arrange(grid_id) %>%
  pull(z_poi_count)

z_reports <- as.data.frame(z_land_data) %>%
  dplyr::arrange(grid_id) %>%
  pull(z_reports)

cat(sprintf("Full data: %d grids x %d dates\n", nrow(z_temp_clean), ncol(z_temp_clean)))


# ============================================================================
# 2. HANDLE MISSING VALUES
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

cat("Observation-level covariates:\n")
z_RH_clean      <- impute_missing_values(survey_df$z_RH)
z_WS_rain_clean <- impute_missing_values(survey_df$z_WS_rain)

stopifnot(!any(is.na(z_temp_clean)))
stopifnot(!any(is.na(z_rain_clean)))
stopifnot(!any(is.na(z_land)))
stopifnot(!any(is.na(z_RH_clean)))
stopifnot(!any(is.na(z_WS_rain_clean)))
cat("All missing values handled\n\n")


# ============================================================================
# 3. SHARED CONSTANTS
# ============================================================================

n_land_covs <- ncol(z_land)
n_dates     <- ncol(z_temp_clean)

CONSTANT     <- 10000
N_MULTIPLIER <- 40
GRAINSIZE    <- 10


# ============================================================================
# 4A. HSGP-SPECIFIC SETUP
# ============================================================================

if (MODEL_TYPE == "HSGP") {

  cat("--- HSGP setup ---\n")

  grid_coords_km       <- coords_2point5km / 1000
  grid_coords_centered <- scale(grid_coords_km, center = TRUE, scale = FALSE)

  S_x <- max(abs(grid_coords_centered[, 1]))
  S_y <- max(abs(grid_coords_centered[, 2]))
  c   <- 1.5
  L_x <- c * S_x
  L_y <- c * S_y
  M_sqrt <- 20

  cat(sprintf("  Domain: X [%.1f, %.1f] km | Y [%.1f, %.1f] km\n",
              min(grid_coords_centered[,1]), max(grid_coords_centered[,1]),
              min(grid_coords_centered[,2]), max(grid_coords_centered[,2])))
  cat(sprintf("  L_x = %.1f km | L_y = %.1f km\n", L_x, L_y))
  cat(sprintf("  M_sqrt = %d | Total basis functions = %d\n", M_sqrt, M_sqrt^2))

  if (!all(abs(grid_coords_centered[, 1]) <= L_x) ||
      !all(abs(grid_coords_centered[, 2]) <= L_y)) {
    stop("Some coordinates fall outside [-L, L] domain!")
  }
  cat("All coordinates within domain\n\n")

  # No reordering needed for HSGP
  ord     <- seq_len(n_grids_total)
  inv_ord <- seq_len(n_grids_total)

  spatial_data <- list(
    grid_coords = grid_coords_centered,
    L_x         = L_x,
    L_y         = L_y,
    M_sqrt      = M_sqrt
  )
}


# ============================================================================
# 4B. NNGP-SPECIFIC SETUP (shared between NNGP and NNGP_centered)
# ============================================================================

if (MODEL_TYPE %in% c("NNGP", "NNGP_centered")) {

  cat(sprintf("--- %s setup ---\n", MODEL_TYPE))

  library(spNNGP)
  source("NNMatrix.R")

  M_nngp <- 10   # Number of nearest neighbours

  cat(sprintf("  N grids = %d | M neighbours = %d\n", nrow(coords_2point5km), M_nngp))

  # NNMatrix() uses x-axis ordering by default.
  # Use search.type = "brute" if your grid has many identical x-coordinates.
  NN.matrix <- NNMatrix(
    coords        = coords_2point5km,
    n.neighbors   = M_nngp,
    n.omp.threads = 2,
    search.type   = "cb"
  )

  ord     <- NN.matrix$ord
  inv_ord <- order(ord)

  cat("  Neighbour structure built.\n\n")

  spatial_data <- list(
    M        = M_nngp,
    NN_ind   = NN.matrix$NN_ind,
    NN_dist  = NN.matrix$NN_dist,
    NN_distM = NN.matrix$NN_distM
  )
}


# ============================================================================
# 5. REORDER GRID-INDEXED DATA (identity for HSGP; NN ordering for NNGP)
# ============================================================================
# ord and inv_ord are set in section 4 for all model types.
# For HSGP, ord = 1:N so reordering is a no-op.

z_temp_ord    <- z_temp_clean[ord, ]
z_rain_ord    <- z_rain_clean[ord, ]
z_land_ord    <- z_land[ord, ]
z_poi_ord     <- z_poi[ord]
z_reports_ord <- z_reports[ord]
area_grid_ord <- area_grid[ord]

# Remap observation-level grid indices to new ordering
site_to_grid_ord <- inv_ord[site_to_grid]
po_grid_idx_ord  <- inv_ord[as.integer(cs_presences$grid_id_1)]
obs_grid_idx_ord <- inv_ord[obs_grid_idx]


# ============================================================================
# 6. BUILD STAN DATA LIST
# ============================================================================
cat("Building stan_data...\n")

stan_data <- c(
  list(
    # Dimensions
    n_grids_total = n_grids_total,
    n_grids_obs   = n_grids_obs,
    n_dates       = n_dates,
    n_obs_y       = nrow(survey_df),
    n_obs_po      = nrow(cs_presences),
    n_land_covs   = n_land_covs,

    # Indices (reordered)
    obs_grid_idx    = obs_grid_idx_ord,
    site_to_grid    = site_to_grid_ord,
    date_y          = as.integer(survey_df$date_1),
    trap_type       = as.integer(survey_df$trap_type_idx),
    po_grid_idx     = po_grid_idx_ord,
    survey_grid_idx = site_to_grid_ord,
    date_po         = as.integer(cs_presences$date_1),

    # Observations
    y    = as.integer(survey_df$Culex),
    ones = rep(1L, nrow(cs_presences)),

    # Covariates (reordered)
    z_temp    = z_temp_ord,
    z_rain    = z_rain_ord,
    z_land    = z_land_ord,
    z_poi     = z_poi_ord,
    z_reports = z_reports_ord,
    area_grid = area_grid_ord,

    # Observation-level covariates (not grid-indexed, unchanged)
    z_RH      = z_RH_clean,
    z_WS_rain = z_WS_rain_clean,

    # Scalars
    CONSTANT     = CONSTANT,
    N_multiplier = N_MULTIPLIER,
    grainsize    = GRAINSIZE
  ),
  spatial_data   # appends either HSGP or NNGP fields
)


# ============================================================================
# 7. INITIAL VALUES
# ============================================================================

init_fun <- function() {

  base_inits <- list(
    # Abundance
    beta0      = rnorm(1, 2, 0.1),
    beta_temp  = rnorm(1, 0, 0.2),
    beta_rain  = rnorm(1, 0, 0.2),
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
    phi = abs(rnorm(1, 2, 0.5)) + 0.5,

    # GP hyperparameters (same names for both HSGP and NNGP)
    alpha_gp = abs(rnorm(1, 0.5, 0.1)) + 0.1,
    rho_gp   = abs(rnorm(1, 30, 5)) + 1.0
  )

  # Spatial field initialisation varies by model type
  spatial_inits <- if (MODEL_TYPE == "HSGP") {
    list(
      beta_gp_raw = rnorm(stan_data$M_sqrt^2, 0, 0.1)
    )
  } else if (MODEL_TYPE == "NNGP") {
    list(
      w_gp = rep(0, n_grids_total)
    )
  } else {  # NNGP_centered
    list(
      w_b1 = rep(base_inits$beta0, n_grids_total)
      # w_b1 initialised at beta0 so the centred residual starts near zero
    )
  }

  c(base_inits, spatial_inits)
}


# ============================================================================
# 8. VALIDATION
# ============================================================================

check_nas <- function(data_list, name = "stan_data") {
  na_check <- sapply(data_list, function(x) {
    if (is.list(x)) return(NA)
    any(is.na(x))
  })
  if (any(na_check, na.rm = TRUE)) {
    problem_vars <- names(data_list)[which(na_check)]
    stop(sprintf("NAs found in %s: %s", name, paste(problem_vars, collapse = ", ")))
  }
  cat(sprintf("No NAs found in %s\n", name))
}

check_nas(stan_data)

cat("\nStan data ready:\n")
cat(sprintf("  Model type    : %s\n", MODEL_TYPE))
cat(sprintf("  Grids (total) : %d\n", stan_data$n_grids_total))
cat(sprintf("  Dates         : %d\n", stan_data$n_dates))
cat(sprintf("  Survey obs    : %d\n", stan_data$n_obs_y))
cat(sprintf("  PO records    : %d\n", stan_data$n_obs_po))
cat(sprintf("  Land covs     : %d\n", stan_data$n_land_covs))
if (MODEL_TYPE == "HSGP") {
  cat(sprintf("  Basis fns     : %d (M_sqrt = %d)\n",
              stan_data$M_sqrt^2, stan_data$M_sqrt))
} else {
  cat(sprintf("  NN neighbours : %d\n", stan_data$M))
}


# ============================================================================
# 9. SAVE
# ============================================================================

out_stem <- sprintf("01_data/processedCovariates/2.5km/stan_data_2.5km_%s", MODEL_TYPE)

saveRDS(stan_data, paste0(out_stem, ".rds"))

save(stan_data, init_fun, obs_grid_idx, ord, inv_ord,
     n_land_covs, n_dates, MODEL_TYPE,
     file = paste0(out_stem, ".RData"))

cat(sprintf("\nSaved to %s(.rds / .RData)\n", out_stem))
