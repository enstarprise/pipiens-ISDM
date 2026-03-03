# =============================================================================
# NNGP PREPROCESSING: BUILD NEAREST-NEIGHBOUR STRUCTURE VIA spNNGP
# =============================================================================
# This script uses the official NNMatrix() wrapper (from the Stan NNGP case
# study by Lu Zhang) which internally calls spNNGP::spConjNNGP() to build the
# neighbour index efficiently.
# https://mc-stan.org/learn-stan/case-studies/nngp.html#latent-nngp-model-for-simulation-study
# https://github.com/LuZhangstat/NNGP_STAN/tree/master

# Source the helper file first:
#   source("NNMatrix.R")   # contains NNMatrix(), i_dist(), get_NN_*() helpers
#
# OUTPUT: NN.matrix, a list with:
#   $ord          : integer vector [N] — ordering of grid cells used internally
#   $coords.ord   : matrix [N x 2]    — coordinates in that order
#   $NN_ind       : matrix [N-1, M]   — neighbour indices (1-based, ordered)
#   $NN_dist      : matrix [N-1, M]   — distances to each neighbour
#   $NN_distM     : matrix [N-1, M*(M-1)/2] — inter-neighbour distances

# These feed directly into Stan. Every grid-indexed data array must be
# reordered by NN.matrix$ord before passing to Stan.
# =============================================================================

library(spNNGP)

source("NNMatrix.R")   # Load the NNMatrix wrapper and helpers
load("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/AIM 1/pipiens-ISDM/01_data/loaded_data_cov_2point5km.RData")
load("01_data/processedCovariates/2.5km/stan_data_init_obsgrid_2.5km.RData")

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
  stan_data$z_land      <- z_land_data_clean[ord, ]
  stan_data$z_poi       <- z_poi[ord]
  stan_data$z_reports   <- z_reports[ord]
  stan_data$area_grid   <- area_grid[ord]
  stan_data$z_temp      <- z_temp[ord, ]
  stan_data$z_rain      <- z_rain[ord, ]

# Remap observation-level grid indices.
# If site_to_grid[k] = g means "observation k is at original grid g",
# after reordering, that grid's new position is inv_ord[g]:
#   stan_data$site_to_grid <- inv_ord[site_to_grid]
#   stan_data$po_grid_idx  <- inv_ord[po_grid_idx]
#   stan_data$obs_grid_idx <- inv_ord[obs_grid_idx]

# -----------------------------------------------------------------------------
# STEP 4: BUILD THE STAN DATA LIST
# -----------------------------------------------------------------------------
# Add the NNGP fields to existing stan_data.
# These REPLACE the old HSGP fields: M_sqrt, L_x, L_y, grid_coords.

# stan_data <- list(
#   # --- Unchanged dimensions ---
#   n_grids_total = N,
#   n_grids_obs   = ...,
#   n_dates       = ...,
#   n_obs_y       = ...,
#   n_obs_po      = ...,
#   n_land_covs   = ...,
#
#   # --- Reordered grid arrays ---
#   z_temp        = z_temp[ord, ],
#   z_rain        = z_rain[ord, ],
#   z_land        = z_land[ord, ],
#   z_poi         = z_poi[ord],
#   z_reports     = z_reports[ord],
#   area_grid     = area_grid[ord],
#
#   # --- Remapped observation indices ---
#   site_to_grid  = inv_ord[site_to_grid],
#   po_grid_idx   = inv_ord[po_grid_idx],
#   obs_grid_idx  = inv_ord[obs_grid_idx],
#   date_y        = date_y,          # unchanged (not grid-indexed)
#   date_po       = date_po,
#   trap_type     = trap_type,
#   z_RH          = z_RH,
#   z_WS_rain     = z_WS_rain,
#   y             = y,
#   ones          = ones,
#
#   # --- Other unchanged scalars ---
#   CONSTANT      = CONSTANT,
#   N_multiplier  = N_multiplier,
#   grainsize     = grainsize,
#
#   # --- NEW: NNGP fields (replaces M_sqrt, L_x, L_y, grid_coords) ---
#   M        = M,
#   NN_ind   = NN.matrix$NN_ind,
#   NN_dist  = NN.matrix$NN_dist,
#   NN_distM = NN.matrix$NN_distM
# )

# -----------------------------------------------------------------------------
# STEP 5: INITIAL VALUES (RECOMMENDED)
# -----------------------------------------------------------------------------
# Initialising w_gp to zero avoids large gradients at the start of sampling.

# n_chains <- 4
# myinits <- lapply(1:n_chains, function(i) list(
#   w_gp     = rep(0, N),
#   alpha_gp = 1,
#   rho_gp   = 5
# ))

# -----------------------------------------------------------------------------
# STEP 6: POST-PROCESSING — MAP SPATIAL FIELD BACK TO ORIGINAL ORDER
# -----------------------------------------------------------------------------
# Stan returns spatial_field in NNGP ordering. To plot against original
# grid coordinates, reorder using inv_ord:
#
# spatial_field_draws <- extract(fit, "spatial_field")$spatial_field
# # [n_draws x N] — reorder columns back to original grid order
# spatial_field_original <- spatial_field_draws[, inv_ord]
