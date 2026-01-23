################################################################################
# MEMORY-EFFICIENT MODEL SETUP
# The issue: 44,627 grids × dates is too large 
# Solution: Only model grids that have observations or are needed
################################################################################

# The problem:  model creates matrices for ALL 44,627 grids
# But only have observations in n_grids_obs grids
# NIMBLE loads all these matrices into memory during compilation

# ============================================================================
# 0. LIBRARY AND DATA LOADING
# ==============================================================================

setwd("/users/2601581k/pipiens-ISDM")

cat("Loading libraries...\n")

required_packages <- c(
  "readr", "writexl", "readxl", "tidyverse", "dplyr", "lubridate", "nimble",
  "sp", "Matrix", "pbapply"
)

# install.packages('spDataLarge', repos='https://nowosad.github.io/drat/', type='source')

install_missing <- function(packages) {
  installed <- packages %in% installed.packages()[,"Package"]
  if (any(!installed)) {
    install.packages(packages[!installed], repos = "https://cloud.r-project.org")
  }
}
install_missing(required_packages)

suppressPackageStartupMessages({
  library(readr)
  library(writexl)
  library(readxl)
  library(tidyverse)
  library(dplyr)
  library(lubridate)
  library(nimble)
  library(Matrix)
  library(sp)
  library(pbapply)
})


cat("Loading data...\n")
loaded_data <- load("01_data/loaded_data.RData")
# loaded_data <- load("loaded_data.RData")


loaded_data <- load("01_data/processedCovariates/area2km/loaded_data.RData")
load("01_data/processedCovariates/area2km/loaded_data.RData")
z_land_data <- load("01_data/processedCovariates/area2km/env_land_data.RData")
load("01_data/processedCovariates/area2km/env_land_data.RData")
z_temp <- load("01_data/processedCovariates/area2km/env_temp_data.RData")
load("01_data/processedCovariates/area2km/env_temp_data.RData")
z_rain <- load("01_data/processedCovariates/area2km/env_rain_data.RData")
load("01_data/processedCovariates/area2km/env_rain_data.RData")



# ============================================================================
# SOLUTION 1: Use only observed grids in the model
# ============================================================================

cat("Subsetting to observed grids only...\n")

# Get list of grids that actually have data
observed_grid_ids <- sort(unique(c(survey_df$grid_id, cs_df$grid_id)))
n_grids_model <- length(observed_grid_ids)

# Create mapping: original grid_id -> model grid index
grid_id_to_model_idx <- setNames(seq_along(observed_grid_ids), observed_grid_ids)

cat(sprintf("  Reducing from %d to %d grids (%.1f%% reduction)\n", 
            n_grids_total, n_grids_model, 
            100 * (1 - n_grids_model/n_grids_total)))

# Subset climate matrices to only observed grids
obs_grid_positions <- match(observed_grid_ids, grid_levels)

z_temp_subset <- z_temp[obs_grid_positions, , drop = FALSE]
z_rain_subset <- z_rain[obs_grid_positions, , drop = FALSE]

cat(sprintf("  Climate matrices reduced to: %d × %d\n", 
            nrow(z_temp_subset), ncol(z_temp_subset)))

# Create rainfall matrix (rows = grids, columns = dates)
# z_rain <- z_climate_data %>%
#   select(grid_id_area2km, date, z_rain) %>%
#   pivot_wider(
#     names_from = date,
#     values_from = z_rain
#   ) %>%
#   arrange(grid_id_area2km) %>%
#   select(-grid_id_area2km) %>%
#   as.matrix()

# Subset land covariates
z_land_subset <- as.data.frame(z_land_data) %>%
  filter(grid_id_area2km %in% observed_grid_ids) %>%
  arrange(grid_id_area2km) %>%
  dplyr::select(starts_with("z_")) %>%
  dplyr::select(-c(z_urban_km2, z_suburban_km2, z_poi_count_grouped, z_reports)) %>%
  as.matrix()


z_poi_subset <- as.data.frame(z_land_data) %>%
  filter(grid_id_area2km %in% observed_grid_ids) %>%
  arrange(grid_id_area2km) %>%
  pull(z_poi_count_grouped)

z_reports_subset <- as.data.frame(z_land_data) %>%
  filter(grid_id_area2km %in% observed_grid_ids) %>%
  arrange(grid_id_area2km) %>%
  pull(z_reports)

z_buildings_subset <- as.data.frame(z_land_data) %>%
  filter(grid_id_area2km %in% observed_grid_ids) %>%
  arrange(grid_id_area2km) %>%
  pull(z_buildings_count)

# Subset area
area_grid_subset <- area_grid[obs_grid_positions]
area_grid_mat_subset <- area_grid_mat[obs_grid_positions, , drop = FALSE]

# Update citizen science indices to new grid numbering
cs_presences_updated <- cs_presences %>%
  mutate(grid_model_idx = grid_id_to_model_idx[as.character(grid_id)])

# Update survey indices
survey_df_updated <- survey_df %>%
  mutate(grid_model_idx = grid_id_to_model_idx[as.character(grid_id)])

# Site to grid mapping (within observed grids only)
site_to_grid_updated <- match(survey_df$grid_id, observed_grid_ids)

### create 

dist_matrix <- as.matrix(dist(grid_coords))  # grid_coords is n_grids x 2 matrix

cat("✓ Data subsetted successfully\n\n")


# ============================================================================
# UPDATED DATA AND CONSTANTS LIST
# ============================================================================

cat("Creating updated data list...\n")

data <- list(
  z_temp = z_temp_subset,
  z_rain = z_rain_subset,
  z_land = z_land_subset,
  z_poi = z_poi_subset,
  z_reports = z_reports_subset,
  # z_buildings = z_buildings_subset,
  z_RH = as.numeric(survey_df$z_RH),
  z_WS_rain = as.numeric(survey_df$z_WS_rain),
  y = as.numeric(survey_df$Culex),
  ones = cs_presences_updated$ones,
  # area_grid_mat = area_grid_mat_subset,
  area_grid = area_grid_subset,
  dist_matrix = dist_matrix
)

cat("✓ Data list created\n")


cat("Creating updated constants...\n")

constants <- list(
  n_land_covs = ncol(z_land_subset),
  n_grids = n_grids_model,  # Only observed grids now
  n_dates = n_dates,
  n_obs_y = nrow(survey_df),
  n_obs_po = nrow(cs_presences_updated),
  
  # Indices (now all within the subset)
  po_grid_idx = cs_presences_updated$grid_model_idx,
  date_y = as.integer(survey_df$date_1),
  date_po = as.integer(cs_presences_updated$date_1),
  site_to_grid = site_to_grid_updated,
  trap_type = as.numeric(survey_df$trap_type_idx),
  
  CONSTANT = 10000
)

cat("✓ Constants updated\n")
cat(sprintf("  Model now uses %d grids instead of %d\n", 
            constants$n_grids, n_grids_total))


# ============================================================================
# UPDATED INITIAL VALUES
# ============================================================================

cat("\nCreating updated initial values function...\n")

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

# Impute missing values
z_temp_inits <- impute_missing_normally(z_temp_subset, verbose = TRUE)
z_rain_inits <- impute_missing_normally(z_rain_subset, verbose = TRUE)
z_RH_inits <- impute_missing_normally(data$z_RH, verbose = TRUE)
z_WS_rain_inits <- impute_missing_normally(data$z_WS_rain, verbose = TRUE)
z_land_inits <- impute_missing_normally(z_land_subset, verbose = TRUE)

# Initialize N for observed grids only
N_init <- matrix(NA_integer_, nrow = n_grids_model, ncol = n_dates)
N_init[,] <- pmax(1L, ceiling(mean(data$y, na.rm = TRUE) + 5L))

for(k in seq_len(nrow(survey_df))) {
  y_obs <- as.integer(survey_df$Culex[k])
  if(is.na(y_obs)) next
  
  s_local <- site_to_grid_updated[k]
  d <- constants$date_y[k]
  N_init[s_local, d] <- max(N_init[s_local, d], y_obs + 5L)
}

inits <- function() {
  
  log_lambda_init <- matrix(
    rnorm(n_grids_model * n_dates, mean = -2, sd = 0.5),
    nrow = n_grids_model,
    ncol = n_dates
  )
  lambda_grid_init <- exp(log_lambda_init)
  
  logit_p_init <- rnorm(n_grids_model, mean = -1, sd = 0.5)
  p_init <- plogis(logit_p_init)
  
  lambda_thinned_init <- matrix(NA_real_, nrow = n_grids_model, ncol = n_dates)
  for(g in 1:n_grids_model) {
    for(t in 1:n_dates) {
      lambda_thinned_init[g, t] <- lambda_grid_init[g, t] * 
        data$area_grid[g] * 
        p_init[g]
    }
  }
  
  background_init <- numeric(n_dates)
  for(t in 1:n_dates) {
    background_init[t] <- sum(lambda_thinned_init[, t]) / constants$n_obs_po
  }
  
  list(
    z_temp = z_temp_inits,
    z_rain = z_rain_inits,
    z_RH = z_RH_inits,
    z_WS_rain = z_WS_rain_inits,
    z_land = z_land_inits,
    
    beta0 = rnorm(1, -2, 0.3),
    beta_rain = rnorm(1, 0, 0.2),
    beta_temp = rnorm(1, 0, 0.2),
    # beta_buildings = rnorm(1, 0, 0.2),
    beta_land = rnorm(constants$n_land_covs, 0, 0.2),
    
    phi = runif(1, 1, 3),
    
    alpha0 = rnorm(1, qlogis(0.10), 0.2),
    alpha_RH = rnorm(1, 0, 0.1),
    alpha_WS_rain = rnorm(1, 0, 0.1),
    alpha_trap_raw = rnorm(4, 0, 0.2),
    sigma_trap = runif(1, 0.3, 0.8),
    
    delta0 = rnorm(1, -1, 0.3),
    delta_poi = rnorm(1, 0, 0.2),
    delta_reports = rnorm(1, 0 ,0.2),
    
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
cat("✓ Initial values created\n\n")


# ============================================================================
# UPDATED MODEL CODE
# ============================================================================

cat("Defining updated NIMBLE model...\n")

model_code <- nimbleCode({
  ## ============================================================
  ## PRIORS
  ## ============================================================
  beta0 ~ dnorm(0, sd = 1)
  beta_temp ~ dnorm(0, sd = 1)
  beta_rain ~ dnorm(0, sd = 1)
  
  for(k in 1:n_land_covs) {
    beta_land[k] ~ dnorm(0, sd = 1)
  }
  
  # Baseline detection from MRR literature
  alpha0 ~ dnorm(qlogis(0.10), sd = 0.5)  
  alpha_RH ~ dnorm(0, sd = 1)
  alpha_WS_rain ~ dnorm(0, sd = 1)
  
  for(i in 1:4) {
    alpha_trap_raw[i] ~ dnorm(0, sd = sigma_trap)
  }
  sigma_trap ~ dgamma(2, 2)
  
  ## Citizen science thinning
  delta0 ~ dnorm(-1, sd = 1)
  delta_poi ~ dnorm(0, sd = 1)
  delta_reports ~ dnorm(0, sd = 1)
  
  ## Dispersion parameters
  phi ~ dgamma(3, 1.5) 
  
  ## ============================================================
  ## SEPARABLE SPACE-TIME PARAMETERS
  ## ============================================================
  
  # Spatial variance (marginal variance of the spatial field)
  sigma_spatial ~ dgamma(2, 1)
  
  # Range parameter (controls spatial correlation decay)
  range_spatial ~ dgamma(5, 1)
  
  # Temporal AR(1) correlation parameter
  # rho_temporal in (-1, 1), but typically positive
  # Using logit transform: rho = 2*ilogit(rho_raw) - 1 maps to (-1,1)
  # Or simpler: rho ~ dunif(0, 1) for positive autocorrelation
  rho_temporal ~ dbeta(10, 2)  # Prior centered near high correlation
  
  ## Construct spatial covariance matrix (Matérn exponential)
  for(g1 in 1:n_grids_total) {
    for(g2 in 1:n_grids_total) {
      # Exponential correlation: rho(d) = exp(-d/range)
      Sigma_spatial[g1, g2] <- sigma_spatial^2 * exp(-dist_matrix[g1, g2] / range_spatial)
    }
  }
  
  ## Base spatial field (time-invariant component)
  spatial_effect_base[1:n_grids_total] ~ dmnorm(
    mean = rep(0, n_grids_total),
    cov = Sigma_spatial[1:n_grids_total, 1:n_grids_total]
  )
  
  ## Temporal innovations (AR(1) process)
  # First time point uses base field
  for(g in 1:n_grids_total) {
    spatial_effect[g, 1] <- spatial_effect_base[g]
  }
  
  # Subsequent time points follow AR(1)
  # spatial_effect[g,t] = rho * spatial_effect[g,t-1] + epsilon[g,t]
  # where epsilon ~ N(0, (1-rho^2)*sigma_spatial^2) to maintain variance
  sigma_temporal <- sigma_spatial * sqrt(1 - rho_temporal^2)
  
  for(t in 2:n_dates) {
    for(g in 1:n_grids_total) {
      epsilon[g, t] ~ dnorm(0, sd = sigma_temporal)
      spatial_effect[g, t] <- rho_temporal * spatial_effect[g, t-1] + epsilon[g, t]
    }
  }
  
  ## ============================================================
  ## DERIVED VARIABLES
  ## ============================================================
  for(i in 1:4) {
    alpha_trap[i] <- alpha_trap_raw[i]
  }
  alpha_trap[5] <- -sum(alpha_trap_raw[1:4])
  
  for(g in 1:n_grids_total) {
    logit_p[g] <- delta0 + delta_poi * z_poi[g] + 
      delta_reports * z_reports[g]
    p[g] <- ilogit(logit_p[g])
    
    for(t in 1:n_dates) {
      ## Include spatial-temporal random effect
      log_lambda[g, t] <- beta0 + 
        beta_temp * z_temp[g, t] +
        beta_rain * z_rain[g, t] +
        beta_buildings * z_buildings[g] +
        beta_ndvi * z_ndvi[g, t] +
        inprod(beta_land[1:n_land_covs], z_land[g, 1:n_land_covs]) +
        spatial_effect[g, t]
      
      lambda_grid[g, t] <- exp(log_lambda[g, t])
      lambda_thinned[g, t] <- lambda_grid[g, t] * area_grid[g] * p[g]
    }
  }
  
  for(t in 1:n_dates) {
    background[t] <- inprod(lambda_thinned[1:n_grids_total, t], 
                            rep(1, n_grids_total)) / n_obs_po
  }
  
  ## ============================================================
  ## LIKELIHOOD
  ## ============================================================
  
  for(i in 1:n_grids_obs) {
    for(t in 1:n_dates) {
      N[i, t] ~ dnegbin(
        prob = phi / (phi + (lambda_grid[obs_grid_idx[i], t])), 
        size = phi
      )
    }
  }
  
  for(k in 1:n_obs_y) {
    logit(p_trap[k]) <- alpha0 +
      alpha_RH * z_RH[k] +
      alpha_WS_rain * z_WS_rain[k] +
      alpha_trap[trap_type[k]]
    
    y[k] ~ dbin(p_trap[k], N[site_to_grid[k], date_y[k]])
  }
  
  for(r in 1:n_obs_po) {
    ones[r] ~ dbern(
      exp(
        log(lambda_thinned[po_grid_idx[r], date_po[r]]) - 
          log(background[date_po[r]])
      ) / CONSTANT
    )
  }
})

cat("✓ Model code defined\n\n")


# ============================================================================
# BUILD MODEL
# ============================================================================

cat("Building NIMBLE model (this may take a few minutes)...\n")

model <- nimbleModel(
  model_code,
  constants = constants,
  data = data,
  inits = initList,
  calculate = FALSE
)
model$initializeInfo()

model$calculate("alpha_trap")
model$calculate("p_trap")
model$calculate("lifted_phi_over__oPphi_plus_lambda_grid_oBg_comma_t_cB_cP_L29")
model$calculate("lifted_exp_oPlog_oPlambda_thinned_oBpo_grid_idx_oBr_cB_comma_date_po_oBr_cB_cB_cP_minus_log_oPbackground_oBdate_po_oBr_cB_cB_cP_cP_over_10000_L34")

model$initializeInfo()

cat("✓ Model built successfully!\n\n")

# cat("Memory saved:\n")
# cat(sprintf("  Before: %d grids × %d dates = %d cells\n", 
#             n_grids_total, n_dates, n_grids_total * n_dates))
# cat(sprintf("  After: %d grids × %d dates = %d cells\n", 
#             n_grids_model, n_dates, n_grids_model * n_dates))
# cat(sprintf("  Reduction: %.1fx smaller\n\n", 
#             (n_grids_total * n_dates) / (n_grids_model * n_dates)))

# Save updated setup
# save(
#   data, constants, inits, initList, model,
#   observed_grid_ids, grid_id_to_model_idx,
#   file = "model_setup_memory_efficient.RData"
# )

# cat("✓ Saved to: model_setup_memory_efficient.RData\n")

cModel <- compileNimble(model)

cat("Model compiled successfully.\n\n")

# ==============================================================================
# PARAMETERS TO MONITOR
# ==============================================================================

# monitors <- c("beta0", "beta_temp", "beta_rain", "beta_buildings", "beta_land", "beta_ndvi",
#               "alpha0", "alpha_RH", "alpha_WS_rain", "sigma_trap", "alpha_trap",
#               "p_trap",
#               "delta0", "delta_poi", "delta_reports",
#               "phi",
#               "background", "lambda_thinned")

# --- Core parameters of interest ---
params_core <- c(
  # Intensity model (main ecological effects)
  "beta0", "beta_temp", "beta_rain", #"beta_ndvi", "beta_buildings",
  "beta_land",  # All land use coefficients
  
  # Detection model (trap survey)
  "alpha0", "alpha_RH", "alpha_WS_rain",
  "alpha_trap",  # All 5 trap effects (including derived one)
  "sigma_trap",  # Variability among trap types
  "p_trap",
  
  # Thinning model (citizen science)
  "delta0", "delta_poi", "delta_reports", "p",
  
  # Dispersion
  "phi"
)

# --- Derived quantities for interpretation ---
params_derived <- c(
  # Spatial thinning probabilities (detection prob for CS)
  "p",  # Vector of length n_grids_total
  
  # Mean detection probabilities by trap type (for summary)
  # Calculate post-hoc from alpha_trap samples
  
  # Background (for model checking)
  "background"  # Check if stable across time
)

# --- Abundance estimates (if you want predictions) ---
params_abundance <- c(
  "N",  # True abundance at surveyed grids
  "lambda_grid"  # Expected intensity (can get predictions anywhere)
)

# --- For model diagnostics (optional, increases file size) ---
params_diagnostic <- c(
  "log_lambda",      # Linear predictor
  "lambda_thinned",  # For checking PO likelihood
  "logit_p",         # Linear predictor for thinning
  "p_trap"           # Individual detection probabilities (large!)
)

monitors_base <- c(params_core, "background", "lambda_thinned")
monitors_minimal <- c(params_core, "background") ## Minimal (for initial runs / convergence checks)
monitors_standard <- c(params_core, "p", "background") ## Standard (for inference)
monitors_full <- c(params_core, params_derived, params_abundance) ## Full (for detailed analysis, large file)
monitors_diagnostic <- c(params_core, params_diagnostic) ## Diagnostic (for troubleshooting convergence)

# ==============================================================================
# SETUP MCMC
# ==============================================================================

# choose monitors
monitors <- monitors_diagnostic

cat("Note: parameters to monitors were:", 
    monitors,".\n\n")

mcmc_conf <- configureMCMC(model,
                           monitors = monitors,
                           enableWAIC = TRUE)
mcmc <- buildMCMC(mcmc_conf)

cat("Compiling MCMC (this will take a LONG time)...\n")
cmcmc <- compileNimble(mcmc, project = cModel)

cat("MCMC compiled successfully.\n\n")

# ==============================================================================
# 12. RUN MCMC CHAINS
# ==============================================================================
cat("Running MCMC chains...\n")
# cat("Configuration: 3 chains, 10000 iterations, 5000 burn-in, thin = 25\n\n")

nchains <- 3
seeds <- c(1, 2, 3)
niterations = 30000
nburnins = 10000
nthin = 20

cat("Configuration:", 
    nchains,"chains, ", niterations, "iterations, ", 
    nburnins, "burnin, ", nthin, "thinning rate \n\n")

samples_list <- pblapply(1:nchains, function(i) {
  result <- runMCMC(
    cmcmc,
    niter = niterations,
    nburnin = nburnins,
    thin = nthin,
    nchains = 1,
    samplesAsCodaMCMC = TRUE,
    WAIC = TRUE,
    setSeed = seeds[i],
    progressBar = TRUE
  )
}, cl = nchains)

cat("\nAll chains completed successfully!\n\n")
saveRDS(samples_list, "/users/2601581k/sharedscratch/onesTrick_pres.RDS")


# ===============================================================================
# POSTERIOR EVALUATION

library(coda)
samples <- lapply(samples_list, function(x) x$samples)

mcmc_list <- mcmc.list(lapply(samples_list, function(x) {
  if(!is.null(x$samples)) return(x$samples)
  # if runMCMC returned the mcmc object directly, adapt:
  if(inherits(x, "mcmc")) return(x)
  stop("Can't find coda samples in an element of samples_weak")
}))


gelman.diag(mcmc_list)

gelman.diag(mcmc_list, multivariate = FALSE)

params <- c(
  "beta0", "beta_temp", "beta_rain", "beta_land[1]", "beta_land[2]", "beta_land[3]", "beta_land[4]", 
  "beta_land[5]", "beta_land[6]", "beta_land[7]", "beta_land[8]", "beta_land[9]",
  "alpha0", "alpha_RH", "alpha_WS_rain","alpha_trap[1]",
  "alpha_trap[2]","alpha_trap[3]","alpha_trap[4]", "alpha_trap[5]",
  "delta0", "delta_reports", "delta_poi", "phi")

mcmc_subset <- mcmc.list(
  lapply(mcmc_list , function(chain) {
    m <- as.matrix(chain)
    # subset the parameters
    m_sub <- m[, params, drop = FALSE]
    mcmc(m_sub)
  })
)

gelman.diag(mcmc_subset, multivariate = FALSE)

effectiveSize(mcmc_subset)
traceplot(mcmc_subset)
autocorr.plot(mcmc_subset)

#===================================================================
# POSTERIOR AREAS AND PLOTTING
#===================================================================
library(bayesplot)

samples <- mcmc.list(lapply(samples_list, function(x) {
  mcmc(x$samples)
}))

params <- c("beta_temp", "beta_rain", 
            "beta_land[1]", "beta_land[2]", 
            "beta_land[3]", "beta_land[4]",
            "beta_land[5]", "beta_land[6]", 
            "beta_land[7]", "beta_land[8]", "beta_land[9]")

param_names <- c(
  "beta_temp" = "3 week prior temperature",
  "beta_rain" = "4 week prior cumulative rainfall",
  "beta_land[1]" = "Wetland",
  "beta_land[2]" = "Freshwater", 
  "beta_land[3]" = "Woodland", 
  "beta_land[4]" = "Grassland/Heather", 
  "beta_land[5]" = "Arable",
  "beta_land[6]" = "Livestock",
  "beta_land[7]" = "Elevation",
  "beta_buildings" = "Building density",
  "beta_ndvi" = "NDVI"
)
color_scheme_set("darkgray")

posterior_intervals <- mcmc_intervals(samples[, params],
                                      point_est = "median") +
  theme_minimal() + 
  scale_y_discrete(labels = param_names)

posterior_intervals <- posterior_intervals +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    axis.text.y = element_text(size = 14)  # 
  )



posterior_intervals


