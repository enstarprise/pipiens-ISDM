################################################################################
# Cx pipiens in Scotland ISDM
# Author: shennice
# Date: December 2025

# This is the full integrated species distribution model for Cx. pipiens.
# The survey data is modelled through an N-mixture model for the count data,
# the citizen science reports is modelled through a Bernoulli poisson point
# process using the ones trick to account for the background reporting.
################################################################################

# ==============================================================================
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


# ============================================================================
# LOAD COVARIATES (z-scored)
# ============================================================================
cat("Loading data...\n")
loaded_data <- load("01_data/loaded_data.RData")
z_land_data <- load("01_data/env_land_data.RData")
z_temp <- load("01_data/env_temp_data.RData")
z_rain_part1 <- load("01_data/env_rain_data_part1.RData")
z_rain_part2 <- load("01_data/env_rain_data_part2.RData")

load("01_data/loaded_data.RData")
load("01_data/env_land_data.RData")
load("01_data/env_temp_data.RData")
load("01_data/env_rain_data_part1.RData")
load("01_data/env_rain_data_part2.RData")
z_rain <- rbind(z_rain_part1, z_rain_part2)


remove(z_rain_part1)
remove(z_rain_part2)
gc()


# organise land covariates
z_land <- as.data.frame(z_land_data) %>%
  dplyr::arrange(grid_id_area2km) %>%
  dplyr::select(starts_with("z_")) %>%
  dplyr::select(-c(z_urban_km2, z_suburban_km2, z_poi_count_grouped, z_reports)) %>%
  as.matrix()

z_poi <- as.data.frame(z_land_data) %>%
  arrange(grid_id_area2km) %>%
  pull(z_poi_count_grouped)

z_reports <- as.data.frame(z_land_data) %>%
  arrange(grid_id_area2km) %>%
  pull(z_reports)


# ============================================================================
# CONVERT CLIMATE DATA TO MATRICES (grids × dates)
# ============================================================================

cat("\nConverting climate data to matrix format...\n")

# Get sorted unique dates for column ordering
# climate_dates_sorted <- sort(unique(z_climate_data$date))

# Create temperature matrix (rows = grids, columns = dates)
# z_temp <- z_climate_data %>%
#   dplyr::select(grid_id_area2km, date, z_temp_min) %>%
#   pivot_wider(
#     names_from = date,
#     values_from = z_temp_min
#   ) %>%
#   arrange(grid_id_area2km) %>%
#   dplyr::select(-grid_id_area2km) %>%
#   as.matrix()

# Create rainfall matrix (rows = grids, columns = dates)
# z_rain <- z_climate_data %>%
#   dplyr::select(grid_id_area2km, date, z_rain) %>%
#   pivot_wider(
#     names_from = date,
#     values_from = z_rain
#   ) %>%
#   arrange(grid_id_area2km) %>%
#   dplyr::select(-grid_id_area2km) %>%
#   as.matrix()

cat(sprintf("  Temperature matrix: %d grids × %d dates\n", nrow(z_temp), ncol(z_temp)))
cat(sprintf("  Rainfall matrix: %d grids × %d dates\n", nrow(z_rain), ncol(z_rain)))
cat(sprintf("  Temperature completeness: %.1f%%\n", 100 * mean(!is.na(z_temp))))
cat(sprintf("  Rainfall completeness: %.1f%%\n", 100 * mean(!is.na(z_rain))))


# Choose numerical constant for stability
# Larger values (1000-10000) help when probabilities are very small
CONSTANT <- 10000

cat("Presence-only records:", nrow(cs_presences), "\n")
cat("Using Bernoulli ones trick with CONSTANT =", CONSTANT, "\n\n")


# ==============================================================================
# 7. MODEL DATA AND CONSTANTS
# ==============================================================================
cat("Preparing model data and constants...\n")

data <- list(
  z_rain = z_rain,
  z_temp = z_temp,
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
cat("Setting up initial values...\n")

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
z_temp_inits <- impute_missing_normally(z_temp, verbose = TRUE)
z_rain_inits <- impute_missing_normally(z_rain, verbose = TRUE)
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
    z_temp = z_temp_inits,
    z_rain = z_rain_inits,
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

cat("Initial values set.\n\n")


# ==============================================================================
# 9. NIMBLE MODEL CODE
# ==============================================================================
cat("Defining NIMBLE model...\n")

model_code <- nimbleCode({
  ## ============================================================
  ## PRIORS
  ## ============================================================
  beta0 ~ dnorm(0, sd = 1)
  beta_temp ~ dnorm(0, sd = 1)
  beta_rain ~ dnorm(0, sd = 1)
  # beta_ndvi ~ dnorm(0, sd = 1)
  
  for(k in 1:n_land_covs) {
    beta_land[k] ~ dnorm(0, sd = 1)
  }
  
  # Baseline detection from MRR literature (KEEP THIS - it's informed!)
  alpha0 ~ dnorm(qlogis(0.10), sd = 0.5)  
  
  # Microclimate effects on detection
  alpha_RH ~ dnorm(0, sd = 1)
  alpha_WS_rain ~ dnorm(0, sd = 1)
  
  # Trap-specific deviations from baseline (sum to zero)
  for(i in 1:4) {
    alpha_trap_raw[i] ~ dnorm(0, sd = sigma_trap)
  }
  sigma_trap ~ dgamma(2, 2)  # Regularizes how much traps can differ
  
  
  ## Citizen science thinning
  delta0 ~ dnorm(-1, sd = 1)
  delta_poi ~ dnorm(0, sd = 1)
  delta_reports ~ dnorm(0, sd = 1)
  
  ## Dispersion parameters
  phi ~ dgamma(3, 1.5) 
  
  ## ============================================================
  ## DERIVED VARIABLES
  ## ============================================================
  # Constrain trap effects to sum to zero
  for(i in 1:4) {
    alpha_trap[i] <- alpha_trap_raw[i]
  }
  alpha_trap[5] <- -sum(alpha_trap_raw[1:4])
  
  for(g in 1:n_grids_total) {
    ## Thinning process (detection probability for CS data)
    logit_p[g] <- delta0 + delta_poi * z_poi[g] + 
      delta_reports * z_reports[g]
    
    p[g] <- ilogit(logit_p[g])
    
    ## Intensity for each grid and time
    for(t in 1:n_dates) {
      log_lambda[g, t] <- beta0 + 
        beta_temp * z_temp[g, t] +
        beta_rain * z_rain[g, t] +
        beta_buildings * z_buildings[g] +
        beta_ndvi * z_ndvi[g, t] +
        inprod(beta_land[1:n_land_covs], z_land[g, 1:n_land_covs])  
      
      lambda_grid[g, t] <- exp(log_lambda[g, t])
      
      ## Thinned intensity (for PO likelihood)
      lambda_thinned[g, t] <- lambda_grid[g, t] * area_grid[g] * p[g]
      
    }
  }
  
  ## Background interval: sum of thinned intensities across ALL cells
  ## Divide by number of PO points (Riemann sum approximation)
  ## This is the normalizing denominator
  for(t in 1:n_dates) {
    background[t] <- inprod(lambda_thinned[1:n_grids_total, t], 
                            rep(1, n_grids_total)) / n_obs_po
  }
  
  
  ## ============================================================
  ## LIKELIHOOD
  ## ============================================================
  
  ## True abundance at observed grids
  for(i in 1:n_grids_obs) {
    for(t in 1:n_dates) {
      N[i, t] ~ dnegbin(
        prob = phi / (phi + (lambda_grid[obs_grid_idx[i], t])), 
        size = phi
      )
    }
  }
  
  ## Survey detection
  for(k in 1:n_obs_y) {
    logit(p_trap[k]) <- alpha0 +                      # informed baseline from MRR
      alpha_RH * z_RH[k] +          # how humidity affects detection
      alpha_WS_rain * z_WS_rain[k] + # how weather affects detection
      alpha_trap[trap_type[k]]      # how trap type differs from average
    
    # observation likelihood
    y[k] ~ dbin(p_trap[k], N[site_to_grid[k], date_y[k]])
  }
  
  ## Presence-only likelihood using Bernoulli ones trick
  ## Following Fidino (2021) / Koshkina et al. (2017)
  ## Each presence contributes: lambda_thinned / po_denominator
  for(r in 1:n_obs_po) {
    ones[r] ~ dbern(
      exp(
        log(lambda_thinned[po_grid_idx[r], date_po[r]]) - 
          log(background[date_po[r]])
      ) / CONSTANT
    )
  }
})
cat("NIMBLE model defined.\n\n")

# ==============================================================================
# 10. MODEL BUILDING AND COMPILATION
# ==============================================================================
cat("Building NIMBLE model...\n")

model <- nimbleModel(model_code,
                     constants = constants,
                     data = data,
                     inits = initList,
                     calculate = FALSE)

model$initializeInfo()

model$initializeInfo()

# Calculate key nodes
cat("Calculating initial nodes...\n")
# model$calculate("p")
model$calculate("alpha_trap")
model$calculate("p_trap")
#model$calculate("lifted_phi_over__oPphi_plus__oPlambda_grid_oBobs_grid_idx_oBi_cB_comma_t_cB_cP_cP_L32")
#model$calculate("lifted_exp_oPlog_oPlambda_thinned_oBpo_grid_idx_oBr_cB_comma_date_po_oBr_cB_cB_cP_minus_log_oPbackground_oBdate_po_oBr_cB_cB_cP_cP_over_10000_L37")
model$initializeInfo()

cat("Compiling model (this may take several minutes)...\n")
cModel <- compileNimble(model)
cat("Model compiled successfully.\n\n")

cat("Configuring MCMC...\n")

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
  "beta0", "beta_temp", "beta_rain", 
  "beta_land",  # All land use coefficients: NOT INCLUDING NDVI
  
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
monitors <- monitors_base

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
niterations = 15000
nburnins = 5000
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
saveRDS(samples_list, "/users/2601581k/sharedscratch/model.RDS")

