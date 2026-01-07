################################################################################
# Cx pipiens in Scotland ISDM
# Author: shennice
# Date: December 2025

# This is the full integrated species distribution model for Cx. pipiens.
# The survey data is modelled through an N-mixture model for the count data,
# the citizen science reports is modelled through a Bernoulli poisson point
# process using the ones trick to account for the background reporting.

# This is using the marginlisation of the latent discrete parameter N.
# Must sum over all possible values, but since it ranges from 0 to infinity, 
# there must be a cap/max...
################################################################################
# Clear workspace
# rm(list = ls())

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

#===============================================================================
# PACKAGE MANAGEMENT
#===============================================================================

# Define required packages
required_packages <- c(
  # Modeling
  "cmdstanr"       #cmdstanr

)

# Suppress startup messages for cleaner output
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(tidyr)
  library(cmdstanr)
})


# ============================================================================
# \PREPARE DATA 
# ============================================================================

#  data preparation script
# creates: stan_data_marginalized.rds
# source("02_model/stan/01_prepare_stan_data_impute.R")

# # load 
# load("01_data/processedCovariates/1km/stan_data_init_obsgrid_imputation_1km.RData")
# load("01_data/processedCovariates/1km/stan_data_init_obsgrid_1km.RData")
# load("01_data/processedCovariates/1km/01_5kpixel_test/5k_pixelSubset.RData"); stan_data <- subset_stan_data

# stan_data <- readRDS("01_data/processedCovariates/1km/stan_data_imputation_1km.rds")
# obs_grid_ids <- readRDS("01_data/processedCovariates/1km/observed_grid_ids_1km.rds")

# load diagnostic function
# source("02_model/stan/02_stan_test_utilities.R")
# diagnose_stan_data(stan_data)

# ============================================================================
# TEST oN SIM DATA 
# ===========================================================================

#test the function alone
# test_marginalization_function()
# 
# # quick test
# test_result <- quick_model_test(
#   stan_file = "02_model/marginalised_model.stan",
#   n_grids = 100,
#   n_dates = 100,
#   n_obs_y = 30,
#   chains = 2,
#   iter_warmup = 100,
#   iter_sampling = 100
# )


# ============================================================================
#  MODEL RUN 
# ============================================================================

stan_data$N_multiplier <- 5
stan_data$grainsize_survey <- 10
stan_data$grainsize_po <- 1
# stan_data$grainsize <- 10
# stan_data$N_max_quantile <- 0.999 # 0.999 for 99.9th percentile

model <- cmdstan_model("02_model/stan/marginalised_model_parallel_v2.stan",
                       cpp_options = list(stan_threads = TRUE))

# fit
fit <- model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 4,  # Use 4 threads per chain
  thin = 10,
  iter_warmup = 100,
  iter_sampling = 250,
  refresh = 25,
  adapt_delta = 0.95,  # inc if  divergences ??
  max_treedepth = 12,
  init = init_fun, # or a standard 0.5, but chain 3 gets stuck
  seed = 123
)

## with reduce sum, takes 11 minutes

# ============================================================================
# DIAGNOSTICS
# ============================================================================

# diag summary
diag <- fit$diagnostic_summary()
print(diag)

#  warnings
if (any(diag$num_divergent > 0)) {
  cat("\n DIVERGENCES DETECTED\n")
  cat("  try adapt_delta 0.99\n")
}

if (any(diag$num_max_treedepth > 0)) {
  cat("\n MAX TREEDEPTH EXCEEDED >0\n")
  cat("  try inc max_treedepth = 15\n")
}

# convergence (Rhat)
params_summary <- fit$summary()
high_rhat <- params_summary %>% filter(rhat > 1.05)

if (nrow(high_rhat) > 0) {
  print(high_rhat[, c("variable", "rhat", "ess_bulk")])
  cat("\n  try running longer chains??? \n")
} else {
  cat("\n✓ all rhat values < 1.05\n")
}

# ============================================================================
# EXTRACT AND SUMMARIZE RESULTS
# ============================================================================

# parameters of interest
main_params <- c(
  "beta0", "beta_temp", "beta_rain", "beta_land",
  "alpha0", "alpha_RH", "alpha_WS_rain", "sigma_trap", "alpha_trap",
  "delta0", "delta_poi", "delta_reports",
  "phi"
)

results <- fit$summary(main_params)
print(results)

# ============================================================================
# SAVE RESULTS
# ============================================================================
# save fit object
fit$save_object("03_posterior/results/stan_fit_marginalized_1km_5kpixels.rds")

# Save summary
write_csv(results, "03_posterior/results/parameter_estimates_1km.csv")

# extract draws 
draws <- fit$draws(format = "df")
saveRDS(draws, "03_posterior/results/posterior_draws_1km.rds")

fit <- readRDS("03_posterior/results/stan_fit_marginalized_1km.rds")
draws <- readRDS("03_posterior/results/posterior_draws_1km.rds")

# ============================================================================
# PPC
# ============================================================================

# extract N_expected
N_expected <- fit$summary("N_expected")

cat("Expected abundance statistics:\n")
cat(sprintf("  Mean: %.2f\n", mean(N_expected$mean)))
cat(sprintf("  SD: %.2f\n", sd(N_expected$mean)))
cat(sprintf("  Range: [%.2f, %.2f]\n", 
            min(N_expected$mean), max(N_expected$mean)))

# ============================================================================
# POSTEIROR VISUALISATION
# ============================================================================
library(bayesplot)
main_params <- c(
  "beta0", "beta_temp", "beta_rain", "beta_land[1]", "beta_land[2]", "beta_land[3]",
  "beta_land[4]", "beta_land[5]", "beta_land[6]", "beta_land[7]", "beta_land[8]",
  "alpha0", "alpha_RH", "alpha_WS_rain", "sigma_trap", "alpha_trap[1]", "alpha_trap[1]",
  "alpha_trap[2]", "alpha_trap[3]", "alpha_trap[4]", "alpha_trap[5]",
  "delta0", "delta_poi", "delta_reports",
  "phi"
)

lambda_params<- c(
  "beta0", "beta_temp", "beta_rain", "beta_land[1]", "beta_land[2]", "beta_land[3]",
  "beta_land[4]", "beta_land[5]", "beta_land[6]", "beta_land[7]", "beta_land[8]"
)

posterior_mat <- fit$draws(lambda_params) %>% posterior::as_draws_matrix()

# posterior density areas
bayesplot::mcmc_areas(
  posterior_mat,
  pars = lambda_params,
  prob = 0.8,        # 80% interval highlight
  prob_outer = 0.95  # full 95% interval
) + ggplot2::ggtitle("Posterior Distributions with Credible Intervals") +
  theme_bw()


