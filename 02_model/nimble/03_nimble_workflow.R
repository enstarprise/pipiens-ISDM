################################################################################
# Cx pipiens in Scotland ISDM
# Author: shennice
# Date: December 2025

# This is the full integrated species distribution model for Cx. pipiens.
# The survey data is modelled through an N-mixture model for the count data,
# the citizen science reports is modelled through a Bernoulli poisson point
# process using the ones trick to account for the background reporting.

################################################################################
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
  # Modeling
  "cmdstanr"       #cmdstanr
  
)

# Suppress startup messages for cleaner output
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(tidyr)
  library(nimble)
  library(coda)
  library(bayesplot)
})


# ============================================================================
# nimble model
# ============================================================================
source("02_model/nimble/01_prepare_nimble_model_obs.R")
source("02_model/nimble/02_nimble_model_exp.R")

model <- nimbleModel(
  code = model_code,
  constants = constants,
  data = data,
  inits = initList,
  calculate = FALSE
)
model$initializeInfo()

cmodel <- compileNimble(model)

conf <- configureMCMC(cmodel)
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = cmodel)

library(pbapply)
nchains <- 3
seeds <- c(1, 2, 3)
samples_list <- pblapply(1:nchains, function(i) {
  runMCMC(
    cmcmc,
    niter = 60000,
    nburnin = 30000,
    thin = 20,
    nchains = 1,
    samplesAsCodaMCMC = TRUE,
    setSeed = seeds[i],
    progressBar = FALSE 
  )
}, cl = nchains)   




# ============================================================================
# POSTERIOR
# ============================================================================

library(coda)

samples <- mcmc.list(
  lapply(samples_list, function(chain_matrix) {
    mcmc(chain_matrix)   # convert each chain's matrix to an mcmc object
  })
)


params <-  c("beta_land[1]", "beta_land[2]", 
             "beta_land[3]", "beta_land[4]",
             "beta_land[5]", "beta_land[6]", 
             "beta_land[7]", "beta_land[8]", 
             "beta0", "beta_rain", "beta_temp",
             "alpha0", "alpha_RH","alpha_WS_rain",
             "delta0","delta_poi", "delta_reports",
             "phi")

rhat <- gelman.diag(samples[, params])
rhat
rhat_trap <-  gelman.diag(samples[, c("alpha_trap[1]", "alpha_trap[2]", 
                                      "alpha_trap[3]", "alpha_trap[4]",
                                      "alpha_trap[5]")])

effectiveSize(samples[, params])
autocorr.plot(samples[,params])

summary(samples[,"phi"])

mcmc_areas(samples[, params],
           point_est = "median")

params_climate <-  c("beta_temp", "beta_rain")
mcmc_areas(samples[, params_climate],
           point_est = "median") +
  theme_minimal()

params_landuse <-  c("beta_land[1]", "beta_land[2]", 
                     "beta_land[3]", "beta_land[4]",
                     "beta_land[5]", "beta_land[6]", 
                     "beta_land[7]", "beta_land[8]", "beta_land[9]")
mcmc_areas(samples[, params_landuse],
           point_est = "median") +
  theme_minimal()


params_trap <- c("alpha_trap[1]", "alpha_trap[2]", 
                 "alpha_trap[3]", "alpha_trap[4]",
                 "alpha_trap[5]", "alpha0", "alpha_RH","alpha_WS_rain")
params_microclimate <- c("alpha0", "alpha_RH","alpha_WS_rain")
mcmc_areas(samples[, params_microclimate],
           point_est = "median") +
  theme_minimal()

params_cs <- c("delta0", "delta_pop",
               "delta_poi")
mcmc_areas(samples[, params_cs],
           point_est = "median") +
  theme_minimal()

params <- c("beta_temp", "beta_rain", 
            "beta_land[1]", "beta_land[2]", 
            "beta_land[3]", "beta_land[4]",
            "beta_land[5]", "beta_land[6]", 
            "beta_land[7]", "beta_buildings")

param_names <- c(
  "beta_temp" = "1 week temperature",
  "beta_rain" = "4 week cumulative rainfall",
  "beta_land[1]" = "Wetland",
  "beta_land[2]" = "Freshwater", 
  "beta_land[3]" = "Woodland", 
  "beta_land[4]" = "Grassland/Heather", 
  "beta_land[5]" = "Arable",
  "beta_land[6]" = "Livestock",
  "beta_land[7]" = "Elevation",
  "beta_buildings" = "Building density"
)
color_scheme_set("blue")
mcmc_areas_plot <- mcmc_areas(samples[, params],
                              point_est = "median") +
  theme_minimal() + 
  scale_y_discrete(labels = param_names)

ggsave(filename = "priorPosteriorPlots/mcmc_areas_lambdaFullGridxTime2.png",
       plot = mcmc_areas_plot,
       device = "png",
       width = 7,
       height = 8,
       units = "in",
       dpi = 300)


subset_samples <- lapply(samples, function(x) {
  x[, params_landuse, drop = FALSE]
})
subset_samples <- as.mcmc.list(subset_samples)
combined_samples <- do.call(rbind, subset_samples)

if (is.matrix(combined_samples)) {
  mcmc_obj <- as.mcmc(combined_samples)
} else {
  mcmc_obj <- combined_samples  # use as-is if already mcmc
}

color_scheme_set("blue")
mcmc_intervals_plot <- mcmc_intervals(mcmc_obj, pars = colnames(mcmc_obj), 
                                      prob = 0.95, prob_outer = 1) +
  scale_y_discrete(labels = param_names) +
  theme_minimal()


ggsave(filename = "priorPosteriorPlots/mcmc_intervals_lambadaFullGridxTime.png",
       plot = mcmc_intervals_plot2,
       device = "png",
       width = 7,
       height = 8,
       units = "in",
       dpi = 300)

