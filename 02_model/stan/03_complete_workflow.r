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
gc()

suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(tidyr)
  library(cmdstanr)
})

# stan_data_GIRDRES_MODELTYPE_GRIDSCOPE_QUADRATIC.RData
load("01_data/processedCovariates/2.5km/stan_data_2.5km_NNGP_all_both.RData")

# ============================================================================
#  MODEL RUN 
# ============================================================================

stan_data$N_max <- 15


test_init <- init_fun()
print("Initial beta_land values:")
print(test_init$beta_land)
print("Should be around ±0.03")


model <- cmdstan_model("02_model/stan/model_nngp.stan",
                       cpp_options = list(stan_threads = TRUE))


fit <- model$sample(
  data = stan_data,
  chains = 3,
  parallel_chains = 3,
  threads_per_chain = 4,
  iter_warmup = 500,
  iter_sampling = 500,
  adapt_delta = 0.95,
  max_treedepth = 12,
  init = init_fun,
  seed = 123,
  refresh = 25,
  output_dir = "02_model/stan/modelRuns/",
  show_messages = TRUE  
)


saveRDS(fit,"02_model/stan/modelRuns/model_nonspatial_poisson.rds")


fit_nb <- readRDS("02_model/stan/modelRuns/model_nonspatial2.rds")
# fit_quad <- readRDS("02_model/stan/modelRuns/model_nonspatial_quad_both.rds")
# fit_hsgp <- readRDS("02_model/stan/modelRuns/model_spatialHSGP.rds")

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

# ============================================================================
# EXTRACT AND SUMMARIZE RESULTS
# ============================================================================

# parameters of interest
main_params <- c(
  "beta0", "beta_temp", "beta_rain", "beta_land",
  "alpha0", "alpha_RH", "alpha_WS_rain", "sigma_trap", "alpha_trap",
  "delta0", "delta_poi", "delta_reports"
)

results <- fit$summary(main_params)
print(results)

# # Save summary
# write_csv(results, "03_posterior/results/parameter_estimates.csv")


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

# ===========================================================================
# MODEL COMPARISON 
# ===========================================================================
library(loo)
library(posterior)

# survey log-lik
log_lik_survey_1 <- fit$draws("log_lik_survey")
log_lik_survey_2 <- fit_nb$draws("log_lik_survey")

# convert to matrix: draws x observations
ll_survey_1 <- as_draws_matrix(log_lik_survey_1)
ll_survey_2 <- as_draws_matrix(log_lik_survey_2)

loo_survey_1 <- loo(ll_survey_1)
loo_survey_2 <- loo(ll_survey_2)

loo_compare(loo_survey_1, loo_survey_2)

plot(loo_survey_1)
plot(loo_survey_2)

loo_survey_1
loo_survey_2

# citizen science
ll_po_1 <- as_draws_matrix(fit$draws("log_lik_po"))
ll_po_2 <- as_draws_matrix(fit_nb$draws("log_lik_po"))

loo_po_1 <- loo(ll_po_1)
loo_po_2 <- loo(ll_po_2)

loo_compare(loo_po_1, loo_po_2)

loo_po_1
loo_po_2


ll_joint_1 <- cbind(ll_survey_1, ll_po_1)
ll_joint_2 <- cbind(ll_survey_2, ll_po_2)

loo_joint_1 <- loo(ll_joint_1)
loo_joint_2 <- loo(ll_joint_2)

loo_compare(loo_joint_1, loo_joint_2)


# ============================================================================
# PPC for SURVEY
# ============================================================================
# Extract posterior predictive samples
y_rep_draws <- fit_nb$draws("y_rep", format = "draws_matrix")

# summary statistics for each  draw
ppc_summary <- function(y_rep_matrix, y_obs) {
  n_draws <- nrow(y_rep_matrix)
  n_obs <- length(y_obs)
  
  #  statistics 
  stats <- matrix(NA, nrow = n_draws, ncol = 4)
  colnames(stats) <- c("mean", "sd", "max", "zeros")
  
  for (i in 1:n_draws) {
    y_rep_i <- y_rep_matrix[i, ]
    stats[i, "mean"] <- mean(y_rep_i)
    stats[i, "sd"] <- sd(y_rep_i)
    stats[i, "max"] <- max(y_rep_i)
    stats[i, "zeros"] <- sum(y_rep_i == 0) / n_obs
  }
  
  # Observed statistics
  obs_stats <- c(
    mean = mean(y_obs),
    sd = sd(y_obs),
    max = max(y_obs),
    zeros = sum(y_obs == 0) / n_obs
  )
  
  return(list(
    posterior_stats = stats,
    observed_stats = obs_stats,
    bayes_p_values = colMeans(t(apply(stats, 1, function(x) x >= obs_stats)))
  ))
}

# Run PPC
ppc_results <- ppc_summary(y_rep_draws, stan_data$y)

# Print results
cat("Posterior Predictive Check Results:\n")
cat("----------------------------------\n")
for (stat in names(ppc_results$observed_stats)) {
  cat(sprintf("%-10s: Observed = %.3f, Posterior mean = %.3f, p-value = %.3f\n",
              stat,
              ppc_results$observed_stats[stat],
              mean(ppc_results$posterior_stats[, stat]),
              ppc_results$bayes_p_values[stat]))
}

# Visual PPC
library(bayesplot)
library(ggplot2)

# Plot 1: Distribution of test statistics
color_scheme_set("brightblue")
ppc_stat(y = stan_data$y, yrep = y_rep_draws[1:500, ], stat = "mean") +
  ggtitle("Posterior Predictive Check: Mean") +
  theme_minimal()

ppc_stat(y = stan_data$y, yrep = y_rep_draws[1:500, ], stat = "sd") +
  ggtitle("Posterior Predictive Check: Standard Deviation")

ppc_stat(y = stan_data$y, yrep = y_rep_draws[1:500, ], stat = function(x) max(x)) +
  ggtitle("Posterior Predictive Check: Maximum")

ppc_stat(y = stan_data$y, yrep = y_rep_draws[1:500, ], stat = function(x) mean(x == 0)) +
  ggtitle("Posterior Predictive Check: Proportion of Zeros")

# Plot 2: Distribution of observations vs predictions
ppc_dens_overlay(y = stan_data$y, yrep = y_rep_draws[1:50, ]) +
  ggtitle("Distribution: Observed vs Posterior Predictive")

# Plot 3: Rootogram (good for count data)
ppc_rootogram(y = stan_data$y, yrep = y_rep_draws[1:500, ]) +
  ggtitle("Rootogram: Observed vs Predicted Counts")

# Plot 4: ECDF comparison
ppc_ecdf_overlay(y = stan_data$y, yrep = y_rep_draws[1:50, ]) +
  ggtitle("ECDF: Observed vs Posterior Predictive")

# ============================================================================
# PPC for N MAX MARGINALISATION
# ============================================================================


# Extract N_max diagnostics
n_near_boundary <- fit$draws("n_near_boundary", format = "draws_matrix")
max_N_ratio <- fit$draws("max_N_ratio", format = "draws_matrix")

# Summary of N_max diagnostics
cat("\nN_max Adequacy Diagnostics:\n")
cat("---------------------------\n")
cat(sprintf("Maximum N/N_max ratio across all draws: %.4f\n", max(max_N_ratio)))
cat(sprintf("Mean N/N_max ratio: %.4f\n", mean(max_N_ratio)))
cat(sprintf("Proportion of draws with N/N_max > 0.9: %.4f\n", 
            mean(n_near_boundary > 0)))

# Check if N_max needs adjustment
if (mean(max_N_ratio) > 0.8) {
  cat("\n⚠️ WARNING: N_max may be too small!\n")
  cat(sprintf("  Mean N/N_max ratio is %.2f (should be < 0.8)\n", mean(max_N_ratio)))
  cat(sprintf("  %.1f%% of draws have simulations near the boundary\n", 
              mean(n_near_boundary > 0) * 100))
  
  # Calculate recommended N_max
  current_N_max <- stan_data$N_max  # You need to track this
  max_simulated_N <- max(fit$draws("N_sim", format = "draws_matrix"))  # Need to store N_sim
  
  cat(sprintf("\nRecommended adjustments:\n"))
  cat(sprintf("  Current N_max: %d\n", current_N_max))
  cat(sprintf("  Maximum simulated N: %d\n", max_simulated_N))
  cat(sprintf("  Recommended N_max: %d (2× max simulated)\n", max_simulated_N * 2))
} else {
  cat("\n✓ N_max appears adequate.\n")
}

# Plot distribution of N/N_max ratios
ggplot(data.frame(ratio = as.numeric(max_N_ratio)), aes(x = ratio)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 0.9, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Distribution of N/N_max Ratios",
       x = "N / N_max ratio",
       y = "Frequency") +
  theme_minimal()


# ============================================================================
# SPATIAL FIELD VISUALISATION
# ============================================================================

# Extract spatial field posterior draws
spatial_field_draws <- fit$draws("spatial_field", format = "draws_matrix")
# Dimensions: [iterations × n_grids_total]

# Calculate posterior summaries
spatial_mean <- colMeans(spatial_field_draws)
spatial_sd <- apply(spatial_field_draws, 2, sd)
spatial_lower <- apply(spatial_field_draws, 2, quantile, probs = 0.025)
spatial_upper <- apply(spatial_field_draws, 2, quantile, probs = 0.975)

# Create a data frame with spatial information
spatial_df <- data.frame(
  grid_id = 1:stan_data$n_grids_total,
  easting = stan_data$grid_coords[, 1],
  northing = stan_data$grid_coords[, 2],
  spatial_mean = spatial_mean,
  spatial_sd = spatial_sd,
  spatial_lower = spatial_lower,
  spatial_upper = spatial_upper
)

# Check range of spatial effects
cat("Spatial field summary:\n")
cat(sprintf("Mean: %.3f [%.3f, %.3f]\n", 
            mean(spatial_mean), min(spatial_mean), max(spatial_mean)))
cat(sprintf("SD: %.3f [%.3f, %.3f]\n",
            mean(spatial_sd), min(spatial_sd), max(spatial_sd)))


library(ggplot2)
library(viridis)
library(sf)
library(sp)

# 1. Mean spatial field
p1 <- ggplot(spatial_df, aes(x = easting, y = northing, fill = spatial_mean)) +
  geom_tile(width = diff(range(spatial_df$easting))/50, 
            height = diff(range(spatial_df$northing))/50) +
  scale_fill_viridis(option = "viridis", name = "Spatial Effect") +
  labs(title = "Posterior Mean Spatial Field",
       subtitle = "GP effect on log-abundance") +
  theme_minimal() +
  theme(legend.position = "right")

# 2. Uncertainty (SD) of spatial field
p2 <- ggplot(spatial_df, aes(x = easting, y = northing, fill = spatial_sd)) +
  geom_tile(width = diff(range(spatial_df$easting))/50, 
            height = diff(range(spatial_df$northing))/50) +
  scale_fill_viridis(option = "plasma", name = "Uncertainty (SD)") +
  labs(title = "Posterior Uncertainty in Spatial Field") +
  theme_minimal() +
  theme(legend.position = "right")

# 3. Credible interval width
spatial_df$ci_width <- spatial_upper - spatial_lower
p3 <- ggplot(spatial_df, aes(x = easting, y = northing, fill = ci_width)) +
  geom_tile(width = diff(range(spatial_df$easting))/50, 
            height = diff(range(spatial_df$northing))/50) +
  scale_fill_viridis(option = "magma", name = "95% CI Width") +
  labs(title = "Width of 95% Credible Intervals") +
  theme_minimal() +
  theme(legend.position = "right")

# Arrange plots
library(gridExtra)
grid.arrange(p1, p2, p3, ncol = 2)


######### Spatial Correlation Analysis

# Extract GP hyperparameters
if ("rho_gp" %in% fit$metadata()$model_params) {
  rho_samples <- as.numeric(fit$draws("rho_gp", format = "draws_matrix"))
  alpha_samples <- as.numeric(fit$draws("alpha_gp", format = "draws_matrix"))
  
  cat("\nGP Hyperparameter Summary:\n")
  cat(sprintf("Range parameter (rho): %.2f [%.2f, %.2f]\n",
              mean(rho_samples), quantile(rho_samples, 0.025), 
              quantile(rho_samples, 0.975)))
  cat(sprintf("Marginal SD (alpha): %.3f [%.3f, %.3f]\n",
              mean(alpha_samples), quantile(alpha_samples, 0.025),
              quantile(alpha_samples, 0.975)))
  
  # Plot hyperparameter posteriors
  par(mfrow = c(1, 2))
  hist(rho_samples, breaks = 30, col = "steelblue", main = "Posterior: Range (rho)",
       xlab = "Spatial range (km)", border = "white")
  abline(v = mean(rho_samples), col = "red", lwd = 2)
  
  hist(alpha_samples, breaks = 30, col = "darkgreen", main = "Posterior: Marginal SD (alpha)",
       xlab = "Spatial standard deviation", border = "white")
  abline(v = mean(alpha_samples), col = "red", lwd = 2)
  par(mfrow = c(1, 1))
}

# Calculate empirical spatial correlation
calculate_spatial_correlation <- function(coords, spatial_effects_matrix) {
  # spatial_effects_matrix: rows = samples/draws, cols = grid cells
  # coords: n_grids_total × 2 matrix
  
  n_grids <- nrow(coords)
  n_draws <- nrow(spatial_effects_matrix)
  
  # Calculate distance matrix
  distances <- as.matrix(dist(coords))
  
  # Bin distances
  max_dist <- max(distances)
  bins <- seq(0, max_dist, length.out = 20)
  bin_centers <- (bins[-1] + bins[-length(bins)]) / 2
  
  # Initialize results
  cor_by_dist <- data.frame(
    distance = bin_centers,
    correlation = rep(NA, length(bin_centers)),
    correlation_sd = rep(NA, length(bin_centers)),
    n_pairs = rep(0, length(bin_centers))
  )
  
  # For each distance bin
  for (i in 1:(length(bins)-1)) {
    # Find pairs in this distance range
    pair_indices <- which(distances > bins[i] & distances <= bins[i+1] & 
                            upper.tri(distances), arr.ind = TRUE)
    
    if (nrow(pair_indices) > 0) {
      # Calculate correlations for each pair across draws
      pair_correlations <- numeric(nrow(pair_indices))
      
      for (j in 1:nrow(pair_indices)) {
        idx1 <- pair_indices[j, 1]
        idx2 <- pair_indices[j, 2]
        
        # Extract the time series for these two grid cells
        ts1 <- spatial_effects_matrix[, idx1]
        ts2 <- spatial_effects_matrix[, idx2]
        
        # Calculate correlation
        pair_correlations[j] <- cor(ts1, ts2)
      }
      
      # Store results
      cor_by_dist$correlation[i] <- mean(pair_correlations, na.rm = TRUE)
      cor_by_dist$correlation_sd[i] <- sd(pair_correlations, na.rm = TRUE)
      cor_by_dist$n_pairs[i] <- length(pair_correlations)
    }
  }
  
  # Remove bins with no pairs
  cor_by_dist <- cor_by_dist[cor_by_dist$n_pairs > 0, ]
  
  return(cor_by_dist)
}


# Calculate for posterior mean
cor_analysis <- calculate_spatial_correlation(
  stan_data$grid_coords,
  matrix(spatial_mean, nrow = 1)  # Just the mean
)

# Plot empirical correlation function
ggplot(cor_analysis, aes(x = distance, y = correlation)) +
  geom_point(aes(size = n_pairs), alpha = 0.6) +
  geom_smooth(method = "loess", color = "red", se = FALSE) +
  labs(title = "Empirical Spatial Correlation",
       x = "Distance (km)", y = "Correlation") +
  theme_minimal()









# ============================================================================
# MODEL DIAGNOSTICS
# ============================================================================


run_model_diagnostics <- function(fit, stan_data) {
  
  cat("MODEL DIAGNOSTICS REPORT\n")
  cat("========================\n\n")
  
  # 1. Check for divergent transitions
  if ("divergent__" %in% fit$metadata()$model_params) {
    divergences <- fit$draws("divergent__", format = "draws_matrix")
    n_divergent <- sum(divergences)
    cat(sprintf("1. Divergent transitions: %d\n", n_divergent))
    if (n_divergent > 0) {
      cat("   ⚠️  Consider increasing adapt_delta\n")
    } else {
      cat("   ✓ No divergent transitions\n")
    }
  }
  
  # 2. R-hat diagnostics
  rhats <- fit$summary()$rhat
  cat(sprintf("\n2. R-hat diagnostics:\n"))
  cat(sprintf("   Max R-hat: %.3f\n", max(rhats, na.rm = TRUE)))
  cat(sprintf("   Parameters with R-hat > 1.01: %d\n", sum(rhats > 1.01, na.rm = TRUE)))
  
  # 3. ESS diagnostics
  ess_bulk <- fit$summary()$ess_bulk
  ess_tail <- fit$summary()$ess_tail
  cat(sprintf("\n3. ESS diagnostics:\n"))
  cat(sprintf("   Min bulk ESS: %.0f\n", min(ess_bulk, na.rm = TRUE)))
  cat(sprintf("   Min tail ESS: %.0f\n", min(ess_tail, na.rm = TRUE)))
  
  # 4. N_max adequacy
  if ("max_N_ratio" %in% fit$metadata()$model_params) {
    max_N_ratio <- fit$draws("max_N_ratio", format = "draws_matrix")
    n_near_boundary <- fit$draws("n_near_boundary", format = "draws_matrix")
    
    cat(sprintf("\n4. N_max adequacy:\n"))
    cat(sprintf("   Max N/N_max ratio: %.3f\n", max(max_N_ratio)))
    cat(sprintf("   Mean N/N_max ratio: %.3f\n", mean(max_N_ratio)))
    cat(sprintf("   Draws near boundary: %.1f%%\n", 
                mean(n_near_boundary > 0) * 100))
    
    if (mean(max_N_ratio) > 0.8) {
      cat("   ⚠️  Consider increasing N_multiplier\n")
    }
  }
  
  # 5. PPC summary
  if ("y_rep" %in% fit$metadata()$model_params) {
    y_rep <- fit$draws("y_rep", format = "draws_matrix")
    y_obs <- stan_data$y
    
    # Calculate Bayes p-values
    calc_bayes_p <- function(stat_func) {
      stat_obs <- stat_func(y_obs)
      stat_rep <- apply(y_rep, 1, stat_func)
      mean(stat_rep >= stat_obs)
    }
    
    p_mean <- calc_bayes_p(mean)
    p_sd <- calc_bayes_p(sd)
    p_max <- calc_bayes_p(max)
    p_zeros <- calc_bayes_p(function(x) mean(x == 0))
    
    cat(sprintf("\n5. Posterior Predictive Checks:\n"))
    cat(sprintf("   Mean: p = %.3f\n", p_mean))
    cat(sprintf("   SD: p = %.3f\n", p_sd))
    cat(sprintf("   Max: p = %.3f\n", p_max))
    cat(sprintf("   Zeros: p = %.3f\n", p_zeros))
    
    if (any(c(p_mean, p_sd, p_max, p_zeros) < 0.05 | 
            c(p_mean, p_sd, p_max, p_zeros) > 0.95)) {
      cat("   ⚠️  Some statistics show poor fit\n")
    }
  }
  
  cat("\n6. Parameter estimates summary:\n")
  print(fit$summary(variables = c("beta0", "beta_temp", "beta_rain", 
                                  "alpha0", "phi", "alpha_gp", "rho_gp")))
}

# Run diagnostics
run_model_diagnostics(fit, stan_data)
