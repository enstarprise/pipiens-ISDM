# ============================================================================
# POSTERIOR PREDICTIVE CHECKS AND IMPUTATION DIAGNOSTICS
# ============================================================================

library(tidyverse)
library(bayesplot)
library(posterior)

source("02_model/stan/01_prepare_stan_data_impute.R")

# ============================================================================
# 1. POSTERIOR PREDICTIVE CHECKS - SURVEY DATA
# ============================================================================

# Extract PPC summary statistics
ppc_stats <- draws %>%
  select(mean_y_obs, sd_y_obs, max_y_obs, prop_zeros_obs,
         mean_y_rep, sd_y_rep, max_y_rep, prop_zeros_rep)

# --- GRAPHICAL PPC: Summary Statistics ---
ppc_plots <- list()

# Mean comparison
ppc_plots$mean <- ggplot(ppc_stats) +
  geom_density(aes(x = mean_y_rep), fill = "skyblue", alpha = 0.6) +
  geom_vline(aes(xintercept = mean_y_obs[1]), color = "red", linewidth = 1.5) +
  labs(title = "PPC: Mean of counts",
       x = "Mean", y = "Density",
       subtitle = "Red line = observed data") +
  theme_minimal()

# Standard deviation comparison
ppc_plots$sd <- ggplot(ppc_stats) +
  geom_density(aes(x = sd_y_rep), fill = "skyblue", alpha = 0.6) +
  geom_vline(aes(xintercept = sd_y_obs[1]), color = "red", linewidth = 1.5) +
  labs(title = "PPC: SD of counts",
       x = "Standard Deviation", y = "Density") +
  theme_minimal()

# Maximum comparison
ppc_plots$max <- ggplot(ppc_stats) +
  geom_histogram(aes(x = max_y_rep), fill = "skyblue", alpha = 0.6, bins = 30) +
  geom_vline(aes(xintercept = max_y_obs[1]), color = "red", linewidth = 1.5) +
  labs(title = "PPC: Maximum count",
       x = "Maximum", y = "Count") +
  theme_minimal()

# Proportion of zeros comparison
ppc_plots$zeros <- ggplot(ppc_stats) +
  geom_density(aes(x = prop_zeros_rep), fill = "skyblue", alpha = 0.6) +
  geom_vline(aes(xintercept = prop_zeros_obs[1]), color = "red", linewidth = 1.5) +
  labs(title = "PPC: Proportion of zeros",
       x = "Proportion zeros", y = "Density") +
  theme_minimal()

# Display all PPC plots
library(patchwork)
(ppc_plots$mean + ppc_plots$sd) / (ppc_plots$max + ppc_plots$zeros)

# --- NUMERIC PPC: Bayesian p-values ---
bayesian_pvals <- tibble(
  statistic = c("mean", "sd", "max", "prop_zeros"),
  p_value = c(
    mean(ppc_stats$mean_y_rep >= ppc_stats$mean_y_obs[1]),
    mean(ppc_stats$sd_y_rep >= ppc_stats$sd_y_obs[1]),
    mean(ppc_stats$max_y_rep >= ppc_stats$max_y_obs[1]),
    mean(ppc_stats$prop_zeros_rep >= ppc_stats$prop_zeros_obs[1])
  )
)

print("Bayesian p-values (should be between 0.05 and 0.95):")
print(bayesian_pvals)

# ============================================================================
# 2. POSTERIOR PREDICTIVE CHECKS - PRESENCE-ONLY DATA
# this won't really make sense for a model that is only indexed by the observed
# grid cells, and not the full spatial structure
# ============================================================================

# Extract PO replications
po_rep_draws <- draws %>%
  select(starts_with("po_rep["))

# Calculate proportion of presences in replications
prop_presence_rep <- rowMeans(po_rep_draws)

# Observed proportion (should be 1.0 since all are presences)
prop_presence_obs <- 1.0

# Plot
ggplot() +
  geom_density(aes(x = prop_presence_rep), fill = "skyblue", alpha = 0.6) +
  geom_vline(xintercept = prop_presence_obs, color = "red", linewidth = 1.5) +
  labs(title = "PPC: Presence-Only Data",
       subtitle = "Proportion of simulated presences",
       x = "Proportion", y = "Density") +
  theme_minimal()

# ============================================================================
# 3. IMPUTATION DIAGNOSTICS - Rhat and Convergence
# ============================================================================

# Extract imputed parameters
temp_imputed <- draws %>% select(starts_with("z_temp_imputed_raw["))
rain_imputed <- draws %>% select(starts_with("z_rain_imputed_raw["))

# Calculate Rhat for imputed values
library(posterior)
rhat_temp <- summarise_draws(temp_imputed, "rhat")
rhat_rain <- summarise_draws(rain_imputed, "rhat")

# Plot Rhat distributions
rhat_plot <- bind_rows(
  rhat_temp %>% mutate(variable_type = "Temperature"),
  rhat_rain %>% mutate(variable_type = "Precipitation")
) %>%
  ggplot(aes(x = rhat, fill = variable_type)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = 1.02, linetype = "dashed", color = "orange", linewidth = 1) +
  geom_vline(xintercept = 1.05, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Rhat values for imputed observations",
       subtitle = "Orange line = 1.01 (ideal), Red line = 1.05 (concerning)",
       x = "Rhat", y = "Count") +
  theme_minimal() +
  facet_wrap(~variable_type)

print(rhat_plot)

# Summary statistics for Rhat
print(summary(rhat_temp$rhat))
cat("\nNumber with Rhat > 1.05:", sum(rhat_temp$rhat > 1.05), "\n")
cat("Number with Rhat > 1.1:", sum(rhat_temp$rhat > 1.1), "\n")

print(summary(rhat_rain$rhat))
cat("\nNumber with Rhat > 1.05:", sum(rhat_rain$rhat > 1.05), "\n")
cat("Number with Rhat > 1.1:", sum(rhat_rain$rhat > 1.1), "\n")

# ============================================================================
# 4. IMPUTATION QUALITY CHECKS
# ============================================================================

# --- Check 1: Compare imputed vs observed distributions ---

# Get posterior means of imputed values
temp_imputed_means <- colMeans(temp_imputed)
rain_imputed_means <- colMeans(rain_imputed)

# Get observed values (need to extract from your data)
# Assuming you have z_temp_stan and z_rain_stan from your data prep
temp_observed <- as.vector(z_temp_stan)[!is.na(as.vector(z_temp_stan))]
rain_observed <- as.vector(z_rain_stan)[!is.na(as.vector(z_rain_stan))]

# Compare distributions
imputation_comparison <- bind_rows(
  tibble(value = temp_observed, type = "Observed", variable = "Temperature"),
  tibble(value = temp_imputed_means, type = "Imputed", variable = "Temperature"),
  tibble(value = rain_observed, type = "Observed", variable = "Precipitation"),
  tibble(value = rain_imputed_means, type = "Imputed", variable = "Precipitation")
)

ggplot(imputation_comparison, aes(x = value, fill = type)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~variable, scales = "free") +
  labs(title = "Imputed vs Observed Value Distributions",
       x = "Standardized Value", y = "Density") +
  theme_minimal()

# --- Check 2: Uncertainty in imputed values ---

# Calculate posterior SD for each imputed value
temp_imputed_sd <- apply(temp_imputed, 2, sd)
rain_imputed_sd <- apply(rain_imputed, 2, sd)

# Plot uncertainty
uncertainty_plot <- bind_rows(
  tibble(sd = temp_imputed_sd, variable = "Temperature"),
  tibble(sd = rain_imputed_sd, variable = "Precipitation")
) %>%
  ggplot(aes(x = sd)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  facet_wrap(~variable, scales = "free") +
  labs(title = "Posterior Uncertainty in Imputed Values",
       subtitle = "Higher SD = more uncertain imputation",
       x = "Posterior SD", y = "Count") +
  theme_minimal()

print(uncertainty_plot)

# --- Check 3: Effective sample size for imputed values ---
ess_temp <- summarise_draws(temp_imputed, "ess_bulk", "ess_tail")
ess_rain <- summarise_draws(rain_imputed, "ess_bulk", "ess_tail")

ess_plot <- bind_rows(
  ess_temp %>% mutate(variable_type = "Temperature"),
  ess_rain %>% mutate(variable_type = "Precipitation")
) %>%
  pivot_longer(cols = c(ess_bulk, ess_tail), 
               names_to = "ess_type", values_to = "ess") %>%
  ggplot(aes(x = ess, fill = ess_type)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = 400, linetype = "dashed", color = "red") +
  facet_wrap(~variable_type) +
  labs(title = "Effective Sample Size for Imputed Values",
       subtitle = "Red line = 400 (minimum recommended)",
       x = "ESS", y = "Count") +
  theme_minimal()

print(ess_plot)

# ============================================================================
# 5. DIAGNOSE HIGH RHAT VALUES (if Rhat > 1.1)
# ============================================================================

if (any(rhat_temp$rhat > 1.1) || any(rhat_rain$rhat > 1.1)) {
  cat("\n=== DIAGNOSING HIGH RHAT VALUES ===\n")
  
  # Identify problematic parameters
  problem_temp <- rhat_temp %>% filter(rhat > 1.1)
  problem_rain <- rhat_rain %>% filter(rhat > 1.1)
  
  cat("\nProblematic temperature imputations:", nrow(problem_temp), "\n")
  cat("Problematic rain imputations:", nrow(problem_rain), "\n")
  
  # Check trace plots for a few problematic parameters
  if (nrow(problem_temp) > 0) {
    # Extract first problematic parameter
    param_name <- problem_temp$variable[1]
    
    # Trace plot
    trace_plot <- mcmc_trace(fit, pars = param_name)
    print(trace_plot)
    
    cat("\nRECOMMENDATIONS:\n")
    cat("1. Increase number of iterations (especially warmup)\n")
    cat("2. Increase adapt_delta: control = list(adapt_delta = 0.95)\n")
    cat("3. Check if imputation model is well-specified\n")
    cat("4. Consider stronger priors on imputation parameters\n")
  }
}

# ============================================================================
# 6. IMPUTATION MODEL PARAMETER CHECKS
# ============================================================================

# Extract imputation model parameters
imputation_params <- draws %>%
  select(gamma_temp0, gamma_rain0, 
         sigma_temp, sigma_rain,
         sigma_temp_time, sigma_rain_time,
         starts_with("gamma_temp_land"),
         starts_with("gamma_rain_land"))

# Summary
cat("\n=== Imputation Model Parameter Summary ===\n")
print(summarise_draws(imputation_params, "mean", "sd", "rhat", "ess_bulk"))

# Check if temporal effects are appropriate
gamma_temp_time <- draws %>% select(starts_with("gamma_temp_time["))
gamma_rain_time <- draws %>% select(starts_with("gamma_rain_time["))

# Plot temporal effects
temporal_effects <- bind_rows(
  tibble(
    date = 1:ncol(gamma_temp_time),
    mean = colMeans(gamma_temp_time),
    lower = apply(gamma_temp_time, 2, quantile, 0.025),
    upper = apply(gamma_temp_time, 2, quantile, 0.975),
    variable = "Temperature"
  ),
  tibble(
    date = 1:ncol(gamma_rain_time),
    mean = colMeans(gamma_rain_time),
    lower = apply(gamma_rain_time, 2, quantile, 0.025),
    upper = apply(gamma_rain_time, 2, quantile, 0.975),
    variable = "Precipitation"
  )
)

ggplot(temporal_effects, aes(x = date, y = mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  facet_wrap(~variable, scales = "free_y") +
  labs(title = "Temporal Effects in Imputation Model",
       x = "Date Index", y = "Effect Size") +
  theme_minimal()

# ============================================================================
# 7. EXPORT SUMMARY REPORT
# ============================================================================

report <- list(
  ppc_bayesian_pvals = bayesian_pvals,
  rhat_temp_summary = summary(rhat_temp$rhat),
  rhat_rain_summary = summary(rhat_rain$rhat),
  n_high_rhat_temp = sum(rhat_temp$rhat > 1.1),
  n_high_rhat_rain = sum(rhat_rain$rhat > 1.1),
  imputation_params_summary = summarise_draws(imputation_params)
)

# Save report
saveRDS(report, "model_diagnostics_report.rds")

cat("\n=== DIAGNOSTICS COMPLETE ===\n")
cat("All plots and summaries have been generated.\n")
cat("Report saved to: model_diagnostics_report.rds\n")