# ============================================================================
# GAMMA PARAMETER CONVERGENCE DIAGNOSTICS & IMPUTATION VALIDATION
# ============================================================================

library(tidyverse)
library(bayesplot)
library(posterior)
library(patchwork)

# Assuming your fitted model is 'fit'
# draws <- as_draws_df(fit)

# ============================================================================
# 1. COMPREHENSIVE CONVERGENCE CHECK FOR GAMMA PARAMETERS
# ============================================================================

cat("=== GAMMA PARAMETER CONVERGENCE DIAGNOSTICS ===\n\n")

# Extract all gamma and sigma parameters
gamma_params <- draws %>%
  select(starts_with("gamma_"), starts_with("sigma_"))

# Get summary with Rhat, ESS
gamma_summary <- summarise_draws(
  gamma_params,
  "mean", "sd", "rhat", "ess_bulk", "ess_tail"
) %>%
  arrange(desc(rhat))

# Print parameters with Rhat > 1.1
high_rhat <- gamma_summary %>% filter(rhat > 1.1)
print(high_rhat)

# Print parameters with Rhat > 1.05
moderate_rhat <- gamma_summary %>% filter(rhat > 1.05, rhat <= 1.1)
print(moderate_rhat)

# Summary statistics
cat("Min:", min(gamma_summary$rhat), "\n")
cat("Median:", median(gamma_summary$rhat), "\n")
cat("Max:", max(gamma_summary$rhat), "\n")
cat("% with Rhat > 1.05:", 
    round(100 * mean(gamma_summary$rhat > 1.05), 1), "%\n")
cat("% with Rhat > 1.1:", 
    round(100 * mean(gamma_summary$rhat > 1.1), 1), "%\n\n")

# ============================================================================
# 2. VISUALIZE RHAT DISTRIBUTION BY PARAMETER TYPE
# ============================================================================

# Categorize parameters
gamma_summary_categorized <- gamma_summary %>%
  mutate(
    param_type = case_when(
      str_detect(variable, "gamma_temp_land") ~ "Temp: Land Effects",
      str_detect(variable, "gamma_rain_land") ~ "Rain: Land Effects",
      str_detect(variable, "gamma_temp_time") ~ "Temp: Time Effects",
      str_detect(variable, "gamma_rain_time") ~ "Rain: Time Effects",
      str_detect(variable, "gamma_temp0") ~ "Temp: Intercept",
      str_detect(variable, "gamma_rain0") ~ "Rain: Intercept",
      str_detect(variable, "sigma_temp") ~ "Temp: Variance",
      str_detect(variable, "sigma_rain") ~ "Rain: Variance",
      TRUE ~ "Other"
    )
  )

# Plot Rhat by parameter type
p_rhat <- ggplot(gamma_summary_categorized, aes(x = rhat, fill = param_type)) +
  geom_histogram(bins = 40, alpha = 0.7) +
  geom_vline(xintercept = 1.01, linetype = "dashed", color = "orange", linewidth = 1) +
  geom_vline(xintercept = 1.05, linetype = "dashed", color = "red", linewidth = 1) +
  facet_wrap(~param_type, scales = "free_y") +
  labs(title = "Rhat Distribution by Parameter Type",
       subtitle = "Orange = 1.01 (ideal), Red = 1.05 (concerning)",
       x = "Rhat", y = "Count") +
  theme_minimal() +
  theme(legend.position = "none")

print(p_rhat)

# ============================================================================
# 3. ESS DIAGNOSTICS (Low ESS often accompanies high Rhat)
# ============================================================================

# Plot ESS
p_ess <- gamma_summary_categorized %>%
  pivot_longer(cols = c(ess_bulk, ess_tail), 
               names_to = "ess_type", values_to = "ess") %>%
  ggplot(aes(x = ess, fill = ess_type)) +
  geom_histogram(bins = 40, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = 400, linetype = "dashed", color = "red", linewidth = 1) +
  facet_wrap(~param_type, scales = "free") +
  labs(title = "Effective Sample Size by Parameter Type",
       subtitle = "Red line = 400 (minimum recommended)",
       x = "ESS", y = "Count") +
  theme_minimal()

print(p_ess)

# Identify parameters with both high Rhat AND low ESS
problematic <- gamma_summary %>%
  filter(rhat > 1.05 | ess_bulk < 400) %>%
  arrange(desc(rhat))

#Parameters with convergence issues (Rhat > 1.05 OR ESS < 400))
print(problematic)

# ============================================================================
# 4. TRACE PLOTS FOR PROBLEMATIC PARAMETERS
# ============================================================================

if (nrow(moderate_rhat) > 0) {
  
  # Plot up to 6 worst parameters
  n_plot <- min(6, nrow(moderate_rhat))
  worst_params <- moderate_rhat$variable[1:n_plot]
  
  trace_plots <- mcmc_trace(fit, pars = worst_params, 
                            facet_args = list(ncol = 2))
  print(trace_plots)
  
  # Pairs plot for correlation diagnosis
  if (n_plot >= 2) {
    cat("\nPairs plot (checking for correlation issues):\n")
    pairs_plot <- mcmc_pairs(fit, pars = worst_params[1:min(4, n_plot)])
    print(pairs_plot)
  }
}

# ============================================================================
# 5. TEMPORAL PATTERN IN CONVERGENCE (Time effects)
# ============================================================================

# Check if time effects show patterns in convergence
temp_time_rhat <- gamma_summary %>%
  filter(str_detect(variable, "gamma_temp_time\\[")) %>%
  mutate(time_index = as.integer(str_extract(variable, "\\d+")))

rain_time_rhat <- gamma_summary %>%
  filter(str_detect(variable, "gamma_rain_time\\[")) %>%
  mutate(time_index = as.integer(str_extract(variable, "\\d+")))

if (nrow(temp_time_rhat) > 0) {
  p_time_convergence <- bind_rows(
    temp_time_rhat %>% mutate(type = "Temperature"),
    rain_time_rhat %>% mutate(type = "Precipitation")
  ) %>%
    ggplot(aes(x = time_index, y = rhat, color = type)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_line(alpha = 0.4) +
    geom_hline(yintercept = 1.05, linetype = "dashed", color = "red") +
    facet_wrap(~type) +
    labs(title = "Convergence of Temporal Effects Over Time",
         subtitle = "Are certain time periods problematic?",
         x = "Time Index", y = "Rhat") +
    theme_minimal()
  
  print(p_time_convergence)
}

# ============================================================================
# 6. VALIDATION: DO IMPUTED VALUES MATCH OBSERVED DISTRIBUTION?
# ============================================================================

cat("\n=== IMPUTATION QUALITY VALIDATION ===\n\n")

# Extract complete temperature and rain matrices from posterior
# (This requires reconstructing them from draws)

# Get posterior means of complete z_temp and z_rain
# First, extract observed and imputed separately

# Observed values (from your data)
temp_obs_values <- as.vector(z_temp_obs)[!is.na(as.vector(z_temp_obs))]
rain_obs_values <- as.vector(z_rain_obs)[!is.na(as.vector(z_rain_obs))]

# Imputed values (posterior means)
temp_imputed_draws <- draws %>% select(starts_with("z_temp_imputed_raw["))
rain_imputed_draws <- draws %>% select(starts_with("z_rain_imputed_raw["))

temp_imputed_values <- colMeans(temp_imputed_draws)
rain_imputed_values <- colMeans(rain_imputed_draws)

# --- Visual Comparison: Q-Q Plots ---
par(mfrow = c(2, 2))

# Temperature Q-Q plot
qqplot(temp_obs_values, temp_imputed_values,
       main = "Q-Q Plot: Temperature",
       xlab = "Observed Quantiles",
       ylab = "Imputed Quantiles",
       pch = 16, col = rgb(0, 0, 1, 0.3))
abline(0, 1, col = "red", lwd = 2)

# Rain Q-Q plot
qqplot(rain_obs_values, rain_imputed_values,
       main = "Q-Q Plot: Precipitation",
       xlab = "Observed Quantiles",
       ylab = "Imputed Quantiles",
       pch = 16, col = rgb(0, 0, 1, 0.3))
abline(0, 1, col = "red", lwd = 2)

# Temperature density comparison
plot(density(temp_obs_values), main = "Temperature Distribution",
     xlab = "Standardized Value", col = "blue", lwd = 2,
     xlim = range(c(temp_obs_values, temp_imputed_values)))
lines(density(temp_imputed_values), col = "red", lwd = 2, lty = 2)
legend("topright", c("Observed", "Imputed"), 
       col = c("blue", "red"), lwd = 2, lty = c(1, 2))

# Rain density comparison
plot(density(rain_obs_values), main = "Precipitation Distribution",
     xlab = "Standardized Value", col = "blue", lwd = 2,
     xlim = range(c(rain_obs_values, rain_imputed_values)))
lines(density(rain_imputed_values), col = "red", lwd = 2, lty = 2)
legend("topright", c("Observed", "Imputed"), 
       col = c("blue", "red"), lwd = 2, lty = c(1, 2))

par(mfrow = c(1, 1))

# --- Statistical Tests ---
cat("\n--- Kolmogorov-Smirnov Tests ---\n")
cat("Temperature:\n")
ks_temp <- ks.test(temp_obs_values, temp_imputed_values)
print(ks_temp)
cat("\nPrecipitation:\n")
ks_rain <- ks.test(rain_obs_values, rain_imputed_values)
print(ks_rain)

# --- Summary Statistics Comparison ---
comparison_table <- tibble(
  Variable = c("Temperature", "Precipitation"),
  Obs_Mean = c(mean(temp_obs_values), mean(rain_obs_values)),
  Imp_Mean = c(mean(temp_imputed_values), mean(rain_imputed_values)),
  Obs_SD = c(sd(temp_obs_values), sd(rain_obs_values)),
  Imp_SD = c(sd(temp_imputed_values), sd(rain_imputed_values)),
  Obs_Min = c(min(temp_obs_values), min(rain_obs_values)),
  Imp_Min = c(min(temp_imputed_values), min(rain_imputed_values)),
  Obs_Max = c(max(temp_obs_values), max(rain_obs_values)),
  Imp_Max = c(max(temp_imputed_values), max(rain_imputed_values))
)

cat("\n--- Summary Statistics Comparison ---\n")
print(comparison_table)

# ============================================================================
# 7. SPATIAL/TEMPORAL PATTERNS IN IMPUTATION QUALITY
# ============================================================================

# Create a data frame mapping imputed values to their grid/date locations
# (This requires knowing which grid-date combinations were missing)

# Reconstruct missing indices
missing_info <- expand_grid(
  grid = 1:n_grids_total,
  date = 1:n_dates
) %>%
  mutate(
    temp_missing = as.logical(temp_missing[cbind(grid, date)]),
    rain_missing = as.logical(rain_missing[cbind(grid, date)])
  )

temp_missing_info <- missing_info %>%
  filter(temp_missing) %>%
  mutate(
    imputed_mean = temp_imputed_values[row_number()],
    imputed_sd = apply(temp_imputed_draws, 2, sd)[row_number()]
  )

rain_missing_info <- missing_info %>%
  filter(rain_missing) %>%
  mutate(
    imputed_mean = rain_imputed_values[row_number()],
    imputed_sd = apply(rain_imputed_draws, 2, sd)[row_number()]
  )

# Plot uncertainty by time
p_temp_uncertainty <- ggplot(temp_missing_info, aes(x = date, y = imputed_sd)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "blue") +
  labs(title = "Temperature Imputation Uncertainty Over Time",
       x = "Date Index", y = "Posterior SD") +
  theme_minimal()

p_rain_uncertainty <- ggplot(rain_missing_info, aes(x = date, y = imputed_sd)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "blue") +
  labs(title = "Precipitation Imputation Uncertainty Over Time",
       x = "Date Index", y = "Posterior SD") +
  theme_minimal()

print(p_temp_uncertainty / p_rain_uncertainty)

# ============================================================================
# 8. RECOMMENDATIONS BASED ON DIAGNOSTICS
# ============================================================================

cat("\n=== DIAGNOSTIC SUMMARY & RECOMMENDATIONS ===\n\n")

# Count issues
n_high_rhat <- sum(gamma_summary$rhat > 1.1)
n_low_ess <- sum(gamma_summary$ess_bulk < 400)
n_issues <- sum(gamma_summary$rhat > 1.05 | gamma_summary$ess_bulk < 400)

cat("Issue Summary:\n")
cat("- Parameters with Rhat > 1.1:", n_high_rhat, "\n")
cat("- Parameters with ESS < 400:", n_low_ess, "\n")
cat("- Parameters with any convergence issue:", n_issues, "\n\n")

if (n_high_rhat > 0) {
  cat("RECOMMENDATIONS:\n\n")
  
  cat("1. INCREASE SAMPLING:\n")
  cat("   iter = 4000, warmup = 2000, chains = 4\n\n")
  
  cat("2. INCREASE adapt_delta:\n")
  cat("   control = list(adapt_delta = 0.99, max_treedepth = 12)\n\n")
  
  cat("3. USE NON-CENTERED PARAMETERIZATION:\n")
  cat("   Switch to the improved model (already done for z_imputed)\n")
  cat("   Consider for gamma_*_time parameters if they're problematic\n\n")
  
  cat("4. STRONGER PRIORS:\n")
  cat("   gamma_*_time_raw ~ normal(0, 0.5)  # instead of normal(0, sigma_*_time)\n")
  cat("   sigma_*_time ~ normal(0.3, 0.15)   # more informative\n\n")
  
  # Identify which parameter types are most problematic
  worst_types <- gamma_summary_categorized %>%
    group_by(param_type) %>%
    summarise(
      max_rhat = max(rhat),
      mean_rhat = mean(rhat),
      n_problems = sum(rhat > 1.05)
    ) %>%
    arrange(desc(max_rhat))
  
  cat("5. MOST PROBLEMATIC PARAMETER TYPES:\n")
  print(worst_types)
  cat("\n")
  
  if (any(str_detect(high_rhat$variable, "time"))) {
    cat("6. TEMPORAL EFFECTS ISSUE DETECTED:\n")
    cat("   Consider:\n")
    cat("   - Using a random walk prior: gamma_*_time[t] ~ normal(gamma_*_time[t-1], sigma)\n")
    cat("   - Adding AR(1) structure\n")
    cat("   - Using a spline instead of independent time effects\n\n")
  }
}

# Imputation quality assessment
cat("IMPUTATION QUALITY:\n")
if (ks_temp$p.value < 0.05 | ks_rain$p.value < 0.05) {
  cat("⚠️  WARNING: Imputed and observed distributions differ significantly\n")
  cat("   - Check if imputation model includes all relevant predictors\n")
  cat("   - Consider adding spatial effects or other covariates\n\n")
} else {
  cat("✓ Imputed values match observed distribution well\n\n")
}

# ============================================================================
# 9. EXPORT DETAILED REPORT
# ============================================================================

convergence_report <- list(
  gamma_summary = gamma_summary,
  high_rhat_params = high_rhat,
  moderate_rhat_params = moderate_rhat,
  low_ess_params = gamma_summary %>% filter(ess_bulk < 400),
  imputation_comparison = comparison_table,
  ks_tests = list(temperature = ks_temp, precipitation = ks_rain),
  temporal_convergence = list(
    temperature = temp_time_rhat,
    precipitation = rain_time_rhat
  )
)

saveRDS(convergence_report, "gamma_convergence_report.rds")

cat("\nDetailed report saved to: gamma_convergence_report.rds\n")
cat("=== DIAGNOSTICS COMPLETE ===\n")