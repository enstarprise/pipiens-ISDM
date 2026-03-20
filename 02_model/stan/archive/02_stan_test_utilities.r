# ============================================================================
# TESTING UTILITIES FOR STAN MARGINALIZED MODEL
# ============================================================================

library(tidyverse)
library(cmdstanr)

# ============================================================================
# 1. SIMULATE SMALL TEST DATA
# ============================================================================

simulate_test_data <- function(
    n_grids = 10,
    n_dates = 5,
    n_obs_y = 50,
    n_obs_po = 100,
    n_land_covs = 3,
    seed = 123
) {
  set.seed(seed)
  
  cat("Simulating test data...\n")
  
  # True parameters
  beta0 <- -1
  beta_temp <- 0.5
  beta_rain <- -0.3
  # beta_buildings <- 0.2
  beta_land <- rnorm(n_land_covs, 0, 0.3)
  phi <- 2
  
  # Covariates
  z_temp <- matrix(rnorm(n_grids * n_dates), n_grids, n_dates)
  z_rain <- matrix(rnorm(n_grids * n_dates), n_grids, n_dates)
  # z_buildings <- rnorm(n_grids)
  z_land <- matrix(rnorm(n_grids * n_land_covs), n_grids, n_land_covs)
  z_poi <- rnorm(n_grids)
  z_reports <- rnorm(n_grids)
  area_grid <- rep(1, n_grids)  # Simplified
  
  # Generate lambda for each grid-date
  lambda_grid <- matrix(NA, n_grids, n_dates)
  for (g in 1:n_grids) {
    for (t in 1:n_dates) {
      log_lambda <- beta0 + 
        beta_temp * z_temp[g, t] +
        beta_rain * z_rain[g, t] +
        # beta_buildings * z_buildings[g] +
        sum(beta_land * z_land[g, ])
      lambda_grid[g, t] <- exp(log_lambda)
    }
  }
  
  # Generate survey observations
  site_to_grid <- sample(1:n_grids, n_obs_y, replace = TRUE)
  date_y <- sample(1:n_dates, n_obs_y, replace = TRUE)
  trap_type <- sample(1:5, n_obs_y, replace = TRUE)
  
  # Detection probability parameters
  alpha0 <- qlogis(0.10)
  alpha_RH <- 0.3
  alpha_WS_rain <- -0.2
  alpha_trap <- c(0.1, -0.1, 0.2, -0.2, 0)  # Sum to zero
  
  z_RH <- rnorm(n_obs_y)
  z_WS_rain <- rnorm(n_obs_y)
  
  # Generate observations
  y <- integer(n_obs_y)
  for (k in 1:n_obs_y) {
    g <- site_to_grid[k]
    t <- date_y[k]
    
    # True abundance
    N_true <- rnbinom(1, mu = lambda_grid[g, t], size = phi)
    
    # Detection probability
    p_trap <- plogis(alpha0 + 
                       alpha_RH * z_RH[k] + 
                       alpha_WS_rain * z_WS_rain[k] + 
                       alpha_trap[trap_type[k]])
    
    # Observed count
    y[k] <- rbinom(1, N_true, p_trap)
  }
  
  # Generate presence-only observations
  # Thinning parameters
  delta0 <- -1
  delta_poi <- 0.3
  delta_reports <- 0.2
  
  # Thinning probabilities
  p_thin <- plogis(delta0 + delta_poi * z_poi + delta_reports * z_reports)
  
  # Generate PO data
  po_grid_idx <- integer(n_obs_po)
  date_po <- integer(n_obs_po)
  
  for (r in 1:n_obs_po) {
    # Sample proportional to thinned intensity
    lambda_thinned <- lambda_grid * (p_thin %*% t(rep(1, n_dates)))
    probs <- c(lambda_thinned)
    probs <- probs / sum(probs)
    
    idx <- sample(1:(n_grids * n_dates), 1, prob = probs)
    po_grid_idx[r] <- ((idx - 1) %% n_grids) + 1
    date_po[r] <- ((idx - 1) %/% n_grids) + 1
  }
  
  # Create Stan data list
  stan_data <- list(
    n_grids_total = n_grids,
    n_grids_obs = n_grids,
    n_dates = n_dates,
    n_obs_y = n_obs_y,
    n_obs_po = n_obs_po,
    n_land_covs = n_land_covs,
    
    obs_grid_idx = 1:n_grids,
    site_to_grid = site_to_grid,
    date_y = date_y,
    trap_type = trap_type,
    po_grid_idx = po_grid_idx,
    date_po = date_po,
    
    y = y,
    ones = rep(1L, n_obs_po),
    
    z_temp = z_temp,
    z_rain = z_rain,
    z_ndvi = matrix(0, n_grids, n_dates),  # Not used
    # z_buildings = z_buildings,
    z_land = z_land,
    z_poi = z_poi,
    z_reports = z_reports,
    z_RH = z_RH,
    z_WS_rain = z_WS_rain,
    
    area_grid = area_grid,
    CONSTANT = 1000,
    N_multiplier = 5
  )
  
  # True parameters for comparison
  true_params <- list(
    beta0 = beta0,
    beta_temp = beta_temp,
    beta_rain = beta_rain,
    # beta_buildings = beta_buildings,
    beta_land = beta_land,
    phi = phi,
    alpha0 = alpha0,
    alpha_RH = alpha_RH,
    alpha_WS_rain = alpha_WS_rain,
    delta0 = delta0,
    delta_poi = delta_poi,
    delta_reports = delta_reports
  )
  
  cat("✓ Test data simulated\n")
  cat(sprintf("  %d grids, %d dates\n", n_grids, n_dates))
  cat(sprintf("  %d survey obs, %d PO records\n", n_obs_y, n_obs_po))
  cat(sprintf("  y range: [%d, %d], mean = %.1f\n", 
              min(y), max(y), mean(y)))
  
  return(list(data = stan_data, truth = true_params))
}

# ============================================================================
# 2. TEST THE MARGINALIZATION FUNCTION
# ============================================================================

test_marginalization_function <- function() {
  cat("\nTesting marginalization function...\n")
  
  # Create standalone Stan program to test the function
  test_code <- "
functions {
  real binomial_negbin_marginal_lpmf(int y_obs, real p_detect, 
                                     real lambda, real phi, int N_max) {
    vector[N_max - y_obs + 1] log_components;
    
    for (n in 1:(N_max - y_obs + 1)) {
      int N = y_obs + n - 1;
      log_components[n] = binomial_lpmf(y_obs | N, p_detect) + 
                          neg_binomial_2_lpmf(N | lambda, phi);
    }
    
    return log_sum_exp(log_components);
  }
}

data {
  int y_test;
  real p_test;
  real lambda_test;
  real phi_test;
  int N_max_test;
}

generated quantities {
  real log_lik;
  // Can't use | syntax in generated quantities, call directly
  log_lik = binomial_negbin_marginal_lpmf(y_test| p_test, 
                                           lambda_test, phi_test, N_max_test);
}
"

# Write and compile
test_file <- tempfile(fileext = ".stan")
writeLines(test_code, test_file)

tryCatch({
  mod <- cmdstan_model(test_file)
  
  # Test with simple values
  test_data <- list(
    y_test = 5L,
    p_test = 0.3,
    lambda_test = 10.0,
    phi_test = 2.0,
    N_max_test = 50L
  )
  
  fit <- mod$sample(
    data = test_data,
    chains = 1,
    iter_warmup = 100,
    iter_sampling = 100,
    refresh = 0
  )
  
  log_lik <- fit$summary("log_lik")
  cat("✓ Function test passed\n")
  cat(sprintf("  log_lik mean: %.4f\n", log_lik$mean))
  
  unlink(test_file)
  return(TRUE)
  
}, error = function(e) {
  cat("✗ Function test failed:\n")
  cat(paste0("  ", e$message, "\n"))
  return(FALSE)
})
}

# ============================================================================
# 3. QUICK MODEL TEST WITH SMALL DATA
# ============================================================================

quick_model_test <- function(stan_file, 
                             n_grids = 5, 
                             n_dates = 3,
                             n_obs_y = 20,
                             chains = 2,
                             iter_warmup = 200,
                             iter_sampling = 200) {
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("QUICK MODEL TEST\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  # Simulate test data
  sim <- simulate_test_data(
    n_grids = n_grids,
    n_dates = n_dates,
    n_obs_y = n_obs_y,
    n_obs_po = n_obs_y * 2
  )
  
  # Compile model
  cat("\nCompiling Stan model...\n")
  mod <- cmdstan_model(stan_file)
  
  # Fit model
  cat("\nFitting model...\n")
  fit <- mod$sample(
    data = sim$data,
    chains = chains,
    parallel_chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = max(1, iter_sampling %/% 4),
    show_messages = FALSE
  )
  
  # Check diagnostics
  cat("\nDiagnostics:\n")
  diag <- fit$diagnostic_summary()
  print(diag)
  
  # Compare estimates to truth
  cat("\nParameter Recovery:\n")
  params_to_check <- c("beta0", "beta_temp", "beta_rain", "phi", 
                       "alpha0", "delta0")
  
  results <- fit$summary(params_to_check)
  results$true_value <- sapply(params_to_check, function(p) sim$truth[[p]])
  results$bias <- results$mean - results$true_value
  results$coverage <- (results$q5 <= results$true_value) & 
    (results$true_value <= results$q95)
  
  print(results[, c("variable", "mean", "true_value", "bias", "rhat", "coverage")])
  
  cat("\n✓ Quick test complete\n")
  
  return(list(fit = fit, sim = sim, results = results))
}

# ============================================================================
# 4. DIAGNOSE COMMON ISSUES
# ============================================================================

diagnose_stan_data <- function(stan_data) {
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("DIAGNOSING STAN DATA\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  issues <- character(0)
  
  # Check for NAs
  cat("Checking for NAs...\n")
  na_vars <- names(stan_data)[sapply(stan_data, function(x) any(is.na(x)))]
  if (length(na_vars) > 0) {
    issues <- c(issues, sprintf("NAs found in: %s", paste(na_vars, collapse = ", ")))
    cat("  ✗ ", issues[length(issues)], "\n", sep = "")
  } else {
    cat("  ✓ No NAs\n")
  }
  
  # Check index bounds
  cat("\nChecking index bounds...\n")
  if (any(stan_data$site_to_grid < 1 | stan_data$site_to_grid > stan_data$n_grids_obs)) {
    issues <- c(issues, "site_to_grid has out-of-bounds indices")
    cat("  ✗ ", issues[length(issues)], "\n", sep = "")
  } else {
    cat("  ✓ site_to_grid OK\n")
  }
  
  if (any(stan_data$date_y < 1 | stan_data$date_y > stan_data$n_dates)) {
    issues <- c(issues, "date_y has out-of-bounds indices")
    cat("  ✗ ", issues[length(issues)], "\n", sep = "")
  } else {
    cat("  ✓ date_y OK\n")
  }
  
  # Check dimensions
  cat("\nChecking dimensions...\n")
  if (length(stan_data$y) != stan_data$n_obs_y) {
    issues <- c(issues, sprintf("y length (%d) != n_obs_y (%d)", 
                                length(stan_data$y), stan_data$n_obs_y))
    cat("  ✗ ", issues[length(issues)], "\n", sep = "")
  } else {
    cat("  ✓ y length matches n_obs_y\n")
  }
  
  if (nrow(stan_data$z_temp) != stan_data$n_grids_total) {
    issues <- c(issues, "z_temp rows != n_grids_total")
    cat("  ✗ ", issues[length(issues)], "\n", sep = "")
  } else {
    cat("  ✓ z_temp dimensions OK\n")
  }
  
  # Check data types
  cat("\nChecking data types...\n")
  if (!is.integer(stan_data$y)) {
    issues <- c(issues, "y is not integer")
    cat("  ✗ ", issues[length(issues)], "\n", sep = "")
  } else {
    cat("  ✓ y is integer\n")
  }
  
  # Check for extreme values
  cat("\nChecking for extreme values...\n")
  if (max(stan_data$y) > 1000) {
    cat("  ⚠ Large y values detected (max = ", max(stan_data$y), ")\n", sep = "")
    cat("    Consider increasing N_multiplier\n")
  }
  
  # Summary
  cat("\n", rep("=", 70), "\n", sep = "")
  if (length(issues) == 0) {
    cat("✓ No issues detected\n")
  } else {
    cat("✗ Found ", length(issues), " issue(s):\n", sep = "")
    for (i in seq_along(issues)) {
      cat("  ", i, ". ", issues[i], "\n", sep = "")
    }
  }
  cat(rep("=", 70), "\n", sep = "")
  
  invisible(issues)
}

# ============================================================================
# 5. EXAMPLE USAGE
# ============================================================================

if (TRUE) {  # Set to TRUE to run examples
  
  # Test 1: Test the marginalization function
  test_marginalization_function()
  
  # Test 2: Simulate and test with small data
  test_result <- quick_model_test(
    stan_file = "02_model/stan/marginalised_model.stan",
    n_grids = 5,
    n_dates = 3,
    n_obs_y = 30
  )
  
  # Test 3: Diagnose your real data
  stan_data <- readRDS("01_data/processedCovaraites/stan_data_marginalized_1km.rds")
  diagnose_stan_data(stan_data)
  
}
