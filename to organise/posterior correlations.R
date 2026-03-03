library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(corrplot)

# Assuming you have a fitted model object called 'fit'
# If you need to load it from a file:
# fit <- readRDS("path/to/your/fitted_model.rds")

# Function to check posterior correlations
check_posterior_correlations <- function(fit, 
                                         params = NULL,
                                         threshold = 0.7,
                                         method = "pearson") {
  
  # Extract posterior draws
  draws <- fit$draws(format = "df")
  
  # If no specific parameters specified, use all main parameters
  if (is.null(params)) {
    params <- c("beta0", "beta_temp", "beta_rain", 
                paste0("beta_land[", 1:10, "]"),  # adjust based on n_land_covs
                "alpha0", "alpha_RH", "alpha_WS_rain",
                paste0("alpha_trap[", 1:5, "]"),
                "delta0", "delta_poi", "delta_reports",
                "phi", "sigma_trap")
    
    # Keep only parameters that exist in the draws
    params <- params[params %in% names(draws)]
  }
  
  # Extract parameter draws
  param_draws <- draws[, params]
  
  # Calculate correlation matrix
  cor_matrix <- cor(param_draws, method = method)
  
  # Print high correlations
  cat("\n=== Moderate-High Correlations (|r| >", threshold, ") ===\n\n")
  high_cors <- which(abs(cor_matrix) > threshold & 
                       abs(cor_matrix) < 1, arr.ind = TRUE)
  
  if (nrow(high_cors) > 0) {
    # Remove duplicates (upper triangle)
    high_cors <- high_cors[high_cors[,1] < high_cors[,2], , drop = FALSE]
    
    if (nrow(high_cors) > 0) {
      for (i in 1:nrow(high_cors)) {
        row_idx <- high_cors[i, 1]
        col_idx <- high_cors[i, 2]
        cor_val <- cor_matrix[row_idx, col_idx]
        cat(sprintf("%s <-> %s: %.3f\n", 
                    rownames(cor_matrix)[row_idx],
                    colnames(cor_matrix)[col_idx],
                    cor_val))
      }
    } else {
      cat("No high correlations found.\n")
    }
  } else {
    cat("No high correlations found.\n")
  }
  
  return(cor_matrix)
}

# Function to visualize correlations
plot_posterior_correlations <- function(cor_matrix, 
                                        title = "Posterior Correlations") {
  
  # Create correlation plot
  corrplot(cor_matrix, 
           method = "color",
           type = "upper",
           order = "hclust",
           tl.col = "black",
           tl.srt = 45,
           tl.cex = 0.7,
           title = title,
           mar = c(0,0,2,0))
  
  # Also create a ggplot version for more control
  library(reshape2)
  cor_melted <- melt(cor_matrix)
  
  p <- ggplot(cor_melted, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1,1),
                         name = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank()) +
    coord_fixed() +
    labs(title = title)
  
  return(p)
}

# Function to check pairs plots for highly correlated parameters
plot_posterior_pairs <- function(fit, params, max_pairs = 6) {
  
  if (length(params) > max_pairs) {
    warning(paste("Too many parameters. Showing only first", max_pairs))
    params <- params[1:max_pairs]
  }
  
  mcmc_pairs(fit$draws(params),
             diag_fun = "dens",
             off_diag_fun = "hex")
}

# Example usage:
# ============================================
# 1. Check all correlations
cor_matrix <- check_posterior_correlations(fit,
                                           threshold = 0.5)

# 2. Check specific parameters
# cor_matrix <- check_posterior_correlations(fit, 
#                                           params = c("beta0", "beta_temp", 
#                                                     "alpha0", "delta0"))

# 3. Visualize correlations
plot_posterior_correlations(cor_matrix)

# 4. Create pairs plot for highly correlated parameters
# plot_posterior_pairs(fit, c("beta0", "alpha0", "delta0"))

# ============================================
# Quick diagnostic function
# ============================================
full_correlation_check <- function(fit) {
  
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("POSTERIOR CORRELATION DIAGNOSTICS\n")
  cat(rep("=", 60), "\n\n", sep = "")
  
  # Get all main parameters
  all_params <- c("beta0", "beta_temp", "beta_rain", 
                  "alpha0", "alpha_RH", "alpha_WS_rain",
                  "delta0", "delta_poi", "delta_reports",
                  "phi", "sigma_trap")
  
  # Check which exist
  draws <- fit$draws(format = "df")
  existing_params <- all_params[all_params %in% names(draws)]
  
  # Add array parameters that exist
  beta_land_params <- grep("^beta_land\\[", names(draws), value = TRUE)
  alpha_trap_params <- grep("^alpha_trap\\[", names(draws), value = TRUE)
  
  all_existing <- c(existing_params, beta_land_params, alpha_trap_params)
  
  # Calculate correlations
  cor_matrix <- check_posterior_correlations(fit, params = all_existing)
  
  # Plot
  cat("\n\nGenerating correlation plot...\n")
  p <- plot_posterior_correlations(cor_matrix)
  print(p)
  
  # Save plot
  ggsave("posterior_correlations.png", p, width = 10, height = 8)
  cat("\nPlot saved to: posterior_correlations.png\n")
  
  # Return correlation matrix
  invisible(cor_matrix)
}

# Run the full check:
# cor_matrix <- full_correlation_check(fit)


# ------------
# Quick Posterior Correlation Check
# ===================================

library(posterior)
library(bayesplot)

# Load your fitted model
# fit <- readRDS("your_model.rds")  # or however you have it stored

# Quick check function
quick_cor_check <- function(fit) {
  
  # Extract draws
  draws_df <- as_draws_df(fit$draws())
  
  # Get parameter names (excluding lp__, chain, iteration, draw)
  param_names <- setdiff(names(draws_df), 
                         c("lp__", ".chain", ".iteration", ".draw"))
  
  # Focus on main parameters (not indexed)
  main_params <- grep("\\[", param_names, value = TRUE, invert = TRUE)
  
  # Calculate correlation matrix
  cor_mat <- cor(draws_df[, main_params])
  
  # Print correlations > 0.7
  cat("\nHigh correlations (|r| > 0.7):\n")
  cat("================================\n")
  for (i in 1:(nrow(cor_mat)-1)) {
    for (j in (i+1):ncol(cor_mat)) {
      if (abs(cor_mat[i,j]) > 0.7) {
        cat(sprintf("%15s <-> %-15s: %6.3f\n", 
                    rownames(cor_mat)[i],
                    colnames(cor_mat)[j],
                    cor_mat[i,j]))
      }
    }
  }
  
  return(cor_mat)
}

# Use it:
# cor_mat <- quick_cor_check(fit)

# Or check specific parameters:
# draws_df <- as_draws_df(fit$draws())
# params_of_interest <- c("beta0", "alpha0", "delta0", "phi")
# cor(draws_df[, params_of_interest])

# Pairs plot for visual inspection:
# mcmc_pairs(fit$draws(c("beta0", "alpha0", "delta0")))


##################################################
# CORRELATION: N and P_DETECT


# Quick N vs p_detect Correlation Check
# ======================================

library(posterior)
library(ggplot2)

quick_N_p_correlation <- function(fit, data_list) {
  
  cat("\nExtracting posterior samples...\n")
  draws_df <- as_draws_df(fit$draws())
  
  # Use median draw as representative
  draw <- draws_df[ceiling(nrow(draws_df)/2), ]
  
  n_obs <- data_list$n_obs_y
  N_expected <- numeric(n_obs)
  p_detect <- numeric(n_obs)
  
  cat("Computing E[N] and p_detect for each observation...\n")
  
  for (k in 1:n_obs) {
    i <- data_list$site_to_grid[k]
    t <- data_list$date_y[k]
    
    # Compute lambda (expected abundance)
    log_lambda <- draw$beta0 + 
      draw$beta_temp * data_list$z_temp[i, t] +
      draw$beta_rain * data_list$z_rain[i, t]
    
    # Add land covariates if they exist
    if (data_list$n_land_covs > 0) {
      for (j in 1:data_list$n_land_covs) {
        beta_land_j <- draw[[paste0("beta_land[", j, "]")]]
        log_lambda <- log_lambda + beta_land_j * data_list$z_land[i, j]
      }
    }
    
    N_expected[k] <- exp(log_lambda)
    
    # Compute detection probability
    trap_type_k <- data_list$trap_type[k]
    alpha_trap_k <- draw[[paste0("alpha_trap[", trap_type_k, "]")]]
    
    p_detect[k] <- plogis(draw$alpha0 + 
                            draw$alpha_RH * data_list$z_RH[k] + 
                            draw$alpha_WS_rain * data_list$z_WS_rain[k] + 
                            alpha_trap_k)
  }
  
  # Calculate correlation
  correlation <- cor(N_expected, p_detect)
  
  cat("\n==============================================\n")
  cat("Correlation between E[N] and p_detect:", round(correlation, 4), "\n")
  cat("==============================================\n\n")
  
  # Summary statistics
  cat("E[N] summary:\n")
  print(summary(N_expected))
  cat("\np_detect summary:\n")
  print(summary(p_detect))
  
  # Create scatter plot
  plot_data <- data.frame(
    N = N_expected,
    p = p_detect,
    y_obs = data_list$y
  )
  
  p <- ggplot(plot_data, aes(x = N, y = p)) +
    geom_point(aes(color = y_obs, size = y_obs), alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
    scale_color_viridis_c(name = "Observed\nCount") +
    scale_size_continuous(name = "Observed\nCount", range = c(1, 5)) +
    labs(
      title = sprintf("E[N] vs p_detect (r = %.3f)", correlation),
      x = "Expected Abundance E[N]",
      y = "Detection Probability"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right"
    )
  
  print(p)
  
  # Additional diagnostic: check if correlation varies by trap type
  cat("\n\nCorrelation by trap type:\n")
  for (trap in 1:5) {
    idx <- data_list$trap_type == trap
    if (sum(idx) > 5) {
      cor_trap <- cor(N_expected[idx], p_detect[idx])
      cat(sprintf("  Trap %d: %.4f (n = %d)\n", 
                  trap, cor_trap, sum(idx)))
    }
  }
  
  return(list(
    correlation = correlation,
    N_expected = N_expected,
    p_detect = p_detect,
    plot = p
  ))
}

# Example usage:
# ===============================================
result <- quick_N_p_correlation(fit, stan_data)
# 
# # Save the plot
# ggsave("N_vs_pdetect.png", result$plot, width = 8, height = 6)

# To check across multiple posterior draws:
# ===============================================
check_correlation_stability <- function(fit, data_list, n_draws = 100) {
  
  draws_df <- as_draws_df(fit$draws())
  n_total <- nrow(draws_df)
  
  # Sample draws
  sample_idx <- seq(1, n_total, length.out = min(n_draws, n_total))
  correlations <- numeric(length(sample_idx))
  
  cat(sprintf("Checking correlation across %d posterior draws...\n", 
              length(sample_idx)))
  
  pb <- txtProgressBar(min = 0, max = length(sample_idx), style = 3)
  
  for (idx in seq_along(sample_idx)) {
    draw <- draws_df[sample_idx[idx], ]
    
    N_expected <- numeric(data_list$n_obs_y)
    p_detect <- numeric(data_list$n_obs_y)
    
    for (k in 1:data_list$n_obs_y) {
      i <- data_list$site_to_grid[k]
      t <- data_list$date_y[k]
      
      log_lambda <- draw$beta0 + 
        draw$beta_temp * data_list$z_temp[i, t] +
        draw$beta_rain * data_list$z_rain[i, t]
      
      if (data_list$n_land_covs > 0) {
        for (j in 1:data_list$n_land_covs) {
          beta_land_j <- draw[[paste0("beta_land[", j, "]")]]
          log_lambda <- log_lambda + beta_land_j * data_list$z_land[i, j]
        }
      }
      
      N_expected[k] <- exp(log_lambda)
      
      trap_type_k <- data_list$trap_type[k]
      alpha_trap_k <- draw[[paste0("alpha_trap[", trap_type_k, "]")]]
      
      p_detect[k] <- plogis(draw$alpha0 + 
                              draw$alpha_RH * data_list$z_RH[k] + 
                              draw$alpha_WS_rain * data_list$z_WS_rain[k] + 
                              alpha_trap_k)
    }
    
    correlations[idx] <- cor(N_expected, p_detect)
    setTxtProgressBar(pb, idx)
  }
  close(pb)
  
  cat("\n\n=== Correlation Stability Across Posterior ===\n")
  cat(sprintf("Mean: %.4f\n", mean(correlations)))
  cat(sprintf("SD: %.4f\n", sd(correlations)))
  cat(sprintf("95%% CI: [%.4f, %.4f]\n", 
              quantile(correlations, 0.025),
              quantile(correlations, 0.975)))
  
  # Plot distribution
  p <- ggplot(data.frame(cor = correlations), aes(x = cor)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = mean(correlations), 
               color = "red", linetype = "dashed", linewidth = 1) +
    labs(
      title = "Correlation between E[N] and p_detect Across Posterior",
      subtitle = sprintf("Mean = %.3f, SD = %.3f", 
                         mean(correlations), sd(correlations)),
      x = "Correlation",
      y = "Count"
    ) +
    theme_minimal()
  
  print(p)
  
  return(correlations)
}

# Usage:
cors <- check_correlation_stability(fit, stan_data, n_draws = 100)
