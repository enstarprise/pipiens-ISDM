library(tidyverse)

# ===============================================================================
# DIAGNOSTIC FUNCTIONS FOR EXPONENTIAL KERNEL (Matérn ν=0.5)
# ===============================================================================

# For exponential kernel, the diagnostic equations are simpler than squared exponential
# Based on Riutort-Mayol et al. (2022) equations for Matérn family

# Minimum length-scale that can be accurately inferred
# For exponential kernel (Matérn ν=0.5):
min_rho_exponential <- function(m, c, L) {
  # m: number of basis functions (one dimension, so sqrt(M_total))
  # c: boundary factor
  # L: boundary (domain size)
  
  # for Matérn ν=0.5, the relationship is approximately:
  # rho_min ≈ (c * L) / (2 * m)
  # to ensure the highest frequency basis function can capture the correlation
  
  (c * L) / (2 * m)
}

# Minimum boundary factor given rho
min_c_exponential <- function(rho, L) {
  # c must be large enough so that:
  # exp(-c) << 1 (negligible boundary effects)
  # For exponential kernel: c ≥ 3 is typically sufficient
  # But also must satisfy: c * L / rho ≥ some minimum
  
  # Conservative choice: ensure correlation at boundary is < 0.01
  # exp(-c*L/rho) < 0.01
  # c*L/rho > 4.6
  max(3, 4.6 * rho / L)
}

# Minimum number of basis functions given rho and c
min_m_exponential <- function(rho, c, L) {
  # From inverting min_rho equation:
  # m ≥ (c * L) / (2 * rho)
  
  ceiling((c * L) / (2 * rho))
}

# ===============================================================================
# DIAGNOSTIC PROCEDURE
# ===============================================================================

run_hsgp_diagnostics <- function(fit, grid_coords_km, M_current, 
                                 c_current = NULL, max_iterations = 5) {
  
  # Compute boundary (domain size)
  L_x <- diff(range(grid_coords_km[, 1]))
  L_y <- diff(range(grid_coords_km[, 2]))
  L <- max(L_x, L_y)  # Use maximum dimension
  
  # Get current M (total basis functions)
  # For 2D: M = m_x * m_y, typically m_x = m_y
  m <- sqrt(M_current)
  
  if (is.null(c_current)) {
    c_current <- 1.5  # default
  }
  
  cat("\n=== HSGP DIAGNOSTIC PROCEDURE ===\n")
  cat("Domain size L:", round(L, 2), "km\n")
  cat("Current M:", M_current, "(m =", m, "per dimension)\n")
  cat("Current c:", c_current, "\n\n")
  
  # ===============================================================================
  # PHASE A: Length-scale diagnostic
  # ===============================================================================
  
  cat("=== PHASE A: Length-scale Diagnostic ===\n\n")
  
  results <- list()
  
  for (k in 1:max_iterations) {
    cat("Iteration", k, ":\n")
    
    # Extract estimated rho from current model
    rho_draws <- fit$draws("rho_spatial", format = "draws_matrix")
    rho_hat <- mean(rho_draws)
    rho_sd <- sd(rho_draws)
    
    cat("  Estimated rho:", round(rho_hat, 3), "±", round(rho_sd, 3), "km\n")
    
    # Compute minimum rho that can be accurately inferred
    rho_min <- min_rho_exponential(m, c_current, L)
    cat("  Minimum rho (given m, c):", round(rho_min, 3), "km\n")
    
    # Diagnostic check: rho_hat + 0.01 >= rho_min?
    diagnostic_pass <- (rho_hat + 0.01) >= rho_min
    cat("  Diagnostic:", ifelse(diagnostic_pass, "PASS ✓", "FAIL ✗"), "\n")
    
    # Compute recommended m and c for current rho estimate
    c_recommended <- min_c_exponential(rho_hat, L)
    m_recommended <- min_m_exponential(rho_hat, c_recommended, L)
    M_recommended <- m_recommended^2  # For 2D
    
    cat("  Recommended c:", round(c_recommended, 3), "\n")
    cat("  Recommended m:", m_recommended, "(M =", M_recommended, ")\n")
    
    # Store results
    results[[k]] <- data.frame(
      iteration = k,
      rho_hat = rho_hat,
      rho_sd = rho_sd,
      rho_min = rho_min,
      m_current = m,
      M_current = M_current,
      c_current = c_current,
      diagnostic_pass = diagnostic_pass,
      m_recommended = m_recommended,
      M_recommended = M_recommended,
      c_recommended = c_recommended
    )
    
    if (diagnostic_pass) {
      cat("\n✓ Phase A diagnostic passed! Proceeding to Phase B...\n\n")
      break
    } else {
      cat("\n✗ Approximation may be inaccurate. Need to increase m and/or c.\n")
      cat("  Suggested: M =", M_recommended, ", c =", round(c_recommended, 2), "\n\n")
      
      if (k == max_iterations) {
        cat("⚠ Maximum iterations reached. Model needs refit with larger M and c.\n\n")
      }
    }
  }
  
  # ===============================================================================
  # PHASE B: Stability check (if diagnostic passed)
  # ===============================================================================
  
  if (diagnostic_pass) {
    cat("=== PHASE B: Stability Check ===\n")
    cat("Current approximation appears adequate.\n")
    cat("To verify stability, you should:\n")
    cat("1. Increase M by ~20% (e.g., M =", ceiling(M_current * 1.2), ")\n")
    cat("2. Refit model and check if rho_hat remains stable\n")
    cat("3. Compare ELPD/LOO scores\n\n")
  }
  
  # ===============================================================================
  # Summary
  # ===============================================================================
  
  cat("\n=== DIAGNOSTIC SUMMARY ===\n")
  
  results_df <- bind_rows(results)
  print(results_df)
  
  # Final recommendation
  cat("\n=== FINAL RECOMMENDATION ===\n")
  if (diagnostic_pass) {
    cat("✓ Current approximation (M =", M_current, ", c =", c_current, ") is adequate\n")
    cat("  Estimated range:", round(rho_hat, 2), "km\n")
    cat("  Effective range (3×rho):", round(3 * rho_hat, 2), "km\n")
  } else {
    final <- tail(results_df, 1)
    cat("✗ Current approximation is INSUFFICIENT\n")
    cat("  Refit model with:\n")
    cat("    M ≥", final$M_recommended, "(m ≥", final$m_recommended, "per dimension)\n")
    cat("    c ≥", round(final$c_recommended, 2), "\n")
  }
  
  return(invisible(results_df))
}

# ===============================================================================
# LOAD YOUR MODEL AND DATA
# ===============================================================================

# Assuming you have:
# - fit: your Stan fit object
# - grid_coords: your BNG coordinates

# Convert to km and center (as you did for fitting)
grid_coords_km <- grid_coords / 1000
grid_coords_centered <- scale(grid_coords_km, center = TRUE, scale = FALSE)

# Run diagnostics
diagnostic_results <- run_hsgp_diagnostics(
  fit = fit,
  grid_coords_km = grid_coords_centered,
  M_current = 64,  # Your current M (e.g., 8×8)
  c_current = 1.5,  # Your current boundary factor
  max_iterations = 3
)

# Save results
write_csv(diagnostic_results, "hsgp_diagnostic_results.csv")





library(ggplot2)

# ===============================================================================
# PLOT 1: Estimated vs Minimum Rho
# ===============================================================================

ggplot(diagnostic_results, aes(x = iteration)) +
  geom_line(aes(y = rho_hat, color = "Estimated ρ"), linewidth = 1) +
  geom_ribbon(aes(ymin = rho_hat - rho_sd, 
                  ymax = rho_hat + rho_sd, 
                  fill = "Estimated ρ"), 
              alpha = 0.3) +
  geom_line(aes(y = rho_min, color = "Minimum ρ (threshold)"), 
            linewidth = 1, linetype = "dashed") +
  geom_point(aes(y = rho_hat, color = "Estimated ρ"), size = 3) +
  geom_point(aes(y = rho_min, color = "Minimum ρ (threshold)"), size = 3) +
  labs(
    title = "HSGP Diagnostic: Length-Scale Check",
    subtitle = "Estimated ρ should exceed minimum ρ for accurate approximation",
    x = "Diagnostic Iteration",
    y = "Length-Scale ρ (km)",
    color = NULL,
    fill = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("hsgp_diagnostic_rho.png", width = 8, height = 6)

# ===============================================================================
# PLOT 2: Recommended vs Current M
# ===============================================================================

ggplot(diagnostic_results, aes(x = iteration)) +
  geom_col(aes(y = M_current, fill = "Current M"), 
           position = "dodge", alpha = 0.6) +
  geom_col(aes(y = M_recommended, fill = "Recommended M"), 
           position = "dodge", alpha = 0.6) +
  geom_text(aes(y = M_current, label = M_current), 
            vjust = -0.5, size = 4) +
  geom_text(aes(y = M_recommended, label = M_recommended), 
            vjust = -0.5, size = 4) +
  labs(
    title = "HSGP Diagnostic: Basis Functions",
    subtitle = "Current M should meet or exceed recommended M",
    x = "Diagnostic Iteration",
    y = "Number of Basis Functions (M)",
    fill = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("hsgp_diagnostic_M.png", width = 8, height = 6)

# ===============================================================================
# PLOT 3: Diagnostic Pass/Fail
# ===============================================================================

diagnostic_results %>%
  mutate(
    status = ifelse(diagnostic_pass, "PASS", "FAIL"),
    status = factor(status, levels = c("FAIL", "PASS"))
  ) %>%
  ggplot(aes(x = iteration, y = 1, fill = status)) +
  geom_tile(height = 0.5, color = "white", linewidth = 2) +
  geom_text(aes(label = status), size = 6, fontface = "bold") +
  scale_fill_manual(values = c("FAIL" = "#d73027", "PASS" = "#1a9850")) +
  labs(
    title = "HSGP Diagnostic Status by Iteration",
    x = "Iteration",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )

ggsave("hsgp_diagnostic_status.png", width = 8, height = 3)
