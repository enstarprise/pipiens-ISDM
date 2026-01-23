
library(nimble)

model_code <- nimbleCode({
  ## ============================================================
  ## PRIORS
  ## ============================================================
  beta0 ~ dnorm(0, sd = 0.5)
  beta_temp ~ dnorm(0, sd = 1)
  beta_rain ~ dnorm(0, sd = 1)
  
  for(k in 1:n_land_covs) {
    beta_land[k] ~ dnorm(0, sd = 1)
  }
  
  # Baseline detection from MRR literature
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
  delta0 ~ dnorm(0, sd = 0.5)
  delta_poi ~ dnorm(0, sd = 1)
  delta_reports ~ dnorm(0, sd = 1)
  
  ## Dispersion parameters
  phi ~ dgamma(3, 1.5) 
  
  ## ============================================================
  ## SPATIAL PARAMETERS (Matérn Exponential)
  ## ============================================================
  # Spatial variance (marginal variance of the spatial field)
  sigma_spatial ~ dgamma(2, 1)     # marginal SD
  
  # Range parameter (distance at which correlation drops to ~0.05)
  # Prior should reflect your domain knowledge about spatial scale
  range_spatial ~ dgamma(5, 1)  # Adjust shape/rate based on your grid spacing
  
  spatial_effect[1:n_grids_total] ~ dmnorm(
    mean = rep(0, n_grids_total),
    cov = Sigma[1:n_grids_total, 1:n_grids_total]
  )
  
  ## ============================================================
  ## DERIVED VARIABLES
  ## ============================================================
  # Constrain trap effects to sum to zero
  for(i in 1:4) {
    alpha_trap[i] <- alpha_trap_raw[i]
  }
  alpha_trap[5] <- -sum(alpha_trap_raw[1:4])
  
  for(i in 1:n_grids_total) {
    for(j in 1:n_grids_total) {
      Sigma[i, j] <- pow(sigma_spatial, 2) *
        exp(-D[i, j] / rho)
    }
  }
  
  
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
        inprod(beta_land[1:n_land_covs], z_land[g, 1:n_land_covs])  +
        spatial_effect[g]  # ADDED SPATIAL RANDOM EFFECT
      
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
    logit(p_trap[k]) <- alpha0 +     # informed baseline from MRR
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
