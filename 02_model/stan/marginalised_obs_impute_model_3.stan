//This is the full integrated species distribution model for Cx. pipiens.
// The survey data is modelled through an N-mixture model for the count data,
// the citizen science reports is modelled through a Bernoulli poisson point
// process using the ones trick to account for the background reporting.

// This is using the marginlisation of the latent discrete parameter N.
// Must sum over all possible values, but since it ranges from 0 to infinity, 
// there must be a cap/max...

// In the transformed data block: the N_multiplier defined in stan data list
// (prepare_stan_data.r)

// https://mc-stan.org/docs/stan-users-guide/latent-discrete.html
// https://mbjoseph.github.io/posts/2020-04-28-a-step-by-step-guide-to-marginalizing-over-discrete-parameters-for-ecologists-using-stan/
  
 
// Optimized integrated species distribution model for Cx. pipiens with imputation
// Improvements: 
// 1. Vectorized operations where possible
// 2. Adaptive N_max with early stopping in marginalization
// 3. Optimized lambda computations

functions {
  // Marginalized likelihood with EARLY STOPPING and ADAPTIVE N_max
  real binomial_negbin_marginal_lpmf(int y_obs, real p_detect, 
                                     real lambda, real phi, int N_base) {
    // Adaptive N_max based on observation
    int N_max = y_obs + N_base;
    vector[N_max - y_obs + 1] log_components;
    real threshold = log(1e-12);  // Stop when contribution becomes negligible
    int n_components = N_max - y_obs + 1;
    
    for (n in 1:n_components) {
      int N = y_obs + n - 1;
      log_components[n] = binomial_lpmf(y_obs | N, p_detect) + 
                          neg_binomial_2_lpmf(N | lambda, phi);
      
      // Early stopping: if we're past observation and contribution is tiny, stop
      if (n > 20 && log_components[n] < threshold) {
        return log_sum_exp(log_components[1:n]);
      }
    }
    
    return log_sum_exp(log_components);
  }
}

data {
  // Dimensions
  int<lower=1> n_grids_total;
  int<lower=1> n_grids_obs;
  int<lower=1> n_dates;
  int<lower=1> n_obs_y;
  int<lower=1> n_obs_po;
  int<lower=1> n_land_covs;
  
  // Indices  
  array[n_grids_obs] int<lower=1, upper=n_grids_total> obs_grid_idx;
  array[n_obs_y] int<lower=1, upper=n_grids_obs> site_to_grid;
  array[n_obs_y] int<lower=1, upper=n_dates> date_y;
  array[n_obs_y] int<lower=1, upper=5> trap_type;
  array[n_obs_po] int<lower=1, upper=n_grids_total> po_grid_idx;
  array[n_obs_po] int<lower=1, upper=n_dates> date_po;
  
  // Observations
  array[n_obs_y] int<lower=0> y;
  array[n_obs_po] int<lower=1> ones;
  
  // Covariates - WITH MISSING DATA INDICATORS
  matrix[n_grids_total, n_dates] z_temp_obs;
  matrix[n_grids_total, n_dates] z_rain_obs;
  array[n_grids_total, n_dates] int<lower=0, upper=1> temp_missing;
  array[n_grids_total, n_dates] int<lower=0, upper=1> rain_missing;
  int<lower=0> n_temp_missing;
  int<lower=0> n_rain_missing;
  
  matrix[n_grids_total, n_land_covs] z_land;
  vector[n_grids_total] z_poi;
  vector[n_grids_total] z_reports;
  vector[n_obs_y] z_RH;
  vector[n_obs_y] z_WS_rain;
  
  // Other
  vector[n_grids_total] area_grid;
  real CONSTANT;
  int<lower=1> N_multiplier;
}

transformed data {
  int N_base = N_multiplier * 10 + 50;  // Base for adaptive N_max
}

parameters {
  // Abundance model
  real beta0;
  real beta_temp;
  real beta_rain;
  vector[n_land_covs] beta_land;
  
  // Detection model
  real alpha0;
  real alpha_RH;
  real alpha_WS_rain;
  vector[4] alpha_trap_raw;
  real log_sigma_trap;
  
  // Thinning model
  real delta0;
  real delta_poi;
  real delta_reports;
  
  // Dispersion
  real<lower=1e-12> phi;
  
  // IMPUTATION MODEL PARAMETERS
  real gamma_temp0;
  vector[n_land_covs] gamma_temp_land;
  vector[n_dates] gamma_temp_time_raw;
  real<lower=1e-6> sigma_temp;
  real<lower=1e-6> sigma_temp_time;
  
  real gamma_rain0;
  vector[n_land_covs] gamma_rain_land;
  vector[n_dates] gamma_rain_time_raw;
  real<lower=1e-6> sigma_rain;
  real<lower=1e-6> sigma_rain_time;
  
  // IMPUTED VALUES
  vector[n_temp_missing] z_temp_imputed_raw;
  vector[n_rain_missing] z_rain_imputed_raw;
}

transformed parameters {
  // Complete temperature and rain matrices
  matrix[n_grids_total, n_dates] z_temp;
  matrix[n_grids_total, n_dates] z_rain;
  
  // Time effects (sum to zero constraint)
  vector[n_dates] gamma_temp_time = gamma_temp_time_raw - mean(gamma_temp_time_raw);
  vector[n_dates] gamma_rain_time = gamma_rain_time_raw - mean(gamma_rain_time_raw);
  
  real<lower=0> sigma_trap = exp(log_sigma_trap);
  vector[5] alpha_trap;
  
  // Trap effects sum to zero
  alpha_trap[1:4] = alpha_trap_raw;
  alpha_trap[5] = -sum(alpha_trap_raw);
  
  // FILL IN COMPLETE MATRICES (observed + imputed)
  {
    int temp_idx = 1;
    int rain_idx = 1;
    
    for (g in 1:n_grids_total) {
      for (t in 1:n_dates) {
        // Temperature
        if (temp_missing[g, t] == 1) {
          z_temp[g, t] = z_temp_imputed_raw[temp_idx];
          temp_idx += 1;
        } else {
          z_temp[g, t] = z_temp_obs[g, t];
        }
        
        // Rain
        if (rain_missing[g, t] == 1) {
          z_rain[g, t] = z_rain_imputed_raw[rain_idx];
          rain_idx += 1;
        } else {
          z_rain[g, t] = z_rain_obs[g, t];
        }
      }
    }
  }
  
  // Pre-compute land effects for all grids
  vector[n_grids_total] land_effect = z_land * beta_land;
  
  // Lambda computations
  matrix[n_grids_obs, n_dates] lambda_obs;
  matrix[n_grids_total, n_dates] lambda_thinned;
  vector[n_dates] background;
  
  // OPTIMIZED: Vectorized lambda computation for observed grids
  for (i in 1:n_grids_obs) {
    int g = obs_grid_idx[i];
    // Vectorize across dates
    vector[n_dates] log_lambda_raw = beta0 + 
                                      beta_temp * to_vector(z_temp[g, ]) +
                                      beta_rain * to_vector(z_rain[g, ]) +
                                      rep_vector(land_effect[g], n_dates);
    
    lambda_obs[i, ] = to_row_vector(exp(fmin(fmax(log_lambda_raw, -10), 10)));
  }
  
  // OPTIMIZED: Vectorized lambda_thinned for ALL grids
  {
    // Pre-compute thinning probabilities
    vector[n_grids_total] p_thin = inv_logit(delta0 + 
                                              delta_poi * z_poi + 
                                              delta_reports * z_reports);
    
    for (g in 1:n_grids_total) {
      // Vectorize across dates
      vector[n_dates] log_lambda_raw = beta0 + 
                                        beta_temp * to_vector(z_temp[g, ]) +
                                        beta_rain * to_vector(z_rain[g, ]) +
                                        rep_vector(land_effect[g], n_dates);
      
      vector[n_dates] lambda_g = exp(fmin(fmax(log_lambda_raw, -10), 10));
      lambda_thinned[g, ] = to_row_vector(lambda_g * area_grid[g] * p_thin[g]);
    }
  }
  
  // Background for PO likelihood
  for (t in 1:n_dates) {
    background[t] = sum(lambda_thinned[, t]) / n_obs_po;
  }
}

model {
  // PRIORS - Abundance model
  beta0 ~ normal(0, 1);
  beta_temp ~ normal(0, 1);
  beta_rain ~ normal(0, 1);
  beta_land ~ normal(0, 1);
  
  // PRIORS - Detection model
  alpha0 ~ normal(logit(0.10), 0.5);
  alpha_RH ~ normal(0, 1);
  alpha_WS_rain ~ normal(0, 1);
  alpha_trap_raw ~ normal(0, sigma_trap);
  log_sigma_trap ~ normal(log(0.5), 0.5); 
  
  // PRIORS - Thinning model
  delta0 ~ normal(-1, 1);
  delta_poi ~ normal(0, 1);
  delta_reports ~ normal(0, 1);
  
  // PRIORS - Dispersion
  phi ~ normal(2, 1);
  
  // PRIORS - Imputation model
  gamma_temp0 ~ normal(0, 1);
  gamma_temp_land ~ normal(0, 1);
  gamma_temp_time_raw ~ normal(0, sigma_temp_time);
  sigma_temp ~ exponential(1);
  sigma_temp_time ~ exponential(1);
  
  gamma_rain0 ~ normal(0, 1);
  gamma_rain_land ~ normal(0, 1);
  gamma_rain_time_raw ~ normal(0, sigma_rain_time);
  sigma_rain ~ exponential(1);
  sigma_rain_time ~ exponential(1);
  
  // IMPUTATION LIKELIHOOD
  // Pre-compute imputation means for efficiency
  matrix[n_grids_total, n_dates] mu_temp;
  matrix[n_grids_total, n_dates] mu_rain;
  
  for (g in 1:n_grids_total) {
    vector[n_dates] temp_base = gamma_temp0 + 
                                 rep_vector(dot_product(gamma_temp_land, z_land[g, ]), n_dates) +
                                 gamma_temp_time;
    vector[n_dates] rain_base = gamma_rain0 + 
                                 rep_vector(dot_product(gamma_rain_land, z_land[g, ]), n_dates) +
                                 gamma_rain_time;
    
    mu_temp[g, ] = to_row_vector(temp_base);
    mu_rain[g, ] = to_row_vector(rain_base);
  }
  
  // Temperature: observed values inform the imputation model
  for (g in 1:n_grids_total) {
    for (t in 1:n_dates) {
      if (temp_missing[g, t] == 0) {
        z_temp_obs[g, t] ~ normal(mu_temp[g, t], sigma_temp);
      }
    }
  }
  
  // Rain: observed values inform the imputation model
  for (g in 1:n_grids_total) {
    for (t in 1:n_dates) {
      if (rain_missing[g, t] == 0) {
        z_rain_obs[g, t] ~ normal(mu_rain[g, t], sigma_rain);
      }
    }
  }
  
  // Prior on imputed values using the imputation model
  {
    int temp_idx = 1;
    int rain_idx = 1;
    
    for (g in 1:n_grids_total) {
      for (t in 1:n_dates) {
        if (temp_missing[g, t] == 1) {
          z_temp_imputed_raw[temp_idx] ~ normal(mu_temp[g, t], sigma_temp);
          temp_idx += 1;
        }
        
        if (rain_missing[g, t] == 1) {
          z_rain_imputed_raw[rain_idx] ~ normal(mu_rain[g, t], sigma_rain);
          rain_idx += 1;
        }
      }
    }
  }
  
  // OPTIMIZED SURVEY LIKELIHOOD with vectorized detection probability
  {
    // Pre-compute detection probabilities
    vector[n_obs_y] lambda_vec;
    vector[n_obs_y] p_trap_vec;
    
    for (k in 1:n_obs_y) {
      int i = site_to_grid[k];
      int t = date_y[k];
      lambda_vec[k] = lambda_obs[i, t];
      p_trap_vec[k] = inv_logit(alpha0 + 
                                 alpha_RH * z_RH[k] + 
                                 alpha_WS_rain * z_WS_rain[k] + 
                                 alpha_trap[trap_type[k]]);
    }
    
    // Marginalized likelihood with optimized function
    for (k in 1:n_obs_y) {
      target += binomial_negbin_marginal_lpmf(y[k] | p_trap_vec[k], 
                                               lambda_vec[k], phi, N_base);
    }
  }
  
  // PRESENCE-ONLY LIKELIHOOD
  for (r in 1:n_obs_po) {
    int g = po_grid_idx[r];
    int t = date_po[r];
    target += log(lambda_thinned[g, t]) - log(background[t]) - log(CONSTANT);
  }
}

generated quantities {
  // -------------------------------------------------
  // EXPECTED ABUNDANCE
  // -------------------------------------------------
  array[n_grids_obs, n_dates] real N_expected;
  
  for (i in 1:n_grids_obs) {
    for (t in 1:n_dates) {
      N_expected[i, t] = lambda_obs[i, t];
    }
  }
  
  // -------------------------------------------------
  // POSTERIOR PREDICTIVE REPLICATIONS
  // -------------------------------------------------
  array[n_obs_y] int y_rep;
  array[n_obs_po] int po_rep;
  
  // -------------------------------------------------
  // PPC SUMMARY STATISTICS (OBSERVED)
  // -------------------------------------------------
  real mean_y_obs = mean(to_vector(y));
  real sd_y_obs   = sd(to_vector(y));
  int  max_y_obs  = max(y);
  
  real prop_zeros_obs;
  {
    int zeros_obs = 0;
    for (k in 1:n_obs_y) {
      if (y[k] == 0) zeros_obs += 1;
    }
    prop_zeros_obs = zeros_obs * 1.0 / n_obs_y;
  }
  
  // -------------------------------------------------
  // PPC SUMMARY STATISTICS (REPLICATED)
  // -------------------------------------------------
  real mean_y_rep;
  real sd_y_rep;
  int  max_y_rep;
  real prop_zeros_rep;
  
  {
    int zeros_rep = 0;
    vector[n_obs_y] y_rep_vec;
    
    for (k in 1:n_obs_y) {
      int i = site_to_grid[k];
      int t = date_y[k];
      
      real p_trap = inv_logit(
        alpha0 +
        alpha_RH * z_RH[k] +
        alpha_WS_rain * z_WS_rain[k] +
        alpha_trap[trap_type[k]]
      );
      
      int N_sim = neg_binomial_2_rng(lambda_obs[i, t], phi);
      if (N_sim > N_base * 2) N_sim = N_base * 2;  // Use N_base instead of N_max
      
      y_rep[k] = binomial_rng(N_sim, p_trap);
      y_rep_vec[k] = y_rep[k];
      
      if (y_rep[k] == 0) zeros_rep += 1;
    }
    
    mean_y_rep  = mean(y_rep_vec);
    sd_y_rep    = sd(y_rep_vec);
    max_y_rep   = max(y_rep);
    prop_zeros_rep = zeros_rep * 1.0 / n_obs_y;
  }
  
  // -------------------------------------------------
  // PRESENCE-ONLY PPC (BINARY)
  // -------------------------------------------------
  for (r in 1:n_obs_po) {
    int g = po_grid_idx[r];
    int t = date_po[r];
    
    real p_thin = inv_logit(
      delta0 +
      delta_poi * z_poi[g] +
      delta_reports * z_reports[g]
    );
    
    // Use pre-computed land effect for efficiency
    real lambda_g = exp(
      fmin(fmax(
        beta0 +
        beta_temp * z_temp[g, t] +
        beta_rain * z_rain[g, t] +
        land_effect[g],
        -10), 10)
    );
    
    int N_sim = neg_binomial_2_rng(lambda_g, phi);
    int N_thinned = binomial_rng(N_sim, p_thin);
    po_rep[r] = N_thinned > 0 ? 1 : 0;
  }
}

