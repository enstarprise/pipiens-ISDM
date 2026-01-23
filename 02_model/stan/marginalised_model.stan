functions {
  // Compute marginalized log-likelihood for one observation
  // Sums over N = y_obs to N_max
  real binomial_negbin_marginal_lpmf(int y_obs, real p_detect, real lambda, real phi, int N_max) {
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
  
  // Covariates
  matrix[n_grids_total, n_dates] z_temp;
  matrix[n_grids_total, n_dates] z_rain;
  // matrix[n_grids_total, n_dates] z_ndvi;
  // vector[n_grids_total] z_buildings;
  matrix[n_grids_total, n_land_covs] z_land;
  vector[n_grids_total] z_poi;
  vector[n_grids_total] z_reports;
  vector[n_obs_y] z_RH;
  vector[n_obs_y] z_WS_rain;
  
  // Other
  vector[n_grids_total] area_grid;
  real CONSTANT;
  
  // Upper bound for marginalization (as multiple of max(y))
  int<lower=1> N_multiplier; // e.g., 5
}

transformed data {
  // Set upper bound for N summation
  int N_max = max(y) * N_multiplier + 100;
}

parameters {
  real beta0;
  real beta_temp;
  real beta_rain;
  // real beta_buildings;
  // real beta_ndvi;
  vector[n_land_covs] beta_land;
  
  real alpha0;
  real alpha_RH;
  real alpha_WS_rain;
  vector[4] alpha_trap_raw;
  real<lower=0> sigma_trap;
  
  real delta0;
  real delta_poi;
  real delta_reports;
  
  real<lower=0> phi;
}

transformed parameters {
  matrix[n_grids_total, n_dates] lambda_grid;
  matrix[n_grids_total, n_dates] lambda_thinned;
  vector[n_dates] background;
  vector[n_grids_total] p_thin;
  vector[5] alpha_trap;
  
  alpha_trap[1:4] = alpha_trap_raw;
  alpha_trap[5] = -sum(alpha_trap_raw);
  
  for (g in 1:n_grids_total) {
    p_thin[g] = inv_logit(delta0 + delta_poi * z_poi[g] + delta_reports * z_reports[g]);
  }
  
  for (g in 1:n_grids_total) {
    for (t in 1:n_dates) {
      real log_lambda = beta0 + 
                        beta_temp * z_temp[g, t] +
                        beta_rain * z_rain[g, t] +
                        //beta_buildings * z_buildings[g] +
                        //beta_ndvi * z_ndvi[g, t] +*/
                        dot_product(beta_land, z_land[g, ]);
      
      lambda_grid[g, t] = exp(log_lambda);
      lambda_thinned[g, t] = lambda_grid[g, t] * area_grid[g] * p_thin[g];
    }
  }
  
  for (t in 1:n_dates) {
    background[t] = sum(lambda_thinned[, t]) / n_obs_po;
  }
}

model {
  // Priors
  beta0 ~ normal(0, 1);
  beta_temp ~ normal(0, 1);
  beta_rain ~ normal(0, 1);
  // beta_buildings ~ normal(0, 1);
  // beta_ndvi ~ normal(0, 1);
  beta_land ~ normal(0, 1);
  
  alpha0 ~ normal(logit(0.10), 0.5);
  alpha_RH ~ normal(0, 1);
  alpha_WS_rain ~ normal(0, 1);
  alpha_trap_raw ~ normal(0, sigma_trap);
  sigma_trap ~ gamma(2, 2);
  
  delta0 ~ normal(-1, 1);
  delta_poi ~ normal(0, 1);
  delta_reports ~ normal(0, 1);
  
  phi ~ gamma(3, 1.5);
  
  // MARGINALIZED LIKELIHOOD: Survey observations
  // For each observation, marginalize over latent N
  for (k in 1:n_obs_y) {
    int grid_idx = obs_grid_idx[site_to_grid[k]];
    int t = date_y[k];
    real lambda = lambda_grid[grid_idx, t];
    real p_trap = inv_logit(alpha0 + 
                             alpha_RH * z_RH[k] + 
                             alpha_WS_rain * z_WS_rain[k] + 
                             alpha_trap[trap_type[k]]);
    
    target += binomial_negbin_marginal_lpmf(y[k] | p_trap, lambda, phi, N_max);
  }
  
  // Presence-only likelihood
  for (r in 1:n_obs_po) {
    int g = po_grid_idx[r];
    int t = date_po[r];
    target += log(lambda_thinned[g, t]) - log(background[t]) - log(CONSTANT);
  }
}

generated quantities {
  array[n_grids_obs, n_dates] real N_expected;
  
  for (i in 1:n_grids_obs) {
    for (t in 1:n_dates) {
      int g = obs_grid_idx[i];
      N_expected[i, t] = lambda_grid[g, t]; // Mean of NegBin
    }
  }
}
