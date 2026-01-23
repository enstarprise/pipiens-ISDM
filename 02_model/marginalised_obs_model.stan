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



functions {
  // Marginalized likelihood: sum over latent N
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
  // Dimensions
  int<lower=1> n_grids_total;     // Total grids in dataset
  int<lower=1> n_grids_obs;       // Grids with survey observations
  int<lower=1> n_dates;
  int<lower=1> n_obs_y;           // Number of survey observations
  int<lower=1> n_obs_po;          // Number of presence-only records
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
  int<lower=1> N_multiplier;
}

transformed data {
  int N_max = max(y) * N_multiplier + 100; 
}

parameters {
  // Abundance model
  real beta0;
  real beta_temp;
  real beta_rain;
  // real beta_buildings;
  // real beta_ndvi;
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
  real<lower=1e-12> phi;  // directly estimate phi

}

transformed parameters {
  // Only compute lambda for observed grids to save computation
  matrix[n_grids_obs, n_dates] lambda_obs;
  
  // But need lambda_thinned for ALL grids for PO likelihood
  matrix[n_grids_total, n_dates] lambda_thinned;
  vector[n_dates] background;
  
  real<lower=0> sigma_trap = exp(log_sigma_trap);
  vector[5] alpha_trap;
  
  // Trap effects sum to zero
  alpha_trap[1:4] = alpha_trap_raw;
  alpha_trap[5] = -sum(alpha_trap_raw);
  
  // Compute lambda only for observed grids (for efficiency)
  for (i in 1:n_grids_obs) {
    int g = obs_grid_idx[i];
    for (t in 1:n_dates) {
      real log_lambda_raw = beta0 + 
                            beta_temp * z_temp[g, t] +
                            beta_rain * z_rain[g, t] +
                            // beta_buildings * z_buildings[g] +
                            // beta_ndvi * z_ndvi[g, t] +
                            dot_product(beta_land, z_land[g, ]);
      // Constrain to reasonable range to prevent overflow/underflow
      real log_lambda = fmin(fmax(log_lambda_raw, -10), 10);
      lambda_obs[i, t] = exp(log_lambda);
    }
  }
  
  // Compute lambda_thinned for ALL grids (needed for PO likelihood)
  for (g in 1:n_grids_total) {
    real p_thin = inv_logit(delta0 + 
                            delta_poi * z_poi[g] + 
                            delta_reports * z_reports[g]);
    
    for (t in 1:n_dates) {
      real log_lambda_raw = beta0 + 
                            beta_temp * z_temp[g, t] +
                            beta_rain * z_rain[g, t] +
                            // beta_buildings * z_buildings[g] +
                            // beta_ndvi * z_ndvi[g, t] +
                            dot_product(beta_land, z_land[g, ]);
      // Constrain to reasonable range
      real log_lambda = fmin(fmax(log_lambda_raw, -10), 10);
      real lambda_g = exp(log_lambda);
      lambda_thinned[g, t] = lambda_g * area_grid[g] * p_thin;
    }
  }
  
  // Background for PO likelihood
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
  log_sigma_trap ~ normal(log(0.5), 0.5); 

  
  delta0 ~ normal(-1, 1);
  delta_poi ~ normal(0, 1);
  delta_reports ~ normal(0, 1);
  
  phi ~ normal(2, 1);  // or lognormal, gamma, etc.
  
  // MARGINALIZED SURVEY LIKELIHOOD
  for (k in 1:n_obs_y) {
    int i = site_to_grid[k];  // Index into observed grids
    int t = date_y[k];
    real lambda = lambda_obs[i, t];
    real p_trap = inv_logit(alpha0 + 
                             alpha_RH * z_RH[k] + 
                             alpha_WS_rain * z_WS_rain[k] + 
                             alpha_trap[trap_type[k]]);
    
    target += binomial_negbin_marginal_lpmf(y[k] | p_trap, lambda, phi, N_max);
  }
  
  // PRESENCE-ONLY LIKELIHOOD
  for (r in 1:n_obs_po) {
    int g = po_grid_idx[r];
    int t = date_po[r];
    target += log(lambda_thinned[g, t]) - log(background[t]) - log(CONSTANT);
  }
}

generated quantities {
  // Expected abundance at observed grids
  array[n_grids_obs, n_dates] real N_expected;
  
  for (i in 1:n_grids_obs) {
    for (t in 1:n_dates) {
      N_expected[i, t] = lambda_obs[i, t];
    }
  }
}
