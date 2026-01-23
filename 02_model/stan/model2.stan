//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

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
  array[n_obs_po] int<lower=1> ones; // all 1s for Bernoulli trick
  
  // Covariates
  matrix[n_grids_total, n_dates] z_temp;
  matrix[n_grids_total, n_dates] z_rain;
  matrix[n_grids_total, n_dates] z_ndvi;
  vector[n_grids_total] z_buildings;
  matrix[n_grids_total, n_land_covs] z_land;
  vector[n_grids_total] z_poi;
  vector[n_grids_total] z_reports;
  vector[n_obs_y] z_RH;
  vector[n_obs_y] z_WS_rain;
  
  // otherss
  vector[n_grids_total] area_grid;
  real CONSTANT; // for PO likelihood
  
  // for marginalization
  int<lower=1> N_max_multiplier; // e.g., 5 means check up to 5*max(y)
}

transformed data {
  int N_max = max(y) * N_max_multiplier + 50; // upper bound for the sum
}

parameters {
  // Abundance model parameters
  real beta0;
  real beta_temp;
  real beta_rain;
  real beta_buildings;
  real beta_ndvi;
  vector[n_land_covs] beta_land;
  
  // Detection parameters
  real alpha0;
  real alpha_RH;
  real alpha_WS_rain;
  vector[4] alpha_trap_raw;
  real<lower=0> sigma_trap;
  
  // Citizen science thinning
  real delta0;
  real delta_poi;
  real delta_reports;
  
  // Dispersion
  real<lower=0> phi;
}

transformed parameters {
  matrix[n_grids_total, n_dates] lambda_grid;
  matrix[n_grids_total, n_dates] lambda_thinned;
  vector[n_dates] background;
  vector[n_grids_total] p; // thinning probability
  vector[5] alpha_trap;
  
  // constrain trap effects to sum to zero
  alpha_trap[1:4] = alpha_trap_raw;
  alpha_trap[5] = -sum(alpha_trap_raw);
  
  // thinning probability for CS data
  for (g in 1:n_grids_total) {
    p[g] = inv_logit(delta0 + delta_poi * z_poi[g] + delta_reports * z_reports[g]);
  }
  
  // exp abundance at each grid n time
  for (g in 1:n_grids_total) {
    for (t in 1:n_dates) {
      real log_lambda = beta0 + 
                        beta_temp * z_temp[g, t] +
                        beta_rain * z_rain[g, t] +
                        beta_buildings * z_buildings[g] +
                        beta_ndvi * z_ndvi[g, t] +
                        dot_product(beta_land, z_land[g, ]);
      
      lambda_grid[g, t] = exp(log_lambda);
      lambda_thinned[g, t] = lambda_grid[g, t] * area_grid[g] * p[g];
    }
  }
  
  // background exp reports for PO likelihood
  for (t in 1:n_dates) {
    background[t] = sum(lambda_thinned[, t]) / n_obs_po;
  }
}

model {
  // Priors
  beta0 ~ normal(0, 1);
  beta_temp ~ normal(0, 1);
  beta_rain ~ normal(0, 1);
  beta_buildings ~ normal(0, 1);
  beta_ndvi ~ normal(0, 1);
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
  
  // MARGINALIZED LIKELIHOOD FOR SURVEY DATA
  // for each observation, marginalize over all possible N values (but can be inf.... so must set a max....)
  for (k in 1:n_obs_y) {
    int grid_idx = obs_grid_idx[site_to_grid[k]];
    int t = date_y[k];
    real lambda = lambda_grid[grid_idx, t];
    real p_trap = inv_logit(alpha0 + 
                             alpha_RH * z_RH[k] + 
                             alpha_WS_rain * z_WS_rain[k] + 
                             alpha_trap[trap_type[k]]);
    
    //  NegBin parameters; using th eNB2 param
    real prob_nb = phi / (phi + lambda);
    
    // marginalize: sum over N from y[k] to N_max
    vector[N_max - y[k] + 1] log_components;
    
    for (n_idx in 1:(N_max - y[k] + 1)) {
      int N = y[k] + n_idx - 1; // N starts at y[k]
      
      // log[y | N, p]
      real log_binom = binomial_lpmf(y[k] | N, p_trap);
      
      // log[N | lambda, phi]
      real log_negbin = neg_binomial_2_lpmf(N | lambda, phi);
      
      log_components[n_idx] = log_binom + log_negbin;
    }
    
    // Marginalized log-likelihood
    target += log_sum_exp(log_components);
  }
  
  // PRESENCE-ONLY LIKELIHOOD (unchanged)
  for (r in 1:n_obs_po) {
    int g = po_grid_idx[r];
    int t = date_po[r];
    
    target += log(lambda_thinned[g, t]) - log(background[t]) - log(CONSTANT);
  }
}

generated quantities {
  // Optional: compute posterior predictive checks or expected N
  array[n_grids_obs, n_dates] real N_expected;
  
  for (i in 1:n_grids_obs) {
    for (t in 1:n_dates) {
      int g = obs_grid_idx[i];
      N_expected[i, t] = lambda_grid[g, t]; // Expected value under NegBin
    }
  }
}
