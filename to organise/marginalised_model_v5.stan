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
  // Compute N_max based on negative binomial quantile
  // This finds N such that P(X <= N | lambda, phi) >= prob
  int compute_N_max_quantile(real lambda, real phi, real prob) {
    int N = 0;
    real log_cum_prob = negative_infinity();
    real target_log_prob = log(prob);
    int max_iter = 10000;  // Safety limit
    
    // Start from mode/mean of negative binomial as initial guess
    int N_start = fmax(0, (int) ceil(lambda));
    N = N_start;
    
    // Compute cumulative probability up to N_start
    for (n in 0:N_start) {
      log_cum_prob = log_sum_exp(log_cum_prob, 
                                 neg_binomial_2_lpmf(n | lambda, phi));
    }
    
    // Continue until we hit target probability
    while (log_cum_prob < target_log_prob && N < max_iter) {
      N += 1;
      log_cum_prob = log_sum_exp(log_cum_prob, 
                                 neg_binomial_2_lpmf(N | lambda, phi));
    }
    
    return N;
  }
  
  // More efficient version: use analytical approximation for initial guess
  int compute_N_max_quantile_fast(real lambda, real phi, real prob) {
    // For negative binomial with mean=lambda, var=lambda + lambda^2/phi
    // Use normal approximation for large lambda
    
    if (lambda > 20) {
      // Normal approximation: use quantile formula
      real variance = lambda + lambda * lambda / phi;
      real sigma = sqrt(variance);
      // For prob=0.999, z ≈ 3.09; for prob=0.9999, z ≈ 3.72
      real z = inv_Phi(prob);
      int N_approx = (int) ceil(lambda + z * sigma);
      return fmax(N_approx, (int) ceil(lambda * 2));  // At least 2x mean
    } else {
      // For small lambda, use exact computation
      int N = 0;
      real log_cum_prob = negative_infinity();
      real target_log_prob = log(prob);
      
      while (log_cum_prob < target_log_prob && N < 1000) {
        log_cum_prob = log_sum_exp(log_cum_prob, 
                                   neg_binomial_2_lpmf(N | lambda, phi));
        N += 1;
      }
      
      return N;
    }
  }
  
  // Marginalized likelihood with adaptive N_max
  real binomial_negbin_marginal_lpmf(int y_obs, real p_detect, 
                                     real lambda, real phi, 
                                     real quantile_prob) {
    // Compute adaptive N_max based on the negative binomial distribution
    int N_max = compute_N_max_quantile_fast(lambda, phi, quantile_prob);
    
    // Ensure N_max >= y_obs (must be able to observe what we saw)
    N_max = fmax(N_max, y_obs);
    
    // Safety check: if N_max is unreasonably large, cap it
    N_max = fmin(N_max, 10000);
    
    int n_terms = N_max - y_obs + 1;
    vector[n_terms] log_components;
    
    // Marginalize over latent N
    for (n in 1:n_terms) {
      int N = y_obs + n - 1;
      log_components[n] = binomial_lpmf(y_obs | N, p_detect) + 
                          neg_binomial_2_lpmf(N | lambda, phi);
    }
    
    return log_sum_exp(log_components);
  }
  
  // Parallelized survey likelihood
  real partial_sum_survey(array[] int y_slice,
                         int start, int end,
                         array[] int site_to_grid,
                         array[] int date_y,
                         array[] int trap_type,
                         vector z_RH,
                         vector z_WS_rain,
                         matrix lambda_base,
                         real alpha0,
                         real alpha_RH,
                         real alpha_WS_rain,
                         vector alpha_trap,
                         real phi,
                         real quantile_prob) {
    real lp = 0;
    for (k in start:end) {
      int i = site_to_grid[k];
      int t = date_y[k];
      real lambda = lambda_base[i, t];
      real p_trap = inv_logit(alpha0 + 
                              alpha_RH * z_RH[k] + 
                              alpha_WS_rain * z_WS_rain[k] + 
                              alpha_trap[trap_type[k]]);
      
      // Each observation gets its own adaptive N_max!
      lp += binomial_negbin_marginal_lpmf(y_slice[k - start + 1] | 
                                          p_trap, lambda, phi, quantile_prob);
    }
    return lp;
  }
  
  // Parallelized PO likelihood
  real partial_sum_po(array[] int ones_slice,
                     int start, int end,
                     array[] int po_grid_idx,
                     array[] int date_po,
                     matrix lambda_thinned,
                     vector background,
                     real CONSTANT) {
    real lp = 0;
    for (r in start:end) {
      int g = po_grid_idx[r];
      int t = date_po[r];
      lp += log(lambda_thinned[g, t]) - log(background[t]) - log(CONSTANT);
    }
    return lp;
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
  array[n_obs_y] int<lower=1, upper=n_grids_total> site_to_grid;
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
  matrix[n_grids_total, n_land_covs] z_land;
  vector[n_grids_total] z_poi;
  vector[n_grids_total] z_reports;
  vector[n_obs_y] z_RH;
  vector[n_obs_y] z_WS_rain;
  
  // Other
  vector[n_grids_total] area_grid;
  real CONSTANT;
  
  // Tuning parameters
  real<lower=0.9, upper=0.99999> quantile_prob;  // Typically 0.999 or 0.9999
  int<lower=1> grainsize_survey;
  int<lower=1> grainsize_po;
}

transformed data {
  print("Model Configuration:");
  print("  n_grids_total: ", n_grids_total);
  print("  n_dates: ", n_dates);
  print("  n_obs_y: ", n_obs_y);
  print("  n_obs_po: ", n_obs_po);
  print("  quantile_prob for adaptive N_max: ", quantile_prob);
  print("  grainsize_survey: ", grainsize_survey);
  print("  grainsize_po: ", grainsize_po);
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
}

transformed parameters {
  matrix[n_grids_total, n_dates] lambda_base;   
  matrix[n_grids_total, n_dates] lambda_thinned;
  vector[n_dates] background;
  
  real<lower=0> sigma_trap = exp(log_sigma_trap);
  vector[5] alpha_trap;
  
  alpha_trap[1:4] = alpha_trap_raw;
  alpha_trap[5] = -sum(alpha_trap_raw);

  // Compute lambda_base with pre-computed invariants
  {
    vector[n_grids_total] land_contrib = z_land * beta_land;
    vector[n_grids_total] p_thin_vec;
    
    for (g in 1:n_grids_total) {
      p_thin_vec[g] = inv_logit(delta0 +
                                delta_poi * z_poi[g] +
                                delta_reports * z_reports[g]);
    }
    
    for (t in 1:n_dates) {
      for (g in 1:n_grids_total) {
        real log_lambda_base = beta0 +
                              beta_temp * z_temp[g, t] +
                              beta_rain * z_rain[g, t] +
                              land_contrib[g];
        
        log_lambda_base = fmin(fmax(log_lambda_base, -10), 10);
        
        lambda_base[g, t] = exp(log_lambda_base);
        lambda_thinned[g, t] = lambda_base[g, t] * area_grid[g] * p_thin_vec[g];
      }
    }
    
    for (t in 1:n_dates) {
      background[t] = sum(lambda_thinned[, t]) / n_obs_po;
    }
  }
}

model {
  // Priors
  beta0 ~ normal(0, 1);
  beta_temp ~ normal(0, 1);
  beta_rain ~ normal(0, 1);
  beta_land ~ normal(0, 1);
  
  alpha0 ~ normal(logit(0.10), 0.5);
  alpha_RH ~ normal(0, 1);
  alpha_WS_rain ~ normal(0, 1);
  alpha_trap_raw ~ normal(0, sigma_trap);
  log_sigma_trap ~ normal(log(0.5), 0.5); 
  
  delta0 ~ normal(-1, 1);
  delta_poi ~ normal(0, 1);
  delta_reports ~ normal(0, 1);
  
  phi ~ gamma(2, 1);
  
  // PARALLELIZED SURVEY LIKELIHOOD with adaptive N_max
  target += reduce_sum(
    partial_sum_survey,
    y,
    grainsize_survey,
    site_to_grid, date_y, trap_type,
    z_RH, z_WS_rain,
    lambda_base,
    alpha0, alpha_RH, alpha_WS_rain, alpha_trap,
    phi, quantile_prob
  );

  // PARALLELIZED PO LIKELIHOOD
  target += reduce_sum(
    partial_sum_po,
    ones,
    grainsize_po,
    po_grid_idx, date_po,
    lambda_thinned,
    background,
    CONSTANT
  );
}

generated quantities {
  array[n_obs_y] int y_rep;
  
  real mean_y_obs = mean(to_vector(y));
  real sd_y_obs   = sd(to_vector(y));
  int  max_y_obs  = max(y);
  real prop_zeros_obs;
  
  {
    int zeros = 0;
    for (k in 1:n_obs_y)
      if (y[k] == 0) zeros += 1;
    prop_zeros_obs = zeros * 1.0 / n_obs_y;
  }

  real mean_y_rep;
  real sd_y_rep;
  int  max_y_rep;
  real prop_zeros_rep;
  
  // NEW: Monitor adaptive N_max behavior
  real mean_N_max = 0;
  real max_N_max_used = 0;
  real min_N_max_used = 10000;

  {
    vector[n_obs_y] y_vec;
    int zeros = 0;
    real sum_N_max = 0;

    for (k in 1:n_obs_y) {
      int g = site_to_grid[k];
      int t = date_y[k];
      
      real lambda = lambda_base[g, t];
      
      // Track N_max that would be used for this observation
      int N_max_k = compute_N_max_quantile_fast(lambda, phi, quantile_prob);
      N_max_k = fmax(N_max_k, y[k]);
      N_max_k = fmin(N_max_k, 10000);
      
      sum_N_max += N_max_k;
      if (N_max_k > max_N_max_used) max_N_max_used = N_max_k;
      if (N_max_k < min_N_max_used) min_N_max_used = N_max_k;

      // Generate replicate
      int N_sim = neg_binomial_2_rng(lambda, phi);

      real p_detect = inv_logit(
        alpha0
        + alpha_RH      * z_RH[k]
        + alpha_WS_rain * z_WS_rain[k]
        + alpha_trap[trap_type[k]]
      );

      y_rep[k] = binomial_rng(N_sim, p_detect);
      y_vec[k] = y_rep[k];

      if (y_rep[k] == 0) zeros += 1;
    }

    mean_y_rep = mean(y_vec);
    sd_y_rep   = sd(y_vec);
    max_y_rep  = max(y_rep);
    prop_zeros_rep = zeros * 1.0 / n_obs_y;
    mean_N_max = sum_N_max / n_obs_y;
  }

  // PO replications
  array[n_obs_po] int po_rep;

  for (r in 1:n_obs_po) {
    int g = po_grid_idx[r];
    int t = date_po[r];

    real p_thin = inv_logit(
      delta0
      + delta_poi     * z_poi[g]
      + delta_reports * z_reports[g]
    );

    int N_sim = neg_binomial_2_rng(lambda_base[g,t], phi);
    int N_thinned = binomial_rng(N_sim, p_thin);

    po_rep[r] = (N_thinned > 0);
  }
}


 