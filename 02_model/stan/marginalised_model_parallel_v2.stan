// paralellising the survey and PO likelihoods instead of sequential
// lambda calulcations are vectorised


functions {
  // Marginalized likelihood
  real binomial_negbin_marginal_lpmf(int y_obs, real p_detect, 
                                     real lambda, real phi, int N_max) {
    int n_terms = N_max - y_obs + 1;
    vector[n_terms] log_components;
    
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
                         int N_max) {
    real lp = 0;
    for (k in start:end) {
      int i = site_to_grid[k];
      int t = date_y[k];
      real lambda = lambda_base[i, t];
      real p_trap = inv_logit(alpha0 + 
                              alpha_RH * z_RH[k] + 
                              alpha_WS_rain * z_WS_rain[k] + 
                              alpha_trap[trap_type[k]]);
      
      lp += binomial_negbin_marginal_lpmf(y_slice[k - start + 1] | 
                                          p_trap, lambda, phi, N_max);
    }
    return lp;
  }
  
  // UNUSED : vectorized lambda computation using matrix operations
  // faster than nested loops?
  matrix compute_lambda_vectorized(matrix z_temp,
                                  matrix z_rain,
                                  matrix z_land,
                                  real beta0,
                                  real beta_temp,
                                  real beta_rain,
                                  vector beta_land) {
    int n_grids = rows(z_temp);
    int n_dates = cols(z_temp);
    
    matrix[n_grids, n_dates] log_lambda;
    
    // vectorized operations are MUCH faster than loops in Stan?
    log_lambda = rep_matrix(beta0, n_grids, n_dates) +
                 beta_temp * z_temp +
                 beta_rain * z_rain +
                 rep_matrix(z_land * beta_land, n_dates);
    
    // Clamp values
    for (g in 1:n_grids) {
      for (t in 1:n_dates) {
        log_lambda[g, t] = fmin(fmax(log_lambda[g, t], -10), 10);
      }
    }
    
    return exp(log_lambda);
  }
  
  // parallelized PO likelihood -- not really needed since PO is small
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
  int<lower=1> N_multiplier;
  int<lower=1> grainsize_survey;
  int<lower=1> grainsize_po;
  //int<lower=0> N_buffer;
}

transformed data {
  int max_y = max(y);
  int N_max = max_y * N_multiplier;
  //int N_max = max_y * N_multiplier + N_buffer;
  
  // pre-compute z_land * beta_land contribution (grid-level, time-invariant)
  //  !! used in transformed parameters
  
  print("Model Configuration:");
  print("  n_grids_total: ", n_grids_total);
  print("  n_dates: ", n_dates);
  print("  n_obs_y: ", n_obs_y);
  print("  n_obs_po: ", n_obs_po);
  print("  N_max: ", N_max);
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

  // OPTIMIZED: Vectorized computation where possible
  {
    // pre-compute land covariate contribution (shared across all dates, time invariant)
    vector[n_grids_total] land_contrib = z_land * beta_land;
    
    // pre-compute thinning probability for each grid
    vector[n_grids_total] p_thin_vec;
    for (g in 1:n_grids_total) {
      p_thin_vec[g] = inv_logit(delta0 +
                                delta_poi * z_poi[g] +
                                delta_reports * z_reports[g]);
    }
    
    // thrn compute lambda_base with more vectorization
    // main bottleneck for large n_grids_total
    for (t in 1:n_dates) {
      for (g in 1:n_grids_total) {
        real log_lambda_base = beta0 +
                              beta_temp * z_temp[g, t] +
                              beta_rain * z_rain[g, t] +
                              land_contrib[g];  // Pre-computed
        
        log_lambda_base = fmin(fmax(log_lambda_base, -10), 10);
        
        lambda_base[g, t] = exp(log_lambda_base);
        lambda_thinned[g, t] = lambda_base[g, t] * area_grid[g] * p_thin_vec[g];
      }
    }
    
    // vectorized background intensity
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
  
  // PARALLELIZED SURVEY LIKELIHOOD
  target += reduce_sum(
    partial_sum_survey,
    y,
    grainsize_survey,
    site_to_grid, date_y, trap_type,
    z_RH, z_WS_rain,
    lambda_base,
    alpha0, alpha_RH, alpha_WS_rain, alpha_trap,
    phi, N_max
  );

  // PARALLELIZED PO LIKELIHOOD
  // helps when n_obs_po is large
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
  int n_near_boundary = 0;
  real max_N_ratio = 0.0;

  {
    vector[n_obs_y] y_vec;
    int zeros = 0;

    for (k in 1:n_obs_y) {
      int g = site_to_grid[k];
      int t = date_y[k];

      int N_sim = neg_binomial_2_rng(lambda_base[g,t], phi);
      
      real N_ratio = N_sim * 1.0 / N_max;
      if (N_ratio > max_N_ratio) max_N_ratio = N_ratio;
      if (N_ratio > 0.9) n_near_boundary += 1;

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
  }

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

