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
  real partial_sum(array[] int y_slice,
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
      
      lp += binomial_negbin_marginal_lpmf(y_slice[k - start + 1] | p_trap, lambda, phi, N_max);
    }
    return lp;
  }
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
  
  // Grainsize for reduce_sum
  int<lower=1> grainsize;
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
  //  lambda for observed grids to save computation
  matrix[n_grids_total, n_dates] lambda_base;   
  
  // but need lambda_thinned for ALL grids for PO likelihood
  matrix[n_grids_total, n_dates] lambda_thinned;
  vector[n_dates] background;
  
  real<lower=0> sigma_trap = exp(log_sigma_trap);
  vector[5] alpha_trap;
  
  // trap effects sum to zero
  alpha_trap[1:4] = alpha_trap_raw;
  alpha_trap[5] = -sum(alpha_trap_raw);

  // compute lambda_thinned for ALL grids (needed for PO likelihood)
  for (g in 1:n_grids_total) {
    real p_thin = inv_logit(delta0 +
                            delta_poi * z_poi[g] +
                            delta_reports * z_reports[g]);

    for (t in 1:n_dates) {
      real log_lambda_base = beta0 +
                            beta_temp * z_temp[g, t] +
                            beta_rain * z_rain[g, t] +
                            dot_product(beta_land, z_land[g, ]);
      log_lambda_base = fmin(fmax(log_lambda_base, -10), 10);

      lambda_base[g, t] = exp(log_lambda_base);               //  store actual lambda
      lambda_thinned[g, t] = lambda_base[g, t] *             // and thenuse it here too
                             area_grid[g] * p_thin;
    }
  }
  
  // nackground for PO likelihood
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
  target += reduce_sum(
    partial_sum,
    y,
    grainsize,   // grainsize: start with 50 or 100, tune later
    site_to_grid, date_y, trap_type,
    z_RH, z_WS_rain,
    lambda_base,
    alpha0, alpha_RH, alpha_WS_rain, alpha_trap,
    phi,N_max
);

  // for (k in 1:n_obs_y) {
  //   int i = site_to_grid[k];  // Index into observed grids
  //   int t = date_y[k];
  //   real lambda = lambda_base[obs_grid_idx[i], t];
  //   real p_trap = inv_logit(alpha0 + 
  //                            alpha_RH * z_RH[k] + 
  //                            alpha_WS_rain * z_WS_rain[k] + 
  //                            alpha_trap[trap_type[k]]);
  //                         
  //   // target += binomial_negbin_marginal_lpmf(y[k] | p_trap, 
  //   // lambda_base[obs_grid_idx[i], t],  phi,N_max);
  // 
  // }
  
  // PRESENCE-ONLY LIKELIHOOD
  for (r in 1:n_obs_po) {
    int g = po_grid_idx[r];
    int t = date_po[r];
    target += log(lambda_thinned[g, t]) - log(background[t]) - log(CONSTANT);
  }
}

generated quantities {
  // -------------------------------------------------
  // EXPECTED ABUNDANCE FOR ALL GRIDS
  // -------------------------------------------------
  array[n_grids_total, n_dates] real N_expected;

  for (g in 1:n_grids_total) {
    for (t in 1:n_dates) {
      // direct copy of stored latent abundance
      N_expected[g,t] = lambda_base[g,t];
    }
  }

  // -------------------------------------------------
  // POSTERIOR PREDICTIVE REPLICATIONS (COUNTS: SURVEY)
  // -------------------------------------------------
  array[n_obs_y] int y_rep;

  // observed statistics
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

  // replicated stats
  real mean_y_rep;
  real sd_y_rep;
  int  max_y_rep;
  real prop_zeros_rep;

  {
    vector[n_obs_y] y_vec;
    int zeros = 0;

    for (k in 1:n_obs_y) {
      int g = site_to_grid[k];
      int t = date_y[k];

      // draw latent abundance from stored lambda
      int N_sim = neg_binomial_2_rng(lambda_base[g,t], phi);

      // trap detection model
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

  // -------------------------------------------------
  // POSTERIOR PREDICTIVE REPLICATIONS (PO PRESENCE-ONLY)
  // -------------------------------------------------
  array[n_obs_po] int po_rep;

  for (r in 1:n_obs_po) {
    int g = po_grid_idx[r];
    int t = date_po[r];

    // thinning / reporting bias model
    real p_thin = inv_logit(
      delta0
      + delta_poi     * z_poi[g]
      + delta_reports * z_reports[g]
    );

    // stored latent abundance
    int N_sim = neg_binomial_2_rng(lambda_base[g,t], phi);

    // thinned into an observed PO presence
    int N_thinned = binomial_rng(N_sim, p_thin);

    po_rep[r] = (N_thinned > 0);
  }

}

