
// poisson

functions {
  // MARGINALISED BINOMIAL-NEGBIN LIKELIHOOD
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
                   int N_max) {
    real lp = 0;
    for (k in 1:size(y_slice)) {
      int idx = start + k - 1;
      int i   = site_to_grid[idx];
      int t   = date_y[idx];
      real p_trap = inv_logit(
        alpha0 +
        alpha_RH       * z_RH[idx] +
        alpha_WS_rain  * z_WS_rain[idx] +
        alpha_trap[trap_type[idx]]
      );
      lp += binomial_poisson_marginal_lpmf(
        y_slice[k] | p_trap, lambda_base[i, t], N_max);
    }
    return lp;
  }

  real binomial_poisson_marginal_lpmf(int y_obs, real p_detect,
                                     real lambda, int N_max) {
    vector[N_max - y_obs + 1] log_components;
    for (n in 1:(N_max - y_obs + 1)) {
      int N_val = y_obs + n - 1;
      log_components[n] = binomial_lpmf(y_obs | N_val, p_detect) +
                          poisson_log_lpmf(N_val | lambda);
    }
    return log_sum_exp(log_components);
  }

}


data {
  // --- dimensions ---
  int<lower=1> n_grids_total;
  int<lower=1> n_grids_obs;
  int<lower=1> n_dates;
  int<lower=1> n_obs_y;
  int<lower=1> n_obs_po;
  int<lower=1> n_land_covs;

  // --- indices ---
  array[n_grids_obs] int<lower=1, upper=n_grids_total> obs_grid_idx;
  array[n_obs_y]     int<lower=1, upper=n_grids_total> site_to_grid;
  array[n_obs_y]     int<lower=1, upper=n_dates>       date_y;
  array[n_obs_y]     int<lower=1, upper=5>             trap_type;
  array[n_obs_po]    int<lower=1, upper=n_grids_total> po_grid_idx;
  array[n_obs_po]    int<lower=1, upper=n_dates>       date_po;

  // --- thee data observations ---
  array[n_obs_y]  int<lower=0> y;
  array[n_obs_po] int<lower=1> ones;

  // --- covariates ---
  matrix[n_grids_total, n_dates]     z_temp;
  matrix[n_grids_total, n_dates]     z_rain;
  matrix[n_grids_total, n_land_covs] z_land;
  vector[n_grids_total]              z_poi;
  vector[n_grids_total]              z_reports;
  vector[n_obs_y]                    z_RH;
  vector[n_obs_y]                    z_WS_rain;

  // --- other ---
  vector[n_grids_total] area_grid;
  real  CONSTANT;
  int<lower=1> N_multiplier;
  int<lower=1> grainsize;
}


transformed data {
  int max_y = max(y);
  int N_max = max_y + max(250, max_y * N_multiplier);
}


parameters {
  // --- Abundance model  ---
  real beta0;
  real beta_temp;
  real beta_rain;
  vector[n_land_covs] beta_land;

  // --- Detection model  ---
  real alpha0;
  real alpha_RH;
  real alpha_WS_rain;
  vector[4] alpha_trap_raw;
  real log_sigma_trap;

  // --- Thinning model  ---
  real delta0;
  real delta_poi;
  real delta_reports;

  // --- Dispersion  ---
  // POISSON
 }


transformed parameters {
  matrix[n_grids_total, n_dates] lambda_base;
  matrix[n_grids_total, n_dates] lambda_thinned;
  vector[n_dates] background;

  real sigma_trap = exp(log_sigma_trap);
  vector[5] alpha_trap;

  // Trap effects sum to zero 
  alpha_trap[1:4] = alpha_trap_raw;
  alpha_trap[5]   = -sum(alpha_trap_raw);

  for (g in 1:n_grids_total) {
    real p_thin = inv_logit(delta0 +
                            delta_poi     * z_poi[g] +
                            delta_reports * z_reports[g]);

    for (t in 1:n_dates) {
      real log_lambda_base = beta0 +
                             beta_temp * z_temp[g, t] +
                             beta_rain * z_rain[g, t] +
                             dot_product(beta_land, z_land[g, ]);

      lambda_base[g, t]    = exp(log_lambda_base);
      lambda_thinned[g, t] = lambda_base[g, t] * area_grid[g] * p_thin;
    }
  }

  for (t in 1:n_dates) {
    background[t] = sum(lambda_thinned[, t]) / n_obs_po;
  }
}


model {
  // --- priors ---
  beta0     ~ normal(0, 2);
  beta_temp ~ normal(0, 1);
  beta_rain ~ normal(0, 1);
  beta_land ~ normal(0, 0.2);

  alpha0        ~ normal(logit(0.10), 0.5);
  alpha_RH      ~ normal(0, 1);
  alpha_WS_rain ~ normal(0, 1);
  alpha_trap_raw ~ normal(0, sigma_trap);
  log_sigma_trap ~ normal(log(0.5), 0.5);

  delta0        ~ normal(-1, 1);
  delta_poi     ~ normal(0, 1);
  delta_reports ~ normal(0, 1);

  // --- survey likelihood ---
  target += reduce_sum(
    partial_sum,
    y,
    grainsize,
    site_to_grid, date_y, trap_type,
    z_RH, z_WS_rain,
    lambda_base,
    alpha0, alpha_RH, alpha_WS_rain, alpha_trap,
    N_max
  );

  // --- presence-only likelihood ---
  for (r in 1:n_obs_po) {
    int g = po_grid_idx[r];
    int t = date_po[r];
    target += log(lambda_thinned[g, t])
            - log(background[t])
            - log(CONSTANT);
  }
}


generated quantities {
  array[n_obs_y] int y_rep;

  int  n_near_boundary = 0;
  real max_N_ratio     = 0.0;

  {
    for (k in 1:n_obs_y) {
      int g = site_to_grid[k];
      int t = date_y[k];

      int N_sim = poisson_rng(lambda_base[g, t]);

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
    }
  }

  // =========================================================================
  // LOG-LIKELIHOOD FOR LOO-CV
  // =========================================================================
  vector[n_obs_y]  log_lik_survey;
  vector[n_obs_po] log_lik_po;

  for (k in 1:n_obs_y) {
    int g = site_to_grid[k];
    int t = date_y[k];
    real p_trap = inv_logit(alpha0 +
                            alpha_RH      * z_RH[k] +
                            alpha_WS_rain * z_WS_rain[k] +
                            alpha_trap[trap_type[k]]);
    log_lik_survey[k] = binomial_poisson_marginal_lpmf(
      y[k] | p_trap, lambda_base[g, t], N_max
    );
  }

  for (r in 1:n_obs_po) {
    int g = po_grid_idx[r];
    int t = date_po[r];
    log_lik_po[r] = log(lambda_thinned[g, t]) -
                    log(background[t]) -
                    log(CONSTANT);
  }
}

