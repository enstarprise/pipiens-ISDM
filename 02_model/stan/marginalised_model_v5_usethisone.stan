//This is the full integrated species distribution model for Cx. pipiens.
// The survey data is modelled through an N-mixture model for the count data,
// the citizen science reports is modelled through a Bernoulli poisson point
// process using the ones trick to account for the background reporting.

// This is using the marginlisation of the latent discrete parameter N.
// Must sum over all possible values, but since it ranges from 0 to infinity, 
// there must be a cap/max...

// This uses a CDF to find the value (x) in a given distribution that 
// corresponds to a given cumulative probability (i.e. the prob that it will
// take a value less than or equal to it/ summing all the prob. up to that 
// value)

// In the transformed data block: the quantile to sum up to, defined in stan
// data list (prepare_stan_data.r)
 
// https://mc-stan.org/docs/stan-users-guide/latent-discrete.html
// https://mbjoseph.github.io/posts/2020-04-28-a-step-by-step-guide-to-marginalizing-over-discrete-parameters-for-ecologists-using-stan/

// this is the updated version from v2... 

      // v2 is another version of this one (v5), but a much faster:
            // v2 (explicit marginalization) using a global N max based on 
            // observed values of the survey data
            
functions {

  // -------------------------------------------------
  // the inverse CDF for NegBin2 using lcdf. the function
  // returns the smallest integer where N ≥ y_obs 
  // such that 
  
  // P(N <= N_max) >= p
      // where p = 0.99, mu = lambda, phi = disp, 
      //    y_obs = obs count
      // i.e. the 99th percentile of the NB dist, but 
      // never smaller than the obs count
      
  // the output is an integer: N_max
  
  // -------------------------------------------------
  int negbin_quantile(real p, real mu, real phi, int y_obs) {
    
    // initialise the smallest possible N
    int n = y_obs; 
    
    // below i.e. log P(N <= n | lambda, phi), the probability mass at or below n
    // the probability that N is less than or equal to some integer n
      // log bc these probabilities can be really sum (avoid numericalk underflow)
    real cdf = neg_binomial_2_lcdf(n | mu, phi);
    
    // increment until the target prob is reached:
        // check if P(N<= n) <p
        // if so, increase n by 1
        // recompute the cdf
        // repeat until it exceeds p
    while (cdf < p) {
      n += 1;
      cdf = neg_binomial_2_lcdf(n | mu, phi); 
    }
    return n; // this isthe smallest integer satisfying P(N<= n) >= p,
    // whch becomes the per observation N_max
  }

  // -------------------------------------------------
  // Marginalised likelihood over latent N
  // -------------------------------------------------
  real binomial_negbin_marginal_lpmf(
      int y_obs,
      real p_detect,
      real lambda,
      real phi,
      int N_max
  ) {
    vector[N_max - y_obs + 1] log_components;

    for (n in 1:(N_max - y_obs + 1)) {
      int N = y_obs + n - 1;
      log_components[n] =
        binomial_lpmf(y_obs | N, p_detect) +
        neg_binomial_2_lpmf(N | lambda, phi);
    }

    return log_sum_exp(log_components);
  }

  // -------------------------------------------------
  // reduce_sum partial
  // -------------------------------------------------
  real partial_sum(
      array[] int y_slice,
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
      real phi
  ) {
    real lp = 0;

    for (k in 1:(end - start + 1)) {
      int k_global = start + k - 1;

      int g = site_to_grid[k_global];
      int t = date_y[k_global];

      real lambda = lambda_base[g, t];

      real p_detect = inv_logit(
        alpha0 +
        alpha_RH      * z_RH[k_global] +
        alpha_WS_rain * z_WS_rain[k_global] +
        alpha_trap[trap_type[k_global]]
      );

      // for each obs, usr the inverse CDF trunc.
      int N_max_k =
        negbin_quantile(0.99, lambda, phi, y_slice[k]);

      lp += binomial_negbin_marginal_lpmf(
        y_slice[k] |
        p_detect,
        lambda,
        phi,
        N_max_k
      );
    }

    return lp;
  }
}

data {
  int<lower=1> n_grids_total;
  int<lower=1> n_dates;
  int<lower=1> n_obs_y;
  int<lower=1> n_obs_po;
  int<lower=1> n_land_covs;

  array[n_obs_y] int<lower=1, upper=n_grids_total> site_to_grid;
  array[n_obs_y] int<lower=1, upper=n_dates> date_y;
  array[n_obs_y] int<lower=1, upper=5> trap_type;

  array[n_obs_po] int<lower=1, upper=n_grids_total> po_grid_idx;
  array[n_obs_po] int<lower=1, upper=n_dates> date_po;

  array[n_obs_y] int<lower=0> y;
  array[n_obs_po] int<lower=1> ones;

  matrix[n_grids_total, n_dates] z_temp;
  matrix[n_grids_total, n_dates] z_rain;
  matrix[n_grids_total, n_land_covs] z_land;

  vector[n_grids_total] z_poi;
  vector[n_grids_total] z_reports;

  vector[n_obs_y] z_RH;
  vector[n_obs_y] z_WS_rain;

  vector[n_grids_total] area_grid;

  real CONSTANT;
  int<lower=1> grainsize;
}

parameters {
  // abundance
  real beta0;
  real beta_temp;
  real beta_rain;
  vector[n_land_covs] beta_land;

  // detection
  real alpha0;
  real alpha_RH;
  real alpha_WS_rain;
  vector[4] alpha_trap_raw;
  real log_sigma_trap;

  // thinning
  real delta0;
  real delta_poi;
  real delta_reports;

  // dispersion
  real<lower=1e-12> phi;
}

transformed parameters {
  matrix[n_grids_total, n_dates] lambda_base;
  matrix[n_grids_total, n_dates] lambda_thinned;
  vector[n_dates] background;

  real<lower=0> sigma_trap = exp(log_sigma_trap);
  vector[5] alpha_trap;

  alpha_trap[1:4] = alpha_trap_raw;
  alpha_trap[5]   = -sum(alpha_trap_raw);

  for (g in 1:n_grids_total) {
    real p_thin = inv_logit(
      delta0 +
      delta_poi     * z_poi[g] +
      delta_reports * z_reports[g]
    );

    for (t in 1:n_dates) {
      real log_lambda =
        beta0 +
        beta_temp * z_temp[g, t] +
        beta_rain * z_rain[g, t] +
        dot_product(beta_land, z_land[g]);

      log_lambda = fmin(fmax(log_lambda, -10), 10);

      lambda_base[g, t] = exp(log_lambda);
      lambda_thinned[g, t] =
        lambda_base[g, t] * area_grid[g] * p_thin;
    }
  }

  for (t in 1:n_dates)
    background[t] = sum(lambda_thinned[, t]) / n_obs_po;
}

model {
  // priors
  beta0 ~ normal(0, 1);
  beta_temp ~ normal(0, 1);
  beta_rain ~ normal(0, 1);
  beta_land ~ normal(0, 1);

  alpha0 ~ normal(logit(0.1), 0.5);
  alpha_RH ~ normal(0, 1);
  alpha_WS_rain ~ normal(0, 1);
  alpha_trap_raw ~ normal(0, sigma_trap);
  log_sigma_trap ~ normal(log(0.5), 0.5);

  delta0 ~ normal(-1, 1);
  delta_poi ~ normal(0, 1);
  delta_reports ~ normal(0, 1);

  phi ~ gamma(2, 1);

  // survey likelihood (latent N marginalised)
  target += reduce_sum(
    partial_sum,
    y,
    grainsize,
    site_to_grid,
    date_y,
    trap_type,
    z_RH,
    z_WS_rain,
    lambda_base,
    alpha0,
    alpha_RH,
    alpha_WS_rain,
    alpha_trap,
    phi
  );

  // presence-only likelihood
  for (r in 1:n_obs_po) {
    int g = po_grid_idx[r];
    int t = date_po[r];
    target += log(lambda_thinned[g, t])
              - log(background[t])
              - log(CONSTANT);
  }
}



