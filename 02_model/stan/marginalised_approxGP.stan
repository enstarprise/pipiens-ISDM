
//This is the full integrated species distribution model for Cx. pipiens.
// The survey data is modelled through an N-mixture model for the count data,
// the citizen science reports is modelled through a Bernoulli poisson point
// process using the ones trick to account for the background reporting.

// This is using the marginlisation of the latent discrete parameter N.
// Must sum over all possible values, but since it ranges from 0 to infinity, 
// there must be a cap/max...

// In the transformed data block: the N_multiplier defined in stan data list
// (prepare_stan_data.r)

// -----------------------------------------------------------------------------
// HILBERT SPACE GAUSSIAN PROCESS (HSGP): SPATIAL RANDOM EFFECT

// This implements an approximation to a Gaussian Process (GP) with a 
// squared-exponential kernel, following the Hilbert Space GP (HSGP) 
//construction, as defined by Riutort-Mayol(2022):
//
// tHE OVERVIEW:
// A GP defines a prior over smooth functions by specifying how much variation
// is allowed at different spatial scales. Rather than forming a large dense 
// covariance matrix over all locations, the HSGP represents the latent 
// spatial field as a weighted sum of fixed sine-wave basis functions defined 
// over a bounded domain.
//
// BUILDING BLOCKS:
// 1. Sine-wave basis functions (PHI):
//    A fixed set of smooth spatial patterns spanning different spatial scales.
//    These are deterministic functions of space and define *how* the latent
//    field can vary.
//
// 2. Basis frequencies (lambda --> w):
//    Each basis function has an associated frequency that determines how wiggly
//    it is. Lower frequencies = smooth, large-scale variation;
//    higher frequencies = fine-scale variation.
//
// 3. Spectral density (spd / spd_2D):
//    The GP kernel is expressed in the frequency domain to determine how much
//    variance each basis function is allowed to have. This encodes the GP
//    smoothness assumptions (controlled by alpha_gp and rho_gp).
//
// 4. Random coefficients (beta_gp_raw):
//    Independent standard normal variables that introduce randomness. When
//    scaled by the GP-derived standard deviations, they become the coefficients
//    of the basis functions.
//
// HOW THE GP IS CONSTRUCTED:
// The latent spatial field is formed by multiplying the basis matrix (PHI)
// by GP-scaled Gaussian coefficients. Low-frequency basis functions are allowed
// large coefficients, while high-frequency functions are strongly shrunk,
// producing a smooth spatial surface.
//
// The resulting spatial_effect is a draw from a Gaussian Process prior
// (approximately), without ever explicitly constructing a covariance matrix.
// -----------------------------------------------------------------------------


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
    for (k in 1:size(y_slice)) {
      int idx = start + k - 1;
      int i = site_to_grid[idx];
      int t = date_y[idx];
      
      real lambda = lambda_base[i, t];
      real p_trap = inv_logit(
        alpha0 +
        alpha_RH * z_RH[idx] +
        alpha_WS_rain * z_WS_rain[idx] +
        alpha_trap[trap_type[idx]]
        );
      
      lp += binomial_negbin_marginal_lpmf(
        y_slice[k] | p_trap, lambda_base[i, t], phi, N_max);
    }
    return lp;
  }
  
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

  // HSGP HELPER FUNCTIONS (sine basis as in authors' code)
  
  // 1. DEFINE MATHEMATICAL BUILDING BLOCKS
  
  // 1.1 the frequency of a basis function: this step identifies which sine waves 
  // are wiggly (they get penalised) and which are smoother. it computes how
  // fast the sine wave oscillates. to get the weighted sum of the sine waves
  // it requires knowing how wiggly each wave is.
      // m is just an index that labels the sine waves (E.G. 1st 2nd 3rd etc.);
      // it is an ordered list up to M basis functions, defined in stan_data
      // each wave is slightly wigglier than the last...
  real lambda(real L, int m) {
    return ((m * pi()) / (2 * L))^2;
  }
  
  // 1.2 this enforces GP behaviour on the sine waves: in a GP we know how 
  // wiggly the waves need to be, so we assign weights to the waves; "how much
  // weight should a sine wave of frequency *w* get"
      // w is the actual freq of the sine wave labelled by 'm'
      // lambda is the squared frequecy 
  real spd(real alpha, real rho, real w) {
        // Matérn ν=3/2 spectral density in 1D
        // S(w) = alpha^2 * 4 * sqrt(3) * rho / (1 + (rho^2 * w^2) / 3)^2
    return alpha^2 * 4 * sqrt(3) * rho / square(1 + (rho^2 * w^2) / 3);
    }
  
  // 2. BUILD THE BASIS FUNCTIONS 
  
  // this is the GP kernel written in freq space: if the underlying process
  // is smooth ( a GP with lengscale rho), how much variance should a wave of a 
  // given frequency get?
  
  // 2.1 this step creates the sine wave, in dimensions [-L, L] as precomputed
  // based on the regions spatial extent. then creates the 2D basis function 
  // (product of 1D sine functions, i.e. sin(x) * sin(y)). now the sine waves
  // are hills with valleys and crests in a multidimensionl space.
  
  // the sine waves are not the GP yet, but a fixed set of smooth spatial
  // patterns
  vector phi_2D(real L_x, real L_y, int m1, int m2, vector x1, vector x2) {
    vector[rows(x1)] fi1 = 1 / sqrt(L_x) * sin(m1 * pi() * (x1 + L_x) / (2 * L_x));
    vector[rows(x1)] fi2 = 1 / sqrt(L_y) * sin(m2 * pi() * (x2 + L_y) / (2 * L_y));
    return fi1 .* fi2;
  }
  
  // 2.2 spectral density for 2D with a single length-scale (isotropic? 
  // same value in different directions...)
  // the GP kernel (alpha, rho) determines how much variance/wiggly the BFs are
  // allowed to have; this encodes the GP smoothness assumptions controlled by
  // alpha and rho
  real spd_2D(real alpha, real rho, real w1, real w2) {
    real kappa = sqrt(3) / rho;
    real w2_sum = w1^2 + w2^2;
    return alpha^2 * 4 * kappa^3 / square(kappa^2 + w2_sum);
  }
  
 // 3. CONVERT GP KERNEL INTO BASIS WEIGHTS
 
 // 3.1 this step computes how much each sine wave is allowed to contribute. it is
 // the mathematical brigde from GP into weighted sine waves...
 // this is basically asking the GP kernel how much a given freuqnecy of a sine
  // wave is allowed
  
  // this loops over all sine waves, computes their frequencies, and asks the
  // GP kernel how strong each one is allowed to be
  
  // this is still the prior? no randomness or functions drawn
  vector compute_spd_sqrt(real alpha_gp, real rho_gp, int M_sqrt, real L_x, real L_y) {
    int M = M_sqrt * M_sqrt;
    vector[M] spd_sqrt;
    
    for (m1 in 1:M_sqrt) {
      for (m2 in 1:M_sqrt) {
        int idx = (m1 - 1) * M_sqrt + m2;
        real w1 = sqrt(lambda(L_x, m1));
        real w2 = sqrt(lambda(L_y, m2));
        spd_sqrt[idx] = sqrt(spd_2D(alpha_gp, rho_gp, w1, w2));
      }
    }
    return spd_sqrt;
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
  matrix[n_grids_total, n_dates]     z_temp;
  matrix[n_grids_total, n_dates]     z_temp_sq;
  matrix[n_grids_total, n_dates]     z_rain;
  matrix[n_grids_total, n_dates]     z_rain_sq;
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
  
  // HSGP parameters
  int<lower=1> M_sqrt;  // Square root of total basis functions
  matrix[n_grids_total, 2] grid_coords;  // Coordinates in km (centered)
  real L_x;  // Domain boundary in x-direction
  real L_y;  // Domain boundary in y-direction
}

transformed data {
  int max_y = max(y);
  int N_max = max_y + max(250, max_y * N_multiplier);
  
  // PRECOMPUTE THE BASIS MATRIX (PHI)
  // the vectorized version to compute all SPD values. this adds the sine 
  // waves together efficiently, and evaluates every basis function at every 
  // data point, and stores the results in a big matrix.
  
  // this is all the possible smooth shpaes the function is allowed to use
  int M = M_sqrt * M_sqrt;  // Total basis functions
  
  matrix[n_grids_total, M] PHI;
  
  for (m1 in 1:M_sqrt) {
    for (m2 in 1:M_sqrt) {
      int idx = (m1 - 1) * M_sqrt + m2;
      PHI[, idx] = phi_2D(L_x, L_y, m1, m2, 
                          grid_coords[, 1], grid_coords[, 2]);
    }
  }
  
  // Debug: check PHI matrix
  // print("PHI dimensions: ", rows(PHI), " x ", cols(PHI));
  // print("PHI min: ", min(PHI));
  // print("PHI max: ", max(PHI));
  // print("PHI mean: ", mean(PHI));
}

parameters {
  // Abundance model
  real beta0;
  real beta_temp;
  real beta_rain;
  real beta_temp2;  
  real beta_rain2;  
  vector[n_land_covs] beta_land;
  
  // Detection model
  real alpha0;
  real alpha_RH;
  real alpha_WS_rain;
  vector[4] alpha_trap_raw; // four parameters that are free to vary, and the 
  // fifth trap depends on the sum of the other 4 (a5 = -sum(a1...a4))
  real log_sigma_trap;
  
  // Thinning model
  real delta0;
  real delta_poi;
  real delta_reports;
  
  // Dispersion
  //real<lower=0> phi;
  real<lower=1e-6> phi;


  // HSGP parameters (time-independent spatial effect)
  // THE ACTUAL GP 
  // the model parameters (random variables) are what makes the sum of the
  // sine waves random but structured. 
      // beta: random coefficients
      // rho: controls smoothness
      // alpha: controls scales
  real<lower=0> alpha_gp;    // Marginal standard deviation
  real<lower=0> rho_gp;      // Length-scale
  vector[M] beta_gp_raw;     // Standard normal coefficients
}

transformed parameters {
  matrix[n_grids_total, n_dates] lambda_base;
  matrix[n_grids_total, n_dates] lambda_thinned;
  vector[n_dates] background;
  
  real sigma_trap = exp(log_sigma_trap);
  vector[5] alpha_trap;
  
  // Trap effects sum to zero
  alpha_trap[1:4] = alpha_trap_raw;
  alpha_trap[5] = -sum(alpha_trap_raw);
  
  // ---- HSGP CONSTRUCTION: DRAW A SPATIAL GAUSSIAN PROCESS (APPROXIMATION) ----
  // Convert GP kernel parameters (alpha_gp, rho_gp) into frequency-specific
  // standard deviations for each sine-wave basis function. 
  // low-frequency (smooth) waves get large variance; 
  // high-frequency (wiggly) waves get strongly shrunk.
  vector[M] spd_sqrt = compute_spd_sqrt(alpha_gp, rho_gp, M_sqrt, L_x, L_y);
  
  // this creates the spatial latent field (spastial effect) by linearly
  // combining fixed sine-wave basis functions (PHI) with independent
  // gaussian coefficients
  // the coefficients beta_gp_raw are scaled by spd_sqrt so that the resulting
  // function has the same smoothness properties as a GP with a squared-
  // exponential kernel
  vector[n_grids_total] spatial_effect = PHI * (spd_sqrt .* beta_gp_raw);
  
  // Debug: check spatial effect
  // print("Spatial effect min: ", min(spatial_effect));
  // print("Spatial effect max: ", max(spatial_effect));
  // print("Spatial effect mean: ", mean(spatial_effect));
  
  
  // Compute lambda values
  for (g in 1:n_grids_total) {
    real p_thin = inv_logit(delta0 +
                            delta_poi * z_poi[g] +
                            delta_reports * z_reports[g]);

    for (t in 1:n_dates) {
      real log_lambda_base = beta0 +
                             beta_temp  * z_temp[g, t] +
                             beta_temp2 * z_temp_sq[g, t] +
                             beta_rain  * z_rain[g, t] +
                             beta_rain2 * z_rain_sq[g, t] +
                            dot_product(beta_land, z_land[g, ]) +
                            spatial_effect[g];  // Same spatial effect for all times

      lambda_base[g, t] = exp(log_lambda_base);
      lambda_thinned[g, t] = lambda_base[g, t] * 
                             area_grid[g] * p_thin;
    }
  }
  
  // Compute background for PO likelihood
  for (t in 1:n_dates) {
    background[t] = (sum(lambda_thinned[, t])) / n_obs_po;
  }
}

model {
  // ---- DEBUG PRINTS ----
  // print("Model initialization:");
  // print("  beta0 = ", beta0);
  // print("  alpha_gp = ", alpha_gp);
  // print("  rho_gp = ", rho_gp);

  // Priors
  beta0 ~ normal(0, 2);
  beta_temp ~ normal(0, 1);
  beta_rain ~ normal(0, 1);
  beta_temp2 ~ normal(0, 1);  // inactive when z_temp_sq = 0: samples prior only
  beta_rain2 ~ normal(0, 1);  // inactive when z_rain_sq = 0: samples prior only
  beta_land ~ normal(0, 0.2);  
  
  alpha0 ~ normal(logit(0.10), 0.5);
  alpha_RH ~ normal(0, 1);
  alpha_WS_rain ~ normal(0, 1);
  alpha_trap_raw ~ normal(0, sigma_trap);
  log_sigma_trap ~ normal(log(0.5), 0.5); 
  
  delta0 ~ normal(-1, 1);
  delta_poi ~ normal(0, 1);
  delta_reports ~ normal(0, 1);
  
  //phi ~ gamma(2, 1);
  phi ~ lognormal(0, 1);
  
  // HSGP priors (as in authors' code)
  beta_gp_raw ~ std_normal();
  rho_gp ~ gamma(4, 1);      // Authors use gamma(4,1) for rho
  alpha_gp ~ inv_gamma(2, 5); // Authors use inv_gamma(2,5) for alpha

  // MARGINALIZED SURVEY LIKELIHOOD
  target += reduce_sum(
    partial_sum,
    y,
    grainsize,
    site_to_grid, date_y, trap_type,
    z_RH, z_WS_rain,
    lambda_base,
    alpha0, alpha_RH, alpha_WS_rain, alpha_trap,
    phi, N_max
  );
  
  // PRESENCE-ONLY LIKELIHOOD
  for (r in 1:n_obs_po) {
    int g = po_grid_idx[r];
    int t = date_po[r];

    target += log(lambda_thinned[g, t])
            - log(background[t])
            - log(CONSTANT);
  }
}

generated quantities {
  // Posterior predictive for counts
  array[n_obs_y] int y_rep;
  
  // Spatial field for visualization
  vector[n_grids_total] spatial_field = PHI * (spd_sqrt .* beta_gp_raw);
  
  // Monitor N_max adequacy
  int n_near_boundary = 0;
  real max_N_ratio = 0.0;
  
  {
    vector[n_obs_y] y_vec;
    int zeros = 0;

    for (k in 1:n_obs_y) {
      int g = site_to_grid[k];
      int t = date_y[k];

      int N_sim = neg_binomial_2_rng(lambda_base[g, t], phi);
      
      // Monitor N_max adequacy
      real N_ratio = N_sim * 1.0 / N_max;
      if (N_ratio > max_N_ratio) max_N_ratio = N_ratio;
      if (N_ratio > 0.9) n_near_boundary += 1;

      real p_detect = inv_logit(
        alpha0
        + alpha_RH * z_RH[k]
        + alpha_WS_rain * z_WS_rain[k]
        + alpha_trap[trap_type[k]]
      );

      y_rep[k] = binomial_rng(N_sim, p_detect);
      y_vec[k] = y_rep[k];

      if (y_rep[k] == 0) zeros += 1;
    }
  }
  
  // =========================================================================
  // LOG-LIKELIHOOD FOR LOO-CV
  // =========================================================================
  
  // Survey observations
  vector[n_obs_y] log_lik_survey;
  
  for (k in 1:n_obs_y) {
    int g = site_to_grid[k];
    int t = date_y[k];
    real lambda = lambda_base[g, t];
    real p_trap = inv_logit(alpha0 + 
                            alpha_RH * z_RH[k] + 
                            alpha_WS_rain * z_WS_rain[k] + 
                            alpha_trap[trap_type[k]]);
    
    // Marginalized log-likelihood
    log_lik_survey[k] = binomial_negbin_marginal_lpmf(
      y[k] | p_trap, lambda, phi, N_max
    );
  }
  
  // Presence-only observations
  vector[n_obs_po] log_lik_po;
  
  for (r in 1:n_obs_po) {
    int g = po_grid_idx[r];
    int t = date_po[r];
    log_lik_po[r] = log(lambda_thinned[g, t]) - 
                    log(background[t]) - 
                    log(CONSTANT);
  }
}


