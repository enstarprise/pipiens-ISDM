// =============================================================================
// NEAREST-NEIGHBOUR GAUSSIAN PROCESS (NNGP): SPATIAL RANDOM EFFECT
//
// latent NNGP following:
//   Zhang (2018) https://mc-stan.org/learn-stan/case-studies/nngp.html
//   Datta et al. (2016), Finley et al. (2017)
//
// The spatial random effect  spatial_effect[g] is given a GP prior whose precision
// matrix is the NNGP sparse approximation:
//
//   C*^{-1} = (I - A*)^T D*^{-1} (I - A*)
//
// KERNEL: Exponential (Matérn ν=0.5)
//   C(s_i, s_j) = sigmasq * exp(-phi * dist(i, j))
//
// where sigmasq = alpha_gp^2 and phi = 1/rho_gp 
// =============================================================================


functions {

  real nngp_w_lpdf(vector w, real sigmasq, real phi,
                   matrix NN_dist, matrix NN_distM,
                   array[,] int NN_ind, int N, int M) {

    vector[N] V;
    vector[N] I_Aw = w;   // will subtract kriging means in-place
    int h;

    for (i in 2:N) {
      // dim: actual number of neighbours for location i.
      // Early locations (i < M+1) have fewer than M predecessors.
      int dim = (i < (M + 1)) ? (i - 1) : M;

      matrix[dim, dim] iNNdistM;  // correlation matrix among neighbours
      vector[dim]      iNNcorr;   // correlation from i to its neighbours
      vector[dim]      v;
      row_vector[dim]  v2;

      // Build correlation matrix among the dim neighbours
      if (dim == 1) {
        iNNdistM[1, 1] = 1.0;
      } else {
        h = 0;
        for (j in 1:(dim - 1)) {
          for (k in (j + 1):dim) {
            h += 1;
            iNNdistM[j, k] = exp(-phi * NN_distM[(i - 1), h]);
            iNNdistM[k, j] = iNNdistM[j, k];
          }
        }
        for (j in 1:dim) {
          iNNdistM[j, j] = 1.0;
        }
      }

      // Cholesky of the neighbour correlation matrix
      matrix[dim, dim] iNNCholL = cholesky_decompose(iNNdistM);

      // Correlation from location i to its neighbours
      iNNcorr = to_vector(exp(-phi * NN_dist[(i - 1), 1:dim]));

      // Forward solve: v = L^{-1} * iNNcorr
      v = mdivide_left_tri_low(iNNCholL, iNNcorr);

      // Conditional variance (in correlation units): 1 - v'v
      V[i] = 1.0 - dot_self(v);

      // Kriging weights: v2 = v' * L^{-1}
      v2 = mdivide_right_tri_low(v', iNNCholL);

      // Residual: w[i] minus its kriging mean from neighbours
      I_Aw[i] -= v2 * w[NN_ind[(i - 1), 1:dim]];
    }

    // Location 1 has no predecessors: marginal variance = sigmasq
    V[1] = 1.0;

    return -0.5 * (
      (1.0 / sigmasq) * dot_product(I_Aw, (I_Aw ./ V)) +
      sum(log(V)) +
      N * log(sigmasq)
    );
  }

  // ---------------------------------------------------------------------------
  //MARGINALISED BINOMIAL-NEGBIN LIKELIHOOD 
  // ---------------------------------------------------------------------------
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
                   real phi_nb,
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
      lp += binomial_negbin_marginal_lpmf(
        y_slice[k] | p_trap, lambda_base[i, t], phi_nb, N_max);
    }
    return lp;
  }

  real binomial_negbin_marginal_lpmf(int y_obs, real p_detect,
                                     real lambda, real phi_nb, int N_max) {
    vector[N_max - y_obs + 1] log_components;
    for (n in 1:(N_max - y_obs + 1)) {
      int N_val = y_obs + n - 1;
      log_components[n] = binomial_lpmf(y_obs | N_val, p_detect) +
                          neg_binomial_2_lpmf(N_val | lambda, phi_nb);
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
  matrix[n_grids_total, n_dates]     z_temp_sq;
  matrix[n_grids_total, n_dates]     z_rain;
  matrix[n_grids_total, n_dates]     z_rain_sq;
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

  // --- nngp neighbour structure ---
  int<lower=1> M;  // Number of nearest neighbours

  // Neighbour indices: row i-1 gives M predecessor indices for location i
  array[n_grids_total - 1, M] int<lower=1, upper=n_grids_total> NN_ind;

  // Distances from each location to its M neighbours
  matrix[n_grids_total - 1, M] NN_dist;

  // Pairwise distances among the M neighbours (lower-tri, stored row-wise)
  matrix[n_grids_total - 1, (M * (M - 1) / 2)] NN_distM;
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
  real beta_temp2;  // quadratic temperature; samples from prior when z_temp_sq = 0
  real beta_rain2;  // quadratic rainfall;    samples from prior when z_rain_sq = 0
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
  real<lower=1e-6> phi;

  // --- NNGP GP hyperparameters ---
  real<lower=0> alpha_gp;  // Marginal SD: sqrt(sigmasq) in the kernel
  real<lower=0> rho_gp;    // Length-scale in km

  // --- spatial field ---
  vector[n_grids_total] spatial_effect;
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
                             beta_temp  * z_temp[g, t] +
                             beta_temp2 * z_temp_sq[g, t] +
                             beta_rain  * z_rain[g, t] +
                             beta_rain2 * z_rain_sq[g, t] +
                             dot_product(beta_land, z_land[g, ]) +
                             spatial_effect[g];  

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
  beta_temp  ~ normal(0, 1);
  beta_rain  ~ normal(0, 1);
  beta_temp2 ~ normal(0, 1);  // inactive when z_temp_sq = 0: samples prior only
  beta_rain2 ~ normal(0, 1);  // inactive when z_rain_sq = 0: samples prior only

  alpha0        ~ normal(logit(0.10), 0.5);
  alpha_RH      ~ normal(0, 1);
  alpha_WS_rain ~ normal(0, 1);
  alpha_trap_raw ~ normal(0, sigma_trap);
  log_sigma_trap ~ normal(log(0.5), 0.5);

  delta0        ~ normal(-1, 1);
  delta_poi     ~ normal(0, 1);
  delta_reports ~ normal(0, 1);

  phi ~ lognormal(0, 1);

  // ---nngp priors ---
  // marginal SD. inv_gamma(2,5) keeps it from exploding 
  alpha_gp ~ inv_gamma(2, 5);

  // length-scale in km: gamma(4,1) centres it around 4km
  // gamma(8, 0.5) centres it around 16km
  rho_gp ~ gamma(8, 0.5);

  // nngp prior on the latent spatial field.
   spatial_effect ~ nngp_w(square(alpha_gp), 1.0 / rho_gp,
                NN_dist, NN_distM, NN_ind,
                n_grids_total, M);

  // --- survey likelihood ---
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

  vector[n_grids_total] spatial_field = spatial_effect;

  int  n_near_boundary = 0;
  real max_N_ratio     = 0.0;

  {
    for (k in 1:n_obs_y) {
      int g = site_to_grid[k];
      int t = date_y[k];

      int N_sim = neg_binomial_2_rng(lambda_base[g, t], phi);

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
    log_lik_survey[k] = binomial_negbin_marginal_lpmf(
      y[k] | p_trap, lambda_base[g, t], phi, N_max
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

