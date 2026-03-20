// =============================================================================
// LATENT NNGP — NON-CENTRED + INTERCEPT-CENTRED PARAMETERISATION
//
// Combines two fixes for HMC geometry problems:
//
// FIX 1 (intercept-centring): breaks beta0 <-> spatial field ridge
//   w_b1[g] = beta0 + alpha_gp * w_raw[g]
//   beta0 anchors the mean of the field; sampled separately with its own prior
//
// FIX 2 (non-centred): breaks alpha_gp <-> spatial field funnel
//   w_raw ~ NNGP(0, 1) on unit scale — alpha_gp does NOT appear in the prior
//   alpha_gp only appears as a scalar multiplier in transformed parameters
//
// Together: HMC sees a flat, well-conditioned geometry for both beta0 and alpha_gp
//
// NOTE: Presence-only likelihood and thinning component removed.
//       Area offset removed (area_grid absorbed into beta0).
// =============================================================================


functions {

  real nngp_w_lpdf(vector w_raw, real sigmasq, real phi,
                   matrix NN_dist, matrix NN_distM,
                   array[,] int NN_ind, int N, int M) {

    vector[N] V;
    vector[N] I_Aw = w_raw;
    int h;

    for (i in 2:N) {
      int dim = (i < (M + 1)) ? (i - 1) : M;

      matrix[dim, dim] iNNdistM;
      vector[dim]      iNNcorr;
      vector[dim]      v;
      row_vector[dim]  v2;

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
        for (j in 1:dim){
          iNNdistM[j, j] = 1.0;
        }
      }

      matrix[dim, dim] iNNCholL = cholesky_decompose(iNNdistM);
      iNNcorr = to_vector(exp(-phi * NN_dist[(i - 1), 1:dim]));

      v       = mdivide_left_tri_low(iNNCholL, iNNcorr);
      V[i]    = 1.0 - dot_self(v);
      v2      = mdivide_right_tri_low(v', iNNCholL);
      I_Aw[i] -= v2 * w_raw[NN_ind[(i - 1), 1:dim]];
    }
    V[1] = 1.0;

    return -0.5 * (
      (1.0 / sigmasq) * dot_product(I_Aw, (I_Aw ./ V)) +
      sum(log(V)) +
      N * log(sigmasq)
    );
  }

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
        alpha_RH      * z_RH[idx] +
        alpha_WS_rain * z_WS_rain[idx] +
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
  int<lower=1> n_grids_total;
  int<lower=1> n_grids_obs;
  int<lower=1> n_dates;
  int<lower=1> n_obs_y;
  int<lower=1> n_land_covs;

  array[n_grids_obs] int<lower=1, upper=n_grids_total> obs_grid_idx;
  array[n_obs_y]     int<lower=1, upper=n_grids_total> site_to_grid;
  array[n_obs_y]     int<lower=1, upper=n_dates>       date_y;
  array[n_obs_y]     int<lower=1, upper=5>             trap_type;

  array[n_obs_y] int<lower=0> y;

  matrix[n_grids_total, n_dates]     z_temp;
  matrix[n_grids_total, n_dates]     z_temp_sq;
  matrix[n_grids_total, n_dates]     z_rain;
  matrix[n_grids_total, n_dates]     z_rain_sq;
  matrix[n_grids_total, n_land_covs] z_land;
  vector[n_obs_y]                    z_RH;
  vector[n_obs_y]                    z_WS_rain;

  int<lower=1> N_multiplier;
  int<lower=1> grainsize;

  int<lower=1> M;
  array[n_grids_total - 1, M] int<lower=1, upper=n_grids_total> NN_ind;
  matrix[n_grids_total - 1, M]                    NN_dist;
  matrix[n_grids_total - 1, (M * (M - 1) / 2)]   NN_distM;
}


transformed data {
  int max_y = max(y);
  int N_max = max_y + max(250, max_y * N_multiplier);
}


parameters {
  // --- Abundance ---
  real beta0;
  real beta_temp;
  real beta_temp2;
  real beta_rain;
  real beta_rain2;
  vector[n_land_covs] beta_land;

  // --- Detection ---
  real alpha0;
  real alpha_RH;
  real alpha_WS_rain;
  vector[4] alpha_trap_raw;
  real log_sigma_trap;

  // --- Dispersion ---
  real<lower=1e-6> phi;

  // --- NNGP hyperparameters ---
  real<lower=0> alpha_gp;
  real<lower=0> rho_gp;

  // --- Unit-scale spatial field (non-centred) ---
  vector[n_grids_total] w_raw;
}


transformed parameters {
  matrix[n_grids_total, n_dates] lambda_base;

  real sigma_trap = exp(log_sigma_trap);
  vector[5] alpha_trap;
  alpha_trap[1:4] = alpha_trap_raw;
  alpha_trap[5]   = -sum(alpha_trap_raw);

  // Recover the intercept-centred scaled spatial field
  vector[n_grids_total] w_b1 = beta0 + alpha_gp * w_raw;

  for (g in 1:n_grids_total) {
    for (t in 1:n_dates) {
      lambda_base[g, t] = exp(
        w_b1[g] +
        beta_temp  * z_temp[g, t] +
        beta_temp2 * z_temp_sq[g, t] +
        beta_rain  * z_rain[g, t] +
        beta_rain2 * z_rain_sq[g, t] +
        dot_product(beta_land, z_land[g, ])
      );
    }
  }
}


model {
  // --- Abundance priors ---
  beta0      ~ normal(0, 2);
  beta_temp  ~ normal(0, 1);
  beta_temp2 ~ normal(0, 1);
  beta_rain  ~ normal(0, 1);
  beta_rain2 ~ normal(0, 1);
  beta_land  ~ normal(0, 1);

  // --- Detection priors ---
  alpha0         ~ normal(logit(0.10), 0.5);
  alpha_RH       ~ normal(0, 1);
  alpha_WS_rain  ~ normal(0, 1);
  alpha_trap_raw ~ normal(0, sigma_trap);
  log_sigma_trap ~ normal(log(0.5), 0.5);

  // --- Dispersion ---
  phi ~ lognormal(0, 1);

  // --- NNGP hyperparameter priors ---
  alpha_gp ~ normal(0, 1);
  rho_gp   ~ gamma(8, 0.5);

  // --- Non-centred NNGP prior on unit-scale field ---
  w_raw ~ nngp_w(1.0, 1.0 / rho_gp,
                 NN_dist, NN_distM, NN_ind,
                 n_grids_total, M);

  // --- Survey likelihood ---
  target += reduce_sum(
    partial_sum, y, grainsize,
    site_to_grid, date_y, trap_type,
    z_RH, z_WS_rain,
    lambda_base,
    alpha0, alpha_RH, alpha_WS_rain, alpha_trap,
    phi, N_max
  );
}


generated quantities {
  array[n_obs_y] int y_rep;

  vector[n_grids_total] spatial_field = w_b1;

  int  n_near_boundary = 0;
  real max_N_ratio     = 0.0;

  {
    for (k in 1:n_obs_y) {
      int g = site_to_grid[k];
      int t = date_y[k];

      int N_sim    = neg_binomial_2_rng(lambda_base[g, t], phi);
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

  vector[n_obs_y] log_lik_survey;

  for (k in 1:n_obs_y) {
    int g = site_to_grid[k];
    int t = date_y[k];
    real p_trap = inv_logit(alpha0 +
                            alpha_RH      * z_RH[k] +
                            alpha_WS_rain * z_WS_rain[k] +
                            alpha_trap[trap_type[k]]);
    log_lik_survey[k] = binomial_negbin_marginal_lpmf(
      y[k] | p_trap, lambda_base[g, t], phi, N_max);
  }
}

