// =============================================================================
// LATENT NNGP — SPATIAL FIELD CENTRED AT THE INTERCEPT
//
// This is a variant of model_nngp.stan that addresses the well-known mixing
// problem in latent NNGP models: the intercept (beta0) and the spatial field
// (w_gp) are highly correlated in the standard parameterisation, causing slow
// convergence and poor chain mixing.
//
// THE FIX:
//   Instead of sampling w_gp ~ GP(0, C*) and adding it to beta0 separately,
//   we sample w_b1 ~ GP(beta0, C*) — a spatial field centred at the intercept.
//   Internally the GP prior subtracts the intercept before evaluating the
//   log-density: w = w_b1 - beta0, then w ~ GP(0, C*).
//
//   In the linear predictor, w_b1 replaces (beta0 + w_gp):
//     log_lambda = w_b1[g] + beta_temp*... + beta_rain*... + beta_land*...
//   so beta0 is absorbed into w_b1 entirely and removed as a separate parameter.
//
// REFERENCE:
//   Zhang (2018) https://mc-stan.org/learn-stan/case-studies/nngp.html
//   (latent NNGP with intercept-centred spatial process)
// =============================================================================


functions {

  // ---------------------------------------------------------------------------
  // NNGP LOG-PRIOR: SPATIAL FIELD CENTRED AT THE INTERCEPT
  //
  // Identical to nngp_w_lpdf() in model_nngp.stan, except:
  //   - Takes w_b1 (the centred field) and intercept as inputs
  //   - Internally computes w = w_b1 - intercept before the GP evaluation
  //   - The GP prior is on the zero-mean residual w, not on w_b1 directly
  //
  // This means w_b1 ~ GP(intercept, C*), i.e. the field fluctuates around
  // the intercept rather than around zero.
  // ---------------------------------------------------------------------------
  real nngp_w_lpdf(vector w_b1, real sigmasq, real phi,
                   matrix NN_dist, matrix NN_distM,
                   array[,] int NN_ind, int N, int M,
                   real intercept) {

    vector[N] w    = w_b1 - intercept;  // zero-mean residual field
    vector[N] I_Aw = w;
    vector[N] V;
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
        for (j in 1:dim) iNNdistM[j, j] = 1.0;
      }

      matrix[dim, dim] iNNCholL = cholesky_decompose(iNNdistM);
      iNNcorr = to_vector(exp(-phi * NN_dist[(i - 1), 1:dim]));

      v     = mdivide_left_tri_low(iNNCholL, iNNcorr);
      V[i]  = 1.0 - dot_self(v);
      v2    = mdivide_right_tri_low(v', iNNCholL);

      // Residual using zero-mean w (not w_b1) for the neighbour values
      I_Aw[i] -= v2 * w[NN_ind[(i - 1), 1:dim]];
    }
    V[1] = 1.0;

    return -0.5 * (
      (1.0 / sigmasq) * dot_product(I_Aw, (I_Aw ./ V)) +
      sum(log(V)) +
      N * log(sigmasq)
    );
  }

  // ---------------------------------------------------------------------------
  // EXISTING: MARGINALISED BINOMIAL-NEGBIN LIKELIHOOD (unchanged)
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
  // Identical to model_nngp.stan — no changes to the data block
  int<lower=1> n_grids_total;
  int<lower=1> n_grids_obs;
  int<lower=1> n_dates;
  int<lower=1> n_obs_y;
  int<lower=1> n_obs_po;
  int<lower=1> n_land_covs;

  array[n_grids_obs] int<lower=1, upper=n_grids_total> obs_grid_idx;
  array[n_obs_y]     int<lower=1, upper=n_grids_total> site_to_grid;
  array[n_obs_y]     int<lower=1, upper=n_dates>       date_y;
  array[n_obs_y]     int<lower=1, upper=5>             trap_type;
  array[n_obs_po]    int<lower=1, upper=n_grids_total> po_grid_idx;
  array[n_obs_po]    int<lower=1, upper=n_dates>       date_po;

  array[n_obs_y]  int<lower=0> y;
  array[n_obs_po] int<lower=1> ones;

  matrix[n_grids_total, n_dates]     z_temp;
  matrix[n_grids_total, n_dates]     z_rain;
  matrix[n_grids_total, n_land_covs] z_land;
  vector[n_grids_total]              z_poi;
  vector[n_grids_total]              z_reports;
  vector[n_obs_y]                    z_RH;
  vector[n_obs_y]                    z_WS_rain;
  matrix[n_grids_total, n_dates] z_temp_sq;
  matrix[n_grids_total, n_dates] z_rain_sq;

  vector[n_grids_total] area_grid;
  real  CONSTANT;
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
  // --- Abundance model ---
  // NOTE: beta0 is REMOVED as a separate parameter.
  // It is now absorbed into w_b1 (the intercept-centred spatial field).
  // beta0 still has a prior (see model block) and is passed to nngp_w_lpdf
  // as the `intercept` argument, but it is no longer in the linear predictor
  // directly — w_b1[g] already contains the intercept contribution.
  real beta0;           // intercept: sampled for its prior, absorbed into w_b1
  real beta_temp;
  real beta_rain;
  real beta_temp2;
  real beta_rain2;
  vector[n_land_covs] beta_land;

  // --- Detection model (unchanged) ---
  real alpha0;
  real alpha_RH;
  real alpha_WS_rain;
  vector[4] alpha_trap_raw;
  real log_sigma_trap;

  // --- Thinning model (unchanged) ---
  real delta0;
  real delta_poi;
  real delta_reports;

  // --- Dispersion (unchanged) ---
  real<lower=1e-6> phi;

  // --- NNGP hyperparameters (unchanged) ---
  real<lower=0> alpha_gp;
  real<lower=0> rho_gp;

  // --- Intercept-centred spatial field ---
  // w_b1[g] = beta0 + w_gp[g], where w_gp is the zero-mean spatial residual.
  // Sampling w_b1 jointly with beta0 (rather than w_gp separately) breaks
  // the strong beta0 <-> w_gp posterior correlation.
  vector[n_grids_total] w_b1;
}


transformed parameters {
  matrix[n_grids_total, n_dates] lambda_base;
  matrix[n_grids_total, n_dates] lambda_thinned;
  vector[n_dates] background;

  real sigma_trap = exp(log_sigma_trap);
  vector[5] alpha_trap;

  alpha_trap[1:4] = alpha_trap_raw;
  alpha_trap[5]   = -sum(alpha_trap_raw);

  for (g in 1:n_grids_total) {
    real p_thin = inv_logit(delta0 +
                            delta_poi     * z_poi[g] +
                            delta_reports * z_reports[g]);

    for (t in 1:n_dates) {
      // w_b1[g] replaces (beta0 + w_gp[g]) from model_nngp.stan.
      // The intercept is now baked into the spatial field itself.
      real log_lambda_base = w_b1[g] +
                       beta_temp  * z_temp[g, t] +
                       beta_temp2 * z_temp_sq[g, t] + 
                       beta_rain  * z_rain[g, t] +
                       beta_rain2 * z_rain_sq[g, t] +   
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
  // --- Abundance priors ---
  // beta0 still gets a prior — it anchors the mean of the spatial field.
  // Without this prior, beta0 and the mean of w_b1 would be unidentified.
  beta0     ~ normal(0, 2);
  beta_temp ~ normal(0, 1);
  beta_rain ~ normal(0, 1);
  beta_temp2 ~ normal(0, 1);
  beta_rain2 ~ normal(0, 1);
  beta_land ~ normal(0, 0.2);

  // --- Detection priors (unchanged) ---
  alpha0        ~ normal(logit(0.10), 0.5);
  alpha_RH      ~ normal(0, 1);
  alpha_WS_rain ~ normal(0, 1);
  alpha_trap_raw ~ normal(0, sigma_trap);
  log_sigma_trap ~ normal(log(0.5), 0.5);

  // --- Thinning priors (unchanged) ---
  delta0        ~ normal(-1, 1);
  delta_poi     ~ normal(0, 1);
  delta_reports ~ normal(0, 1);

  // --- Dispersion prior (unchanged) ---
  phi ~ lognormal(0, 1);

  // --- NNGP hyperparameter priors (unchanged) ---
  alpha_gp ~ normal(0, 1); 
  rho_gp   ~ gamma(8, 0.5);

  // --- Intercept-centred NNGP prior on w_b1 ---
  // The key change: beta0 is passed as the `intercept` argument.
  // Inside nngp_w_lpdf, w = w_b1 - beta0 is computed, and the GP prior
  // is evaluated on that zero-mean residual. This means w_b1 fluctuates
  // around beta0 rather than around 0, which decorrelates beta0 from w_b1
  // and substantially improves HMC mixing.
  w_b1 ~ nngp_w(square(alpha_gp), 1.0 / rho_gp,
                NN_dist, NN_distM, NN_ind,
                n_grids_total, M,
                beta0);   // <-- the only difference from model_nngp.stan

  // --- Survey likelihood (unchanged) ---
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

  // --- Presence-only likelihood (unchanged) ---
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

  // Spatial field for visualisation.
  // w_b1 is the intercept-centred field; subtract beta0 to recover the
  // zero-mean spatial residual if you want to compare with model_nngp.stan:
  //   spatial_residual = w_b1 - beta0
  vector[n_grids_total] spatial_field = w_b1;

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
  // LOG-LIKELIHOOD FOR LOO-CV (unchanged)
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
