## ---------------------------------------------------------
## 0. Setup
## ---------------------------------------------------------

library(cmdstanr)
library(posterior)
library(spdep)

fit <- readRDS("~/sharedscratch/expGP_fit_M60.rds")
load("~/pipiens-ISDM/stan/2.5km/stan_data.RData")

set.seed(1)

n_grids <- stan_data$n_grids_total
n_dates <- stan_data$n_dates


## ---------------------------------------------------------
## 1. Extract posterior draws
## ---------------------------------------------------------

draws_basic <- fit$draws(
  c("beta0", "beta_temp", "beta_rain", "beta_land", "lambda_base")
)

draws_df <- as_draws_df(draws_basic)
n_draws  <- ndraws(draws_basic)

beta_land_mat <- as.matrix(fit$draws("beta_land"))      # draws × n_land_covs
lambda_mat    <- as_draws_matrix(fit$draws("lambda_base"))


## ---------------------------------------------------------
## 2. Fixed-effects predictor (one grid, one time)
## ---------------------------------------------------------

log_lambda_hat_fixed <- function(beta0, beta_temp, beta_rain, beta_land,
                                 z_temp_gt, z_rain_gt, z_land_g) {
  beta0 +
    beta_temp * z_temp_gt +
    beta_rain * z_rain_gt +
    sum(beta_land * z_land_g)
}


## ---------------------------------------------------------
## 3. Compute latent sresiduals
## ---------------------------------------------------------

residual_array <- array(
  NA_real_,
  dim = c(n_draws, n_grids, n_dates)
)

for (s in 1:n_draws) {
  
  beta0_s      <- draws_df$beta0[s]
  beta_temp_s <- draws_df$beta_temp[s]
  beta_rain_s <- draws_df$beta_rain[s]
  beta_land_s <- beta_land_mat[s, ]
  
  for (g in 1:n_grids) {
    for (t in 1:n_dates) {
      
      log_lambda_fixed <-
        log_lambda_hat_fixed(
          beta0_s,
          beta_temp_s,
          beta_rain_s,
          beta_land_s,
          stan_data$z_temp[g, t],
          stan_data$z_rain[g, t],
          stan_data$z_land[g, ]
        )
      
      colname <- sprintf("lambda_base[%d,%d]", g, t)
      log_lambda_full <- log(lambda_mat[s, colname])
      
      residual_array[s, g, t] <- log_lambda_full - log_lambda_fixed
    }
  }
}


## ---------------------------------------------------------
## 4. Remove posterior-mean GP effect
##    (this isolates *unexplained* spatial structure)
## ---------------------------------------------------------

# Posterior mean spatial effect per grid (averaged over time)
gp_mean <- apply(residual_array, c(2, 3), mean)   # grids × dates
gp_mean <- rowMeans(gp_mean)                      # grids

# True misspecification residuals
true_residuals <- array(NA_real_, dim = dim(residual_array))

for (s in 1:n_draws) {
  for (g in 1:n_grids) {
    true_residuals[s, g, ] <- residual_array[s, g, ] - gp_mean[g]
  }
}

# Aggregate to grid level
true_residual_grid <- apply(true_residuals, c(1, 2), mean)
# dimensions: draws × grids


## ---------------------------------------------------------
## 5. Spatial weights
## ---------------------------------------------------------

coords <- as.matrix(stan_data$grid_coords)

# Distance-based neighbors (choose threshold carefully!)
nb <- dnearneigh(coords, 0, 50)  # e.g. 50 km
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)


## ---------------------------------------------------------
## 6. Moran’s I per posterior draw
## ---------------------------------------------------------

moran_I <- numeric(n_draws)

for (s in 1:n_draws) {
  moran_I[s] <- moran(
    true_residual_grid[s, ],
    lw,
    zero.policy = TRUE
  )$I
}


## ---------------------------------------------------------
## 7. Summaries
## ---------------------------------------------------------

cat("\nPosterior Moran's I summary:\n")
print(quantile(moran_I, c(0.025, 0.5, 0.975)))

cat("\nP(Moran's I > 0):\n")
print(mean(moran_I > 0))


## ---------------------------------------------------------
## 8. Optional diagnostic plots
## ---------------------------------------------------------

hist(moran_I, breaks = 40,
     main = "Posterior distribution of Moran's I\n(latent-scale residuals)",
     xlab = "Moran's I")
abline(v = 0, col = "red", lwd = 2)
