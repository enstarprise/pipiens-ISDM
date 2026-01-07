################################################################################
# Cx pipiens in Scotland - Variogram
# Author: shennice
# Date: December 2025

# This is the posterior checking of lambda variogram to create a model
# with spatial dependencies
################################################################################

library(readr)
library(gstat)

# this is just the lambda value at the first time point
lambda_df <- read_csv("02_model/lambda_variogram_df.csv")

vg <- variogram(
  lambda ~ 1,
  data = lambda_df,
  locations = ~ x + y
)

plot(vg)


#> there cannot be trends in the relationship between the 
#> coords and the response:
par(mfrow = c(1,2))
plot(combined_data$longitude,combined_data$Culex)
plot(combined_data$latitude,combined_data$Culex)

lm_long <- lm(lambda_df$lambda ~ lambda_df$x)
summary(lm_long)

lm_lat <- lm(lambda_df$lambda ~ lambda_df$y)
summary(lm_lat)


# ===========================================================================
# --- Fit and empirical variogram

vg <- variogram(lambda ~ 1, data = lambda_df, locations = ~ x + y)

# --- Fit an exponential variogram model
vgm_fit <- fit.variogram(
  vg,
  model = vgm(psill = var(lambda_df$lambda), model = "Exp", range = max(dist(cbind(lambda_df$x, lambda_df$y))) / 3, nugget = 0)
)
vgm_fit
plot(vg, vgm_fit)

max_dist <- max(vg$dist)
plot(vg, vgm_fit, xlim = c(0, max_dist * 1.5))

# Create a smooth distance sequence up to 2 × the max empirical distance
d_seq <- seq(0, max_dist * 2, length.out = 200)

# Predict model semivariance
gamma_pred <- variogramLine(vgm_fit, maxdist = max_dist * 2)$gamma

# Plot
plot(vg$dist, vg$gamma,
     pch = 1,
     xlab = "distance",
     ylab = "semivariance",
     xlim = c(0, max_dist * 2)
)

lines(d_seq, gamma_pred, col = "blue", lwd = 2)


max_dist <- max(vg$dist)




# --- Fit other variogram shapes
vgm_fit_gauss <- fit.variogram(vg, model = vgm("Gau"))
vgm_fit_sph   <- fit.variogram(vg, model = vgm("Sph"))

plot(vg, vgm_fit_gauss)
plot(vg, vgm_fit_sph)
plot(vg, vgm_fit_sph, xlim = c(0, max_dist * 1.5))

