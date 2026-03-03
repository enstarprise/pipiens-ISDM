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
  model = vgm(psill = var(lambda_df$lambda), 
              model = "Exp", 
              range = max(dist(cbind(lambda_df$x, lambda_df$y))) / 3, 
              nugget = 0)
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
vgm_fit_gauss <- fit.variogram(vg, 
                               model = vgm("Gau"))
vgm_fit_gauss

vgm_fit_sph   <- fit.variogram(vg, model = vgm("Sph"))
vgm_fit_sph

plot(vg, vgm_fit_gauss)
plot(vg, vgm_fit_sph)
plot(vg, vgm_fit_sph, xlim = c(0, max_dist * 1.5))



# ----------- 1/2 VS 3/2 MATERN COVARIANCE FUNCTIONS -------

vg <- variogram(lambda ~ 1, data = lambda_df, locations = ~ x + y)

vgm_exp <- fit.variogram(
  vg,
  model = vgm(
    psill = var(lambda_df$lambda),
    model = "Exp",
    range = max(dist(cbind(lambda_df$x, lambda_df$y))) / 3,
    nugget = 0
  )
)

vgm_m32 <- fit.variogram(
  vg,
  model = vgm(
    psill = var(lambda_df$lambda),
    model = "Mat",
    range = vgm_exp$range[2],  # initialize near Exp
    nugget = 0,
    kappa = 1.5                # ν = 3/2
  )
)


vgm_exp_n <- fit.variogram(
  vg,
  vgm(model = "Exp")
)

vgm_m32_n <- fit.variogram(
  vg,
  vgm(model = "Mat", kappa = 1.5)
)



# Try plotting with explicit x and y
plot(vg$dist, vg$gamma, 
     pch = 16, col = "grey40",
     xlab = "Distance", ylab = "Semivariance",
     main = "Empirical variogram with fitted models")

# Calculate variogram lines
vl_exp <- variogramLine(vgm_exp, maxdist = max(vg$dist))
vl_m32 <- variogramLine(vgm_m32, maxdist = max(vg$dist))

# Add lines to the plot
lines(vl_exp$dist, vl_exp$gamma, col = "red", lwd = 2)
lines(vl_m32$dist, vl_m32$gamma, col = "blue", lwd = 2)

# Add legend
legend(
  "bottomright",
  legend = c("Exponential (ν = 1/2)", "Matérn ν = 3/2"),
  col = c("red", "blue"),
  lwd = 2
)



# Clear any existing plots
graphics.off()

# Try plotting with explicit x and y
plot(vg$dist, vg$gamma, 
     pch = 16, col = "grey40",
     xlab = "Distance", ylab = "Semivariance",
     main = "Empirical variogram with fitted models")




#===========================================================================
# lambda for model with HSGP spatial random effect


spatial_df_full <- read_csv("spatial_df_full.csv")

library(sf)
library(dplyr)

grid <- st_read("01_data/grids/grid_clipped_2point5km.gpkg")
# grid <- read_csv("01_data/grids/grid_clipped_2point5km_wkt.csv")

load("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/AIM 1/pipiens-ISDM/01_data/loaded_data_cov_2point5km.RData")


#  grid_id_lookup named vector to a dataframe
lookup_table<- data.frame(
  grid_id = as.integer(names(grid_id_lookup)),  # The actual grid IDs
  grid_index = as.integer(grid_id_lookup)       # The indices (1, 2, 3, ...)
)


merged_data <- spatial_df_full %>%
  left_join(lookup_table, by = c("grid_id" = "grid_index"))


grid_coords <- grid %>%
  mutate(
    x = st_coordinates(st_centroid(.))[, "X"],
    y = st_coordinates(st_centroid(.))[, "Y"]
  ) %>%
  st_drop_geometry() %>%  # Remove geometry column to get a regular dataframe
  select(grid_id, x, y)  # Adjust column names as needed


final_df <- merged_data %>%
  left_join(grid_coords, by = c("grid_id.y" = "grid_id"))


spatial_df_full <- final_df



vg <- variogram(
  spatial_effect ~ 1,
  data = spatial_df_full,
  locations = ~ x + y
)

plot(vg)

vgm_fit <- fit.variogram(
  vg,
  model = vgm(psill = var(spatial_df_full$spatial_effect), 
              model = "Exp", 
              range = max(dist(cbind(spatial_df_full$x, spatial_df_full$y))) / 3, 
              nugget = 0)
)
vgm_fit
plot(vg, vgm_fit)

plot(vgm_fit)


vgm_fit_gauss <- fit.variogram(vg, 
                               model = vgm("Gau"))
vgm_fit_gauss 
