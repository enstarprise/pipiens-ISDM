
library(readr)
library(ggplot2)
library(patchwork)

load("01_data/processedCovariates/2.5km/stan_data_2.5km_HSGP_all_both.RData")
draws <- readRDS("03_posterior/savedObjects/nngp_M10_temprain_draws.rds")
temperature <- read_csv("01_data/covariates/mean_temp_21day_2point5km_new.csv")
rainfall <- read_csv("01_data/covariates/rainfall_28day_2point5km_new.csv")



# GET SCALING PARAMETERS FROM STAN DATA DIRECTLY

# Extract z-score vectors from Stan data matrix
z_rain_vals <- as.vector(stan_data$z_rain)
z_temp_vals <- as.vector(stan_data$z_temp)

#  mean and SD from the actual z-scored values used in the model
rain_mean_z <- mean(z_rain_vals)
rain_sd_z   <- sd(z_rain_vals)

temp_mean_z <- mean(z_temp_vals)
temp_sd_z   <- sd(z_temp_vals)


# GENERATE PREDICTION SEQUENCES IN Z-SCORE SPACE


rain_seq_z <- seq(
  quantile(z_rain_vals, 0.01),
  quantile(z_rain_vals, 0.99),
  length.out = 200
)

temp_seq_z <- seq(
  quantile(z_temp_vals, 0.01),
  quantile(z_temp_vals, 0.99),
  length.out = 200
)

# BACK TRANSFORM TO ORIGINAL SCALE FOR X AXIS
#  original raw variable means and SDs for labelling
rain_mean_orig <- mean(rainfall$rainfall_28d, na.rm = TRUE)
rain_sd_orig   <- sd(rainfall$rainfall_28d, na.rm = TRUE)

temp_mean_orig <- mean(temperature$mean_temp_21d_celsius, na.rm = TRUE)
temp_sd_orig   <- sd(temperature$mean_temp_21d_celsius, na.rm = TRUE)

rain_seq_orig <- rain_seq_z * rain_sd_orig + rain_mean_orig
temp_seq_orig <- temp_seq_z * temp_sd_orig + temp_mean_orig


# COMPUTE POSTERIOR PREDICTED VALUES
compute_marginal <- function(draws, beta_name, beta2_name, z_seq) {
  
  beta  <- draws[[beta_name]]
  beta2 <- draws[[beta2_name]]
  
  pred_matrix <- exp(outer(beta, z_seq) + outer(beta2, z_seq^2))
  
  data.frame(
    z       = z_seq,
    mean    = apply(pred_matrix, 2, mean),
    lower   = apply(pred_matrix, 2, quantile, 0.025),
    upper   = apply(pred_matrix, 2, quantile, 0.975),
    lower50 = apply(pred_matrix, 2, quantile, 0.25),
    upper50 = apply(pred_matrix, 2, quantile, 0.75)
  )
}

rain_pred <- compute_marginal(draws, "beta_rain", "beta_rain2", rain_seq_z)
temp_pred <- compute_marginal(draws, "beta_temp", "beta_temp2", temp_seq_z)

# Add back-transformed x axis
rain_pred$x_orig <- rain_seq_orig
temp_pred$x_orig <- temp_seq_orig


# PLOT

plot_marginal <- function(pred_df, x_label, colour = "#2166ac") {
  ggplot(pred_df, aes(x = x_orig)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                fill = colour, alpha = 0.05) +
    geom_ribbon(aes(ymin = lower50, ymax = upper50),
                fill = colour, alpha = 0.10) +
    geom_line(aes(y = mean),
              colour = colour, linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    labs(
      x = x_label,
      y = "Marginal effect on intensity"
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.title  = element_text(size = 12),
      axis.text   = element_text(size = 10),
      plot.margin = margin(10, 15, 10, 10)
    )
}

p_rain <- plot_marginal(rain_pred,
                        x_label = "28-day rainfall (mm)",
                        colour  = "#2166ac")

p_temp <- plot_marginal(temp_pred,
                        x_label = "Temperature (°C)",
                        colour  = "#d6604d")

p_combined <- p_rain + p_temp +
  plot_annotation(
    title = expression(italic("Cx. pipiens") ~ "intensity — marginal climate effects"),
    theme = theme(plot.title = element_text(size = 13, hjust = 0.5))
  )

print(p_combined)

ggsave("marginal_climate_effects.png",
       plot  = p_combined,
       width = 10,
       height = 5,
       dpi   = 300)
