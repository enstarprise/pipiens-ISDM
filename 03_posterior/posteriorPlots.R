
library(dplyr)
library(bayesplot)
library(ggplot2)
library(posterior)
library(viridis)
library(tidyverse)
library(Cairo)
library(readr)

posterior_mat <- readRDS("03_posterior/savedObjects/posterior_mat_nngp_M10.rds")

posterior_mat_detection <- readRDS("03_posterior/savedObjects/posterior_mat_detection_nngp.rds")

# ============================================================================
# POSTEIROR VISUALISATION
# ============================================================================

main_params <- c(
  "beta0", "beta_temp", "beta_temp2", "beta_rain", "beta_rain2",
  "beta_land[1]", "beta_land[2]", "beta_land[3]",
  "beta_land[4]", "beta_land[5]", "beta_land[6]",
  "beta_land[7]", "beta_land[8]",
  "alpha0", "alpha_RH", "alpha_WS_rain", 
  "sigma_trap", "alpha_trap[1]", "alpha_trap[1]",
  "alpha_trap[2]", "alpha_trap[3]", "alpha_trap[4]", "alpha_trap[5]",
   # "delta0", "delta_poi", "delta_reports",
  "phi"
)


main_params <- c(
  "beta_land[1]", "beta_land[2]", "beta_land[3]",
  "beta_land[4]", "beta_land[5]", "beta_land[6]", 
  "beta_land[7]", "beta_land[8]",
  "beta_rain", "beta_rain2",
  "beta_temp", "beta_temp2" 
)

param_names <- c(
  "beta_land[1]" = "Wetland",
  "beta_land[2]" = "Woodland",
  "beta_land[3]" = "Freshwater",
  "beta_land[4]" = "Grassland",
  "beta_land[5]" = "Arable",
  "beta_land[6]" = "Livestock\nDensity",
  "beta_land[7]" = "Elevation",
  "beta_land[8]" = "Buildings",
  "beta_temp" = "Temperature\n(3 week avg)",
  "beta_temp2" = "Temperature²",
  "beta_rain" = "Rainfall\n(4 week\ncumulative)",
  "beta_rain2" = "Rainfall²"
)

param_labels <- param_names[main_params]

color_scheme_set("darkgray")

intervals_plot_1 <- mcmc_intervals(
  posterior_mat,
  pars = main_params,
  prob = 0.9,
  prob_outer = 0.95,
  point_est = "mean"
) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  scale_y_discrete(labels = param_labels) +
  scale_x_continuous(
    breaks = seq(-2, 4, by = 0.5)  # adjust range and interval as needed
  ) +
  labs(x = "Posterior coefficient estimates", y = "Covariates") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),  # ensure x-axis labels are readable
    legend.position = "none"
  )

intervals_plot_1

ggsave(filename = "03_posterior/priorPosteriorPlots/lambda_params_survey_nngp.png",
       plot = intervals_plot_1,
       device = "png",
       width = 7,
       height = 8,
       units = "in",
       dpi = 300,
       background = "transparent")




# =======================================================================

lambda_params <- c(
  "beta_land[1]", "beta_land[2]", "beta_land[3]",
  "beta_land[4]", "beta_land[5]", "beta_land[6]", 
  "beta_land[7]", "beta_land[8]",
  "beta_rain", "beta_rain2",
  "beta_temp","beta_temp2" 
)

param_names <- c(
  "beta_land[1]" = "Wetland",
  "beta_land[2]" = "Woodland",
  "beta_land[3]" = "Freshwater",
  "beta_land[4]" = "Grassland",
  "beta_land[5]" = "Arable",
  "beta_land[6]" = "Livestock\nDensity",
  "beta_land[7]" = "Elevation",
  "beta_land[8]" = "Buildings",
  "beta_temp" = "Temperature\n(3 week avg)",
  "beta_temp2" = "Temperature²",
  "beta_rain" = "Rainfall\n(4 week\ncumulative)",
  "beta_rain2" = "Rainfall²"
)



# Calculate summary statistics
posterior_summary <- as.data.frame(posterior_mat) %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  group_by(parameter) %>%
  summarise(
    median = median(value),
    lower_80 = quantile(value, 0.10),
    upper_80 = quantile(value, 0.90),
    lower_90 = quantile(value, 0.05),   # add this
    upper_90 = quantile(value, 0.95),   # add this
    lower_95 = quantile(value, 0.025),
    upper_95 = quantile(value, 0.975)
  ) %>%
  mutate(
    parameter = factor(parameter, levels = lambda_params),
    order = match(parameter, lambda_params)
  )

# Create plot with ALL BLACK color and SQUARE shapes
intervals_plot_2 <- ggplot(posterior_summary, 
                           aes(y = reorder(parameter, order),
                               x = median)) +  # No color or shape mapping needed
  
  # 95% CI - Thin black
  geom_errorbar(aes(xmin = lower_95, xmax = upper_95),
                width = 0.15,
                size = 0.3,
                color = "black",  # Black color
                alpha = 0.4,
                orientation = "y") +
  
  # 80% CI - Thinner black  
  geom_errorbar(aes(xmin = lower_90, xmax = upper_90),
                width = 0,
                size = 0.8,
                color = "black",  # Black color
                orientation = "y") +
  
  # Points with SQUARE shape (shape = 15) and BLACK color
  geom_point(size = 3.5, 
             shape = 15,  # 15 = filled square
             color = "black") +  # Black color
  
  # Dashed vertical reference line at 0
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             color = "gray50", 
             size = 0.5,
             alpha = 0.7) +
  
  # Apply parameter names as y-axis labels
  scale_y_discrete(labels = param_names) +
  
  # Add more x-axis labels
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 8),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  
  # Add axis titles
  labs(
    y = "Parameters",
    x = "Posterior Estimate",
    #title = "Ecological Parameter Estimates with 80% and 95% Credible Intervals"
  ) +
  
  # Use theme_bw() as in your original
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 12, face = "bold", 
                                margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", 
                                margin = margin(r = 10)),
    plot.title = element_text(size = 14, face = "bold", 
                              hjust = 0.5, margin = margin(b = 15)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    # Remove legend completely since all are black squares
    legend.position = "none"
  )

print(intervals_plot_2)

ggsave(filename = "03_posterior/priorPosteriorPlots/key_params_nngp_M10.png",
       plot = intervals_plot_2,
       device = "png",
       width = 7,
       height = 8,
       units = "in",
       dpi = 300)

#================================================

detection_params <- c(
  "alpha_RH", "alpha_WS_rain", "alpha_trap[1]", "alpha_trap[2]",
  "alpha_trap[3]", "alpha_trap[4]", "alpha_trap[5]",
  "delta_poi", "delta_reports"
)

detection_names <- c(
  "alpha_RH" = "Relative\nHumidity", 
  "alpha_WS_rain" = "Rainfall", 
  "alpha_trap[1]" = "BGP", 
  "alpha_trap[2]" = "BGS",
  "alpha_trap[3]" = "GT", 
  "alpha_trap[4]" = "LT", 
  "alpha_trap[5]" = "MM", 
  "delta_poi" = "Points of\nInterest", 
  "delta_reports" = "Reporting\nIntensity"
)



intervals_plot_2 <- mcmc_intervals(
  posterior_mat_detection,
  pars = detection_params,
  prob = 0.9,
  prob_outer = 0.95,
  point_est = "mean"
) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  scale_y_discrete(labels = detection_names) +
  scale_x_continuous(
    breaks = seq(-6, 6, by = 2)  # adjust range and interval as needed
  ) +
  labs(x = "Posterior coefficient estimates", y = "Covariates") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),  # ensure x-axis labels are readable
    legend.position = "none"
  )

intervals_plot_2








######
# posterior_mat_detection <- fit$draws(detection_params) %>% 
#   posterior::as_draws_matrix()

# saveRDS(posterior_mat_detection, "~/pipiens-ISDM/stan/2.5km/posterior_mat_detection.rds")

# Define colors and shapes for each parameter individually
param_styles <- data.frame(
  parameter = detection_params,
  color = c(rep("#2E86AB", 7), rep("#A23B72", 2)),  # First 7 survey params, last 2 CS params
  shape = c(rep(16, 7), rep(17, 2))  # 16 = circle, 17 = triangle
)

posterior_summary <- as.data.frame(posterior_mat_detection) %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  group_by(parameter) %>%
  summarise(
    median = median(value),
    lower_80 = quantile(value, 0.10),
    upper_80 = quantile(value, 0.90),
    lower_90 = quantile(value, 0.05),   # add this
    upper_90 = quantile(value, 0.95),   # add this
    lower_95 = quantile(value, 0.025),
    upper_95 = quantile(value, 0.975)
  ) %>%
  left_join(param_styles, by = "parameter") %>%
  mutate(
    parameter = factor(parameter, levels = detection_params),
    order = match(parameter, detection_params)
  )

# Create plot like the first one but with color-coded symbols
intervals_plot_3_thin <- ggplot(posterior_summary, 
                                aes(y = reorder(parameter, order),
                                    x = median)) +  # No aesthetic mapping
  
  # 95% CI - Thin black
  geom_errorbar(aes(xmin = lower_95, xmax = upper_95),
                width = 0.15,
                size = 0.3,
                color = "black",
                alpha = 0.4,
                orientation = "y") +
  
  # 80% CI - Thinner black  
  geom_errorbar(aes(xmin = lower_90, xmax = upper_95),
                width = 0,
                size = 0.8,
                color = "black",
                orientation = "y") +
  
  # Points with individual colors and shapes
  geom_point(aes(color = parameter, shape = parameter),  # Map both to parameter
             size = 3.5) +
  
  # Dashed vertical reference line at 0
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             color = "gray50", 
             size = 0.5,
             alpha = 0.7) +
  
  # Apply parameter names as y-axis labels
  scale_y_discrete(labels = detection_names) +
  
  # Manual color and shape scales
  scale_color_manual(
    name = "Parameter Type",
    values = setNames(param_styles$color, param_styles$parameter),
    labels = c(rep("Survey Parameters", 7), rep("CS Parameters", 2))
  ) +
  scale_shape_manual(
    name = "Parameter Type",
    values = setNames(param_styles$shape, param_styles$parameter),
    labels = c(rep("Survey Parameters", 7), rep("CS Parameters", 2))
  ) +
  
  # Add more x-axis labels
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 8),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  
  # Add axis titles
  labs(
    y = "Parameters",
    x = "Posterior Estimate"
  ) +
  # Use theme_minimal() with customizations
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 10, face = "bold", 
                                margin = margin(t = 10)),
    axis.title.y = element_text(size = 10, face = "bold", 
                                margin = margin(r = 10)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none"
  ) 

print(intervals_plot_3_thin)

ggsave(filename = "03_posterior/priorPosteriorPlots/detection_params_nngp.png",
       plot = intervals_plot_3_thin,
       device = "png",
       width = 5,
       height = 6,
       units = "in",
       dpi = 300)






####### RAINFALL AND TEMP NON LINEAR PLOTS
library(cmdstanr)
library(posterior)
library(ggplot2)
library(dplyr)
library(patchwork)

load("01_data/processedCovariates/2.5km/stan_data_2.5km_HSGP_all_both.RData")
draws <- readRDS("03_posterior/savedObjects/hsgp_temprain_draws.rds")

# temperature <- read_csv("01_data/covariates/mean_temp_21day_2point5km_new.csv")
# rainfall <- read_csv("01_data/covariates/rainfall_28day_2point5km_new.csv")


temperature[is.na(temperature)] <- 0
temp_mean <- mean(temperature$mean_temp_21d_celsius, na.rm = TRUE)
temp_sd   <- sd(temperature$mean_temp_21d_celsius, na.rm = TRUE)

temp_seq_orig <- seq(
  min(temperature$mean_temp_21d_celsius, na.rm = TRUE),
  max(temperature$mean_temp_21d_celsius, na.rm = TRUE),
  length.out = 200
)

rainfall[is.na(rainfall)] <- 0
rain_mean <- mean(rainfall$rainfall_28d, na.rm = TRUE)
rain_sd   <- 150
rain_sd   <- sd(rainfall$rainfall_28d, na.rm = TRUE)

# rain_seq_orig <- seq(
#   quantile(rainfall$rainfall_28d, 0.01, na.rm = TRUE),
#   quantile(rainfall$rainfall_28d, 0.95, na.rm = TRUE),
#   length.out = 200
# )

rain_seq_orig <- seq(
  min(rainfall$rainfall_28d, na.rm = TRUE),
  max(rainfall$rainfall_28d, na.rm = TRUE),
  length.out = 200
)

# Work purely in z-score space — no back transformation
rain_seq_z <- seq(
  quantile(as.vector(stan_data$z_rain), 0.01),
  quantile(as.vector(stan_data$z_rain), 0.99),
  length.out = 200
)

# Use z-scores directly as x_orig for now
rain_pred <- compute_marginal(draws, "beta_rain", "beta_rain2", rain_seq_z)
rain_pred$x_orig <- rain_seq_z



# Back-transform to z-score scale for prediction
rain_seq_z <- (rain_seq_orig - rain_mean) / rain_sd
temp_seq_z <- (temp_seq_orig - temp_mean) / temp_sd


compute_marginal <- function(draws, beta_name, beta2_name, z_seq) {
  
  beta  <- draws[[beta_name]]
  beta2 <- draws[[beta2_name]]
  
  # Matrix: rows = draws, cols = z values
  pred_matrix <- outer(beta, z_seq) + outer(beta2, z_seq^2)
  
  # Summarise across draws
  data.frame(
    z        = z_seq,
    mean     = apply(pred_matrix, 2, mean),
    lower    = apply(pred_matrix, 2, quantile, 0.025),
    upper    = apply(pred_matrix, 2, quantile, 0.975),
    lower50  = apply(pred_matrix, 2, quantile, 0.25),
    upper50  = apply(pred_matrix, 2, quantile, 0.75)
  )
}

rain_pred <- compute_marginal(draws, "beta_rain", "beta_rain2", rain_seq_z)
temp_pred <- compute_marginal(draws, "beta_temp", "beta_temp2", temp_seq_z)

# Add original scale x axis
rain_pred$x_orig <- rain_seq_orig
temp_pred$x_orig <- temp_seq_orig


plot_marginal <- function(pred_df, x_label, colour = "#2166ac") {
  ggplot(pred_df, aes(x = x_orig)) +
    #geom_ribbon(aes(ymin = lower, ymax = upper),
     #           fill = colour, alpha = 0.15) +
    geom_ribbon(aes(ymin = lower50, ymax = upper50),
                fill = colour, alpha = 0.1) +
    geom_line(aes(y = mean),
              colour = colour, linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    labs(
      x = x_label,
      y = "Marginal effect on log intensity"
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.title  = element_text(size = 12),
      axis.text   = element_text(size = 10),
      plot.margin = margin(10, 15, 10, 10)
    )
}

p_rain <- plot_marginal(rain_pred,
                        x_label  = "28-day rainfall (mm)",
                        colour   = "#2166ac")

p_temp <- plot_marginal(temp_pred,
                        x_label  = "Temperature (°C)",
                        colour   = "#d6604d")

# Combine with patchwork
library(patchwork)
p_combined <- p_rain + p_temp +
  plot_annotation(
    title   = expression(italic("Cx. pipiens") ~ "intensity ~ marginal climate effects"),
    theme   = theme(plot.title = element_text(size = 13, hjust = 0.5))
  )

print(p_combined)

ggsave("marginal_climate_effects.png",
       plot   = p_combined,
       width  = 10,
       height = 5,
       dpi    = 300)


## MARGINAL EFFECT ON INTENSITY (NOT LOGGED)
compute_marginal <- function(draws, beta_name, beta2_name, z_seq) {
  
  beta  <- draws[[beta_name]]
  beta2 <- draws[[beta2_name]]
  
  # Matrix: rows = draws, cols = z values
  pred_matrix <- exp(outer(beta, z_seq) + outer(beta2, z_seq^2))
  
  # Summarise across draws
  data.frame(
    z        = z_seq,
    mean     = apply(pred_matrix, 2, mean),
    lower    = apply(pred_matrix, 2, quantile, 0.025),
    upper    = apply(pred_matrix, 2, quantile, 0.975),
    lower50  = apply(pred_matrix, 2, quantile, 0.25),
    upper50  = apply(pred_matrix, 2, quantile, 0.75)
  )
}

rain_pred <- compute_marginal(draws, "beta_rain", "beta_rain2", rain_seq_z)
temp_pred <- compute_marginal(draws, "beta_temp", "beta_temp2", temp_seq_z)

# Add original scale x axis
rain_pred$x_orig <- rain_seq_orig
temp_pred$x_orig <- temp_seq_orig


plot_marginal <- function(pred_df, x_label, colour = "#2166ac") {
  ggplot(pred_df, aes(x = x_orig)) +
    #geom_ribbon(aes(ymin = lower, ymax = upper),
     #           fill = colour, alpha = 0.15) +
    geom_ribbon(aes(ymin = lower50, ymax = upper50),
                fill = colour, alpha = 0.30) +
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
                        x_label  = "28-day rainfall (mm)",
                        colour   = "#2166ac")

p_temp <- plot_marginal(temp_pred,
                        x_label  = "Temperature (°C)",
                        colour   = "#d6604d")

library(patchwork)
p_combined <- p_rain + p_temp +
  plot_annotation(
    title   = expression(italic("Cx. pipiens") ~ "intensity ~ marginal climate effects"),
    theme   = theme(plot.title = element_text(size = 13, hjust = 0.5))
  )

print(p_combined)



######### SURVEY ONLY MODEL COMPARISON #############

library(ggplot2)
library(tidybayes)
library(dplyr)

isdm_posterior_mat <- readRDS("03_posterior/savedObjects/posterior_mat_nngp.rds")
isdm_posterior_mat <- as.data.frame(isdm_posterior_mat)

posterior_mat <- readRDS("02_model/stan/modelRuns/model_survey_nngp.rds")
posterior_mat <- as.data.frame(posterior_mat)


# Prepare data for plotting
model1_long <- posterior_mat %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  mutate(model = "Model 1")

model2_long <- isdm_posterior_mat %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  mutate(model = "Model 2")

combined_post <- bind_rows(model1_long, model2_long)

# Create overlaid density plots
ggplot(combined_post, aes(x = value, fill = model)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal() +
  labs(title = "Posterior Distributions Comparison",
       x = "Coefficient Value", y = "Density") +
  scale_fill_manual(values = c("steelblue", "coral"))




param_names <- c(
  "beta_land[1]" = "Wetland",
  "beta_land[2]" = "Woodland",
  "beta_land[3]" = "Freshwater",
  "beta_land[4]" = "Grassland",
  "beta_land[5]" = "Arable",
  "beta_land[6]" = "Livestock\nDensity",
  "beta_land[7]" = "Elevation",
  "beta_land[8]" = "Buildings",
  "beta_temp" = "Temperature\n(3 week avg)",
  "beta_temp2" = "Temperature²",
  "beta_rain" = "Rainfall\n(4 week\ncumulative)",
  "beta_rain2" = "Rainfall²"
)

library(dplyr)
library(tidyr)
library(ggplot2)
library(tidybayes)

# Rename columns in your posterior matrices
rename_posterior <- function(post_matrix, name_mapping) {
  # Find which columns exist in your data
  existing_cols <- intersect(names(post_matrix), names(name_mapping))
  
  # Rename them
  post_matrix %>%
    rename(any_of(name_mapping))
}

# Apply renaming to both models
model1_post_named <- rename_posterior(model1_post, param_names)
model2_post_named <- rename_posterior(model2_post, param_names)


# Calculate summaries
model1_summary <- posterior_mat %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  group_by(parameter) %>%
  summarise(
    mean = mean(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    model = "Survey only"
  )

model2_summary <- isdm_posterior_mat %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  group_by(parameter) %>%
  summarise(
    mean = mean(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    model = "Survey & Citizen Science"
  )

combined_summary <- bind_rows(model1_summary, model2_summary)

# Create interval plot
ggplot(combined_summary, 
       aes(x = parameter, y = mean, ymin = lower, ymax = upper, 
           color = model)) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Posterior Intervals Comparison",
       x = "Parameter", y = "Coefficient Value") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)

