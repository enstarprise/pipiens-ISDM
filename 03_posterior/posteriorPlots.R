
library(dplyr)
library(bayesplot)
library(ggplot2)
library(posterior)
library(viridis)
library(tidyverse)
library(Cairo)
library(readr)

posterior_mat <- readRDS("03_posterior/savedObjects/posterior_mat_nngp.rds")

posterior_mat_detection <- readRDS("03_posterior/savedObjects/posterior_mat_detection.rds")


# ============================================================================
# POSTEIROR VISUALISATION
# ============================================================================

main_params <- c(
  "beta0", "beta_temp", "beta_rain", "beta_land[1]", "beta_land[2]", "beta_land[3]",
  "beta_land[4]", "beta_land[5]", "beta_land[6]", "beta_land[7]", "beta_land[8]",
  "alpha0", "alpha_RH", "alpha_WS_rain", "sigma_trap", "alpha_trap[1]", "alpha_trap[1]",
  "alpha_trap[2]", "alpha_trap[3]", "alpha_trap[4]", "alpha_trap[5]",
  "delta0", "delta_poi", "delta_reports",
  "phi"
)


main_params <- c(
  "beta_temp", "beta_rain", "beta_land[1]", "beta_land[2]", "beta_land[3]",
  "beta_land[4]", "beta_land[5]", "beta_land[6]", "beta_land[7]", "beta_land[8]"
)

param_names <- c(
  "beta_temp" = "3 week prior\nmin temperature",
  "beta_rain" = "4 week prior\ncumulative rainfall", 
  "beta_land[1]" = "Wetland",
  "beta_land[2]" = "Woodland",
  "beta_land[3]" = "Freshwater",
  "beta_land[4]" = "Grassland",
  "beta_land[5]" = "Arable",
  "beta_land[6]" = "Livestock\nDensity",
  "beta_land[7]" = "Elevation",
  "beta_land[8]" = "Building\nDensity"
)

param_labels <- param_names[main_params]

intervals_plot_1 <- mcmc_intervals(
  posterior_mat,
  pars = main_params,
  prob = 0.8,
  prob_outer = 0.95
) + 
  scale_y_discrete(labels = param_labels) +  # correctly map labels
  ggtitle("Posterior Distributions with Credible Intervals") +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    axis.text.y = element_text(size = 12)
  )

intervals_plot_1

ggsave(filename = "03_posterior/priorPosteriorPlots/lambda_params_nngp.png",
       plot = intervals_plot_1,
       device = "png",
       width = 7,
       height = 8,
       units = "in",
       dpi = 300,
       background = "transparent")




# =======================================================================

lambda_params <- c(
  "beta_temp", "beta_rain", "beta_land[1]", "beta_land[2]", "beta_land[3]",
  "beta_land[4]", "beta_land[5]", "beta_land[6]", "beta_land[7]", "beta_land[8]"
)

param_names <- c(
  "beta_temp" = "3 week prior\nmin temperature",
  "beta_rain" = "4 week prior\ncumulative rainfall", 
  "beta_land[1]" = "Wetland",
  "beta_land[2]" = "Woodland",
  "beta_land[3]" = "Freshwater",
  "beta_land[4]" = "Grassland",
  "beta_land[5]" = "Arable",
  "beta_land[6]" = "Livestock\nDensity",
  "beta_land[7]" = "Elevation",
  "beta_land[8]" = "Building\nDensity"
)


# Calculate summary statistics
posterior_summary <- as.data.frame(posterior_mat) %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  group_by(parameter) %>%
  summarise(
    median = median(value),
    lower_80 = quantile(value, 0.10),
    upper_80 = quantile(value, 0.90),
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
  geom_errorbar(aes(xmin = lower_80, xmax = upper_80),
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

ggsave(filename = "~/pipiens-ISDM/stan/2.5km/key_params_maternGP144.png",
       plot = intervals_plot_2,
       device = "png",
       width = 5,
       height = 6,
       units = "in",
       dpi = 300)

#================================================

detection_params <- c(
  "alpha_RH", "alpha_WS_rain", "alpha_trap[1]", "alpha_trap[2]",
  "alpha_trap[3]", "alpha_trap[4]", "alpha_trap[5]",
  "delta_poi", "delta_reports"
)

param_names <- c(
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

posterior_mat_detection <- fit$draws(detection_params) %>% 
  posterior::as_draws_matrix()

saveRDS(posterior_mat_detection, "~/pipiens-ISDM/stan/2.5km/posterior_mat_detection.rds")

# Define colors and shapes for each parameter individually
param_styles <- data.frame(
  parameter = detection_params,
  color = c(rep("#2E86AB", 7), rep("#A23B72", 2)),  # First 7 survey params, last 2 CS params
  shape = c(rep(16, 7), rep(17, 2))  # 16 = circle, 17 = triangle
)

posterior_summary <- as.data.frame(posterior_mat) %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  group_by(parameter) %>%
  summarise(
    median = median(value),
    lower_80 = quantile(value, 0.10),
    upper_80 = quantile(value, 0.90),
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
  geom_errorbar(aes(xmin = lower_80, xmax = upper_80),
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
  scale_y_discrete(labels = param_names) +
  
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

ggsave(filename = "~/pipiens-ISDM/stan/2.5km/detection_params_maternGP144.png",
       plot = intervals_plot_3_thin,
       device = "png",
       width = 5,
       height = 6,
       units = "in",
       dpi = 300)


