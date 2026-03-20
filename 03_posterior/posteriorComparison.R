library(ggplot2)
library(tidybayes)
library(dplyr)
library(tidyr)

# Load data ----------------------------------------------------------------
# Model 1: Survey only
fit <- readRDS("02_model/stan/modelRuns/model_survey_nngp.rds")
main_params <- c(
  "beta0", "beta_temp", "beta_temp2", "beta_rain", "beta_rain2",
  "beta_land[1]", "beta_land[2]", "beta_land[3]",
  "beta_land[4]", "beta_land[5]", "beta_land[6]", "beta_land[7]", "beta_land[8]",
  "alpha0", "alpha_RH", "alpha_WS_rain", "sigma_trap", "alpha_trap[1]", "alpha_trap[1]",
  "alpha_trap[2]", "alpha_trap[3]", "alpha_trap[4]", "alpha_trap[5]",
  # "delta0", "delta_poi", "delta_reports",
  "phi"
)

lambda_params<- c(
  "beta_temp", "beta_temp2", "beta_rain", "beta_rain2",
  "beta_land[1]", "beta_land[2]", "beta_land[3]",
  "beta_land[4]", "beta_land[5]", "beta_land[6]", "beta_land[7]", "beta_land[8]"
)

posterior_mat <- fit$draws(lambda_params) %>% posterior::as_draws_matrix()
posterior_mat <- as.data.frame(posterior_mat)

# Model 2: Survey & Citizen Science
isdm_posterior_mat <- readRDS("03_posterior/savedObjects/posterior_mat_nngp_M10.rds")
isdm_posterior_mat <- as.data.frame(isdm_posterior_mat)


# Parameter naming ---------------------------------------------------------
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

# Rename columns in both dataframes
model1_named <- posterior_mat
model2_named <- isdm_posterior_mat

# Rename using dplyr::rename with all parameters
for(old_name in names(param_names)) {
  if(old_name %in% names(model1_named)) {
    names(model1_named)[names(model1_named) == old_name] <- param_names[old_name]
  }
  if(old_name %in% names(model2_named)) {
    names(model2_named)[names(model2_named) == old_name] <- param_names[old_name]
  }
}

# Check if renaming worked
print("Columns in model1 after renaming:")
print(names(model1_named))
print("Columns in model2 after renaming:")
print(names(model2_named))

# Calculate summaries ------------------------------------------------------
model1_summary <- model1_named %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  group_by(parameter) %>%
  summarise(
    mean = mean(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    model = "Survey only",
    .groups = "drop"
  )

model2_summary <- model2_named %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  group_by(parameter) %>%
  summarise(
    mean = mean(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    model = "Survey & Citizen Science",
    .groups = "drop"
  )

combined_summary <- bind_rows(model1_summary, model2_summary)

# Check what parameters we have
print("Unique parameters in combined_summary:")
print(unique(combined_summary$parameter))

# Define parameter order for display (using the new names)
parameter_order <- c(
  "Wetland", "Woodland", "Freshwater", "Grassland", "Arable",
  "Livestock\nDensity", "Elevation", "Buildings",
  "Temperature\n(3 week avg)", "Temperature²",
  "Rainfall\n(4 week\ncumulative)", "Rainfall²"
)

# Only keep parameters that exist in the data
existing_params <- intersect(parameter_order, unique(combined_summary$parameter))
parameter_order <- parameter_order[parameter_order %in% existing_params]

combined_summary <- combined_summary %>%
  filter(parameter %in% existing_params) %>%
  mutate(parameter = factor(parameter, levels = rev(parameter_order)))  # rev for coord_flip

# Plot 1: Interval plot ----------------------------------------------------
ggplot(combined_summary, 
       aes(x = parameter, y = mean, ymin = lower, ymax = upper, 
           color = model, shape = model)) +
  geom_pointrange(position = position_dodge(width = 0.6), size = 0.8) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_line(color = "gray90", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  labs(
    title = "Posterior Coefficient Comparison",
    subtitle = "Survey-only model vs. Integrated model with citizen science data",
    x = "", y = "Coefficient Value"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
  scale_color_manual(values = c("Survey only" = "#2E86AB", 
                                "Survey & Citizen Science" = "#A23B72")) +
  scale_shape_manual(values = c("Survey only" = 16, 
                                "Survey & Citizen Science" = 17))

# Plot 2: Density plots ----------------------------------------------------
# Prepare long format data for densities
model1_long <- model1_named %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  mutate(model = "Survey only") %>%
  filter(parameter %in% existing_params)

model2_long <- model2_named %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  mutate(model = "Survey & Citizen Science") %>%
  filter(parameter %in% existing_params)

combined_long <- bind_rows(model1_long, model2_long) %>%
  mutate(parameter = factor(parameter, levels = parameter_order))

# Density plot
ggplot(combined_long, aes(x = value, fill = model)) +
  geom_density(alpha = 0.5, adjust = 1.2) +
  facet_wrap(~parameter, scales = "free", ncol = 3) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "gray95", color = NA),
    panel.spacing = unit(1, "lines")
  ) +
  labs(
    title = "Posterior Density Comparison",
    subtitle = "Survey-only model vs. Integrated model with citizen science data",
    x = "Coefficient Value", y = "Density"
  ) +
  scale_fill_manual(values = c("Survey only" = "#2E86AB", 
                               "Survey & Citizen Science" = "#A23B72")) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)

# Calculate and print comparison metrics -----------------------------------
comparison_metrics <- map_dfr(existing_params, function(p) {
  p1 <- model1_named[[p]]
  p2 <- model2_named[[p]]
  diff_samples <- p1 - p2
  
  # Remove line breaks for table display
  p_clean <- gsub("\n", " ", p)
  
  data.frame(
    Parameter = p_clean,
    `Survey only mean` = sprintf("%.3f", mean(p1)),
    `Survey only CI` = sprintf("[%.3f, %.3f]", quantile(p1, 0.025), quantile(p1, 0.975)),
    `Integrated mean` = sprintf("%.3f", mean(p2)),
    `Integrated CI` = sprintf("[%.3f, %.3f]", quantile(p2, 0.025), quantile(p2, 0.975)),
    `Difference` = sprintf("%.3f", mean(diff_samples)),
    `P(Integrated > Survey)` = sprintf("%.3f", mean(diff_samples < 0)),  # Note: <0 because Survey - Integrated
    `Uncertainty ratio` = sprintf("%.2f", 
                                  (quantile(p1, 0.975) - quantile(p1, 0.025)) /
                                    (quantile(p2, 0.975) - quantile(p2, 0.025))),
    check.names = FALSE
  )
})

print("Comparison Metrics:")
print(comparison_metrics)

# Optional: Save the plot
ggsave("posterior_comparison.png", width = 10, height = 8, dpi = 300)