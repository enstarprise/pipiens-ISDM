library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)
library(readr)
library(sf)
library(tibble)


load("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/AIM 1/pipiens-ISDM/01_data/loaded_data_cov_2point5km.RData")

draws_resid <- fit$draws("residual_y")

resid_df <- as_draws_df(draws_resid) %>%
  mutate(sample_id = row_number()) %>%
  pivot_longer(
    cols = starts_with("residual_y"),
    names_to = "obs_id",
    names_pattern = "residual_y\\[(\\d+)\\]",
    values_to = "residual"
  ) %>%
  mutate(obs_id = as.integer(obs_id))

survey_index_df <- tibble(
  obs_id = seq_along(stan_data$site_to_grid),
  grid_id = stan_data$site_to_grid,  # Stan grid ID
  date_id = stan_data$date_y
)

resid_df <- resid_df %>%
  left_join(survey_index_df, by = "obs_id")


grid_id_lookup_df <- tibble(
  original_grid_id = as.integer(names(grid_id_lookup)),
  grid_id_stan     = as.integer(grid_id_lookup)
)

resid_df <- resid_df %>%
  left_join(grid_id_lookup_df, by = c("grid_id" = "grid_id_stan"))


centroids_obs_df <- centroids_2point5km %>%
  filter(grid_id %in% unique(resid_df$original_grid_id)) %>%
  rename(geometry = geom)

resid_df <- resid_df %>%
  left_join(centroids_obs_df, by = c("original_grid_id" = "grid_id")) %>%
  mutate(
    x = st_coordinates(geometry)[,1],
    y = st_coordinates(geometry)[,2]
  )

resid_summary_df <- resid_df %>%
  select(sample_id, original_grid_id, grid_id, date_id, residual, x, y)

write_csv(resid_summary_df, "survey_residual_summary_with_dates.csv")


grid_coords_survey <- resid_summary_df %>%
  select(original_grid_id, x, y) %>%
  distinct(original_grid_id, .keep_all = TRUE)

write_csv(grid_coords_survey, "grid_coords_survey.csv")



# Aggregate residuals across dates for each sample and grid
resid_grid_mean_df <- resid_df %>%
  group_by(sample_id, grid_id) %>%
  summarise(
    residual_mean = mean(residual),  # mean over all dates
    .groups = "drop"
  )

# Join with the original grid IDs
resid_grid_mean_df <- resid_grid_mean_df %>%
  left_join(grid_id_lookup_df, by = c("grid_id" = "grid_id_stan"))

# Add coordinates
resid_grid_mean_df <- resid_grid_mean_df %>%
  left_join(
    centroids_2point5km %>% select(original_grid_id, geometry),
    by = "original_grid_id"
  ) %>%
  mutate(
    x = st_coordinates(geometry)[,1],
    y = st_coordinates(geometry)[,2]
  ) %>%
  select(sample_id, original_grid_id, residual_mean, x, y)

# Save
write_csv(resid_grid_mean_df, "residual_summary_grid_mean.csv")


library(dplyr)
library(sf)
library(dplyr)
library(sf)

# 1. select only sample 1
resid_grid_sample1_df <- resid_grid_mean_df %>%
  filter(sample_id == 1)

# Save coordinates and residuals
write_csv(resid_grid_sample1_df, "residual_grid_sample1.csv")


# =================================================================

fit <- readRDS("residuals_fit.rds")
draws <- fit$draws("lambda_base")


lambda_df <- as_draws_df(draws) %>%
  mutate(sample_id = row_number()) %>%
  pivot_longer(
    cols = starts_with("lambda_base"),
    names_to = c("grid_id", "date_id"),
    names_pattern = "lambda_base\\[(\\d+),(\\d+)\\]",
    values_to = "lambda_base"
  ) %>%
  mutate(
    grid_id = as.integer(grid_id),
    date_id = as.integer(date_id)
  ) %>%
  select(sample_id, grid_id, date_id, lambda_base)


# get and subset grid coords
load("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/AIM 1/pipiens-ISDM/01_data/loaded_data_cov_2point5km.RData")
# grid <- st_read("01_data/grids/grid_clipped_2point5km.gpkg")

grid_id_lookup_df <- tibble(
  grid_id = seq_along(grid_id_lookup),   
  original_grid_id = grid_id_lookup
)

lambda_df <- lambda_df %>%
  left_join(grid_id_lookup_df, by = "grid_id")

coords_obs_df <- as.data.frame(coords_obs)

stopifnot(
  nrow(coords_obs_df) == length(observed_grid_ids),
  nrow(coords_obs_df) == max(lambda_df$grid_id)
)

write_csv(coords_obs_df, "grid_coords_observed.csv")

write_csv(
  lambda_df %>% select(sample_id, grid_id, date_id, lambda_base),
  "lambda_base_posterior.csv"
)

# ============================================================================
# FROM MARS MVLS 

lambda_df <- read_csv("lambda_base_posterior_base_model.csv") # from MARS MVLS 

load("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/AIM 1/pipiens-ISDM/01_data/loaded_data_cov_2point5km.RData")


#  grid_id_lookup named vector to a dataframe
lookup_table<- data.frame(
  grid_id = as.integer(names(grid_id_lookup)),  # The actual grid IDs
  grid_index = as.integer(grid_id_lookup)       # The indices (1, 2, 3, ...)
)


merged_data <- lambda_df %>%
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


lambda_df <- final_df

write_csv(grid_coords, "grid_coords.csv")

write_csv(
  lambda_df %>% select(sample_id, grid_id, date_id, lambda_base),
  "lambda_base_posterior.csv"
)

# =============================================================================

# Extract as data frame (long format)
lambda_resid_df <- fit$draws("lambda_residuals", format = "df")

# iteration, chain, draw, lambda_residuals[1], lambda_residuals[2]
lambda_resid_long <- lambda_resid_df %>%
  select(-.chain, -.iteration, -.draw) %>%
  mutate(draw_id = row_number()) %>%
  pivot_longer(
    cols = starts_with("lambda_residuals"),
    names_to = "obs_id",
    values_to = "lambda_residual"
  ) %>%
  mutate(obs_id = as.integer(str_extract(obs_id, "\\d+")))

write.csv(lambda_resid_long, "lambda_residuals_full.csv", row.names = FALSE)


# filter to just one draw (e.g., draw 1)
lambda_resid_one_draw <- lambda_resid_df %>%
  filter(.draw == 1) %>%
  select(-.chain, -.iteration, -.draw) %>%
  pivot_longer(
    cols = starts_with("lambda_residuals"),
    names_to = "obs_id",
    values_to = "lambda_residual"
  ) %>%
  mutate(obs_id = as.integer(str_extract(obs_id, "\\d+")))

write.csv(lambda_resid_one_draw, "lambda_residuals_one_draw.csv", row.names = FALSE)


# summarize across all draws
lambda_resid_summary <- lambda_resid_df %>%
  select(-.chain, -.iteration, -.draw) %>%
  pivot_longer(
    cols = starts_with("lambda_residuals"),
    names_to = "obs_id",
    values_to = "lambda_residual"
  ) %>%
  mutate(obs_id = as.integer(str_extract(obs_id, "\\d+"))) %>%
  group_by(obs_id) %>%
  summarise(
    lambda_residual_mean = mean(lambda_residual),
    lambda_residual_median = median(lambda_residual),
    lambda_residual_sd = sd(lambda_residual),
    lambda_residual_q025 = quantile(lambda_residual, 0.025),
    lambda_residual_q975 = quantile(lambda_residual, 0.975)
  )

write.csv(lambda_resid_summary, "lambda_residuals_summary.csv", row.names = FALSE)