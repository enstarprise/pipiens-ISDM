library(readr)
library(posterior)
library(cmdstanr)
library(bayesplot)
library(tidyverse)
library(sf)
library(dplyr)
library(ggplot2)
library(tmap)
library(terra)
library(tidyterra)


rm(list = ls())

plot_df <- read_csv("nngp_mapping_df_N10_M10.csv") 
load("01_data/processedCovariates/2.5km/stan_data_2.5km_NNGP_all_both.Rdata")
load("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/AIM 1/pipiens-ISDM/01_data/loaded_data_cov_2point5km.RData")

grid_sf <- merge(grid, plot_df, by = "grid_id", all.x = TRUE)
grid_sf <- st_set_crs(grid_sf, 4326)
grid_sf <- st_transform(grid_sf, 27700)


# Quick check
sum(!is.na(grid_sf$mean_abundance))
sum( is.na(grid_sf$mean_abundance))
  
# --------------------------------------
# RELATIVE ABUNDANCE

plot_df <- plot_df %>%
  mutate(
    # Proportion of total
    mean_prop   = mean_abundance / sum(mean_abundance, na.rm = TRUE),
    median_prop = median_abundance / sum(median_abundance, na.rm = TRUE),
    
    # Scaled 0-1
    mean_scaled   = (mean_abundance - min(mean_abundance, na.rm = TRUE)) /
      (max(mean_abundance, na.rm = TRUE) - min(mean_abundance, na.rm = TRUE)),
    median_scaled = (median_abundance - min(median_abundance, na.rm = TRUE)) /
      (max(median_abundance, na.rm = TRUE) - min(median_abundance, na.rm = TRUE)),
    
    # Z-score
    mean_z   = (mean_abundance - mean(mean_abundance, na.rm = TRUE)) /
      sd(mean_abundance, na.rm = TRUE),
    median_z = (median_abundance - mean(median_abundance, na.rm = TRUE)) /
      sd(median_abundance, na.rm = TRUE),
    
    # Percentile rank
    mean_percentile   = rank(mean_abundance,   na.last = "keep") / sum(!is.na(mean_abundance)),
    median_percentile = rank(median_abundance, na.last = "keep") / sum(!is.na(median_abundance))
  )


# --------------------------------------
# RASTER MAPPING

grid_sf <- st_make_valid(grid_sf)

# Rerun the relative abundance mutate on plot_df
plot_df <- plot_df %>%
  mutate(
    mean_scaled = (mean_abundance - min(mean_abundance, na.rm = TRUE)) /
      (max(mean_abundance, na.rm = TRUE) - min(mean_abundance, na.rm = TRUE))
  )

plot_df <- plot_df %>%
  mutate(
    mean_abundance_log   = log1p(mean_abundance),
    mean_abundance_sqrt  = sqrt(mean_abundance),
    mean_abundance_gamma = mean_abundance ^ 0.4,   # tweak 0.4 to taste
    mean_abundance_cbrt  = mean_abundance ^ (1/3)
  )

# Rebuild grid_sf with the correct CRS pipeline
grid_sf <- merge(grid, plot_df, by = "grid_id", all.x = TRUE)
grid_sf <- st_set_crs(grid_sf, 4326)   # assign true CRS first
grid_sf <- st_transform(grid_sf, 27700) # then reproject
grid_sf <- st_make_valid(grid_sf)

# Check
names(grid_sf)
sum(!is.na(grid_sf$mean_scaled))

# Visual comparison of transforms
plot_df %>%
  select(mean_abundance, mean_abundance_log, 
         mean_abundance_sqrt, mean_abundance_gamma) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 50) +
  facet_wrap(~name, scales = "free") +
  theme_minimal()


#### PLOT MEAN ABUNDANCE -------------------

map.2 <- st_make_grid(grid_sf, cellsize = 1000, what = "centers")
crs(map.2); crs(grid_sf)

map.2 <- st_sf(geometry = map.2)
map.2 <- map.2[grid_sf, ]  # clip to study area extent

map.2$mean_abundance <- st_join(map.2, grid_sf)$mean_abundance

map2.r <- rast(
  data.frame(
    x = st_coordinates(map.2)[, 1],
    y = st_coordinates(map.2)[, 2],
    z = map.2$mean_abundance
  ),
  crs = st_crs(grid_sf)$wkt
)

abundance_map <- ggplot() +
  geom_spatraster(data = map2.r) +
  scale_fill_viridis_c(
    option   = "inferno",
    name     = "Mean\nabundance",
    na.value = "transparent",
    trans    = "sqrt"       # swap for "sqrt" or remove entirely
  ) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background  = element_rect(fill = "transparent", colour = NA),
    legend.position  = "right"
  )

abundance_map

ggsave(filename = "04_outputs/maps/nngp_mean_intensity_1kmres_M10_N10.png", 
       plot = abundance_map,
       device = "png",
       width = 10, 
       height = 12,
       dpi = 300,
       units = "in")




#### RELATIVE INTENSITY ---------------
map.3 <- st_make_grid(grid_sf, cellsize = 1000, what = "centers")
crs(map.3); crs(grid_sf)

map.3 <- st_sf(geometry = map.3)
map.3 <- map.2[grid_sf, ]  # clip to study area extent

map.3$mean_scaled <- st_join(map.3, grid_sf)$mean_scaled

map3.r <- rast(
  data.frame(
    x = st_coordinates(map.2)[, 1],
    y = st_coordinates(map.2)[, 2],
    z = map.3$mean_scaled
  ),
  crs = st_crs(grid_sf)$wkt
)



rel_abundance_map <- ggplot() +
  geom_spatraster(data = map3.r) +
  scale_fill_viridis_c(
    option   = "magma",
    name     = "Relative\nintensity",
    na.value = "transparent"
  ) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background  = element_rect(fill = "transparent", colour = NA),
    legend.position  = "right"
  )

rel_abundance_map

ggsave(filename = "04_outputs/maps/nngp_rel_intensity_1kmres.png", 
       plot = abundance_map,
       device = "png",
       width = 10, 
       height = 12,
       dpi = 300,
       units = "in")


## EXTRAS: GGPLOT AND TMAP --------

# e.g. percentile rank
tm_shape(grid_sf) +
  tm_fill("mean_percentile",
          style   = "cont",
          palette = "magma",
          title   = "Percentile rank\n(mean abundance)") +
  tm_layout(main.title = "Relative Abundance (Percentile)",
            frame = FALSE)


# tmap v4 syntax
tm_shape(grid_sf) +
  tm_fill("mean_percentile",
          fill.scale = tm_scale_continuous(values = "magma"),
          fill.legend = tm_legend(title = "Percentile rank\n(mean abundance)")) +
  tm_title("Relative Abundance (Percentile)")



# ggplot
ggplot(grid_sf) +
  geom_sf(aes(fill = mean_scaled), colour = NA) +
  scale_fill_viridis_c(option = "magma", name = "Relative\nabundance") +
  labs(title = "Relative Abundance") +
  theme_minimal()

# tmap v4
tm_shape(grid_sf) +
  tm_fill("mean_scaled",
          fill.scale  = tm_scale_continuous(values = "magma"),
          fill.legend = tm_legend(title = "Relative abundance")) +
  tm_title("Relative Abundance")

# More polished version
tmap_mode("plot")
tmap_options(max.raster = c(plot = 1e6))  # Increase raster resolution

tm_shape(grid_sf) +
  tm_fill("mean_scaled",
          fill.scale = tm_scale_continuous(values = "viridis"),
          fill.legend = tm_legend(title = "Relative abundance",
                                  position = tm_pos_in("left", "bottom"),
                                  frame = FALSE)) +
  tm_layout(legend.frame = FALSE,
            legend.bg.color = "white",
            legend.bg.alpha = 0.8,
            legend.outside = FALSE,
            inner.margins = c(0.02, 0.02, 0.05, 0.02)) +
  tm_layout(legend.frame = FALSE,
            legend.bg.color = "white",
            legend.bg.alpha = 0.8,
            inner.margins = 0.05)


