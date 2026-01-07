# adult to larvae assocaitions -- raw data

library(lubridate)
library(dplyr)
library(ggplot2)
library(tidyr)

larvae_df <- larvae_df %>%
  mutate(week_start = floor_date(date, "week", week_start = 1)) # 1 = monday
survey_df <- survey_df %>%
  mutate(week_start = floor_date(date, "week", week_start = 1))

process_df <- function(df, name) {
  df %>%
    mutate(
      week_start = floor_date(date, "week"), # Sunday as week start
      source = name
    ) %>%
    group_by(week_start, source, Site_name) %>%
    summarize(count = sum(Culex), .groups = "drop") # Sum counts by week
}

# process all dataframes
combined <- bind_rows(
  process_df(survey_df, "adults"),
  process_df(larvae_df, "larvae")
)

# --- count data ------
hist(larvae_df$Culex)
larvae_prop_zeros <- mean(larvae_df$Culex == 0, na.rm = TRUE) 

ggsave(filename = "descriptivePlots/larvae_count.png", 
       plot = ggplot(data = larvae_df, aes(x = Culex)) +
         geom_bar() +
         labs(title = "Culex pipiens larvae sampled", y = "Frequency")+
         theme_minimal(),
       device = "png",
       width = 7,
       height = 8,
       units = "in",
       dpi = 300
)

ggsave(filename = "descriptivePlots/adult_count.png", 
       plot = ggplot(data = survey_df, aes(x = Culex)) +
         geom_bar() +
         labs(title = "Culex pipiens adults trapped", y = "Frequency")+
         theme_minimal(),
       device = "png",
       width = 7,
       height = 8,
       units = "in",
       dpi = 300
)

# --- larvae covariates ----
#> interested in covariates associated with breeding sites

ggsave(filename = "descriptivePlots/larvae_habitat_category.png", 
       plot = ggplot(data = larvae_df, aes(x = Habitat_category)) +
  geom_bar(fill = "grey80", col = "grey60") +
  labs(title = "Habitat Category", y = "Frequency", x = "Habitat category") +  # 10 are artificial
  theme_minimal(),
  device = "png",
  width = 7,
  height = 8,
  units = "in",
  dpi = 300
  )

ggsave(filename = "descriptivePlots/larvae_habitat_type.png", 
       plot = ggplot(data = larvae_df, aes(x = Habitat_type)) +
  geom_bar(fill = "grey80", col = "grey60") +
  labs(title = "Habitat type", y = "Frequency", x = "Habitat type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)),
device = "png",
width = 7,
height = 8,
units = "in",
dpi = 300
)

ggsave(filename = "descriptivePlots/larvae_habitat_size.png",
       plot = ggplot(data = larvae_df, aes(x = Habitat_size)) +
  geom_bar(fill = "grey80", col = "grey60") +
  labs(title = "Habitat size", y = "Frequency", x = "Habitat size") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)),
  device = "png",
  width = 7,
  height = 8,
  units = "in",
  dpi = 300
)

ggsave(filename = "descriptivePlots/larvae_veg_cover.png",
       plot = ggplot(data = larvae_df, aes(x = Vegetation_cover_scaled)) +
  geom_histogram(fill = "grey80", col = "grey60") +
  labs(title = "Vegetation cover (scaled)", y = "Frequency", x = "Vegetation cover (scaled)") +
  theme_minimal(),
  device = "png",
  width = 7,
  height = 8,
  units = "in",
  dpi = 300)

# --- adult covariates ---
#> interested in environmental and climate factors leading to the abundance
#> 
#> land use variabels; 5km grid cells 
#> survey_df columns 15 - 24
#> 
library(ggplot2)
library(tidyr)

## LANDCOVER
scotland <- st_read("gridExtractions/scotlandBoundary/scotland_boundary.shp") %>%
  vect()
grid_clipped <- st_read("grids/grid_clipped_1km.gpkg")
grid_clipped <- st_transform(grid_clipped, st_crs(scotland))
st_crs(grid_clipped)

grid_terra <- vect(grid_clipped)
grid_terra <- crop(grid_terra,  ext(scotland))

landcover_1km <- read_csv("csvs/landcover_1km.csv")


proportions_with_original <- landcover_1km %>%
  mutate(total_area = rowSums(dplyr::select(., -grid_id))) %>%
  mutate(across(-c(grid_id, total_area), ~ . / total_area, .names = "prop_{.col}")) %>%
  dplyr::select(-total_area)


landcover_long <- pivot_longer(proportions_with_original,
                        cols = 24:45,
                        names_to = "variable",
                        values_to = "value")

ggplot(landcover_long) +
  geom_boxplot(aes(x = variable, y = value)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate labels 90 degrees


# create a grouping variable
landcover_long_grouped <- landcover_long %>%
  mutate(
    landcover_group = case_when(
      variable %in% c("prop_Coniferous_woodland_km2", "prop_Deciduous_woodland_km2") ~ "Woodland",
      variable %in% c("prop_Acid_grassland_km2", "prop_Neutral_grassland_km2", "prop_Heather_grassland_km2", "prop_Improved_grassland_km2", "prop_Heather_km2") ~ "Grassland",
      variable %in% c("prop_Saltmarsh_km2", "prop_Bog_km2", "prop_Fen_km2") ~ "Wetland",
      variable %in% c("prop_Supralittoral_rock_km2", "prop_Supralittoral_sediment_km2", "prop_Littoral_rock_km2", "prop_Littoral_sediment_km2") ~ "Coastal",
      variable == "prop_Arable_km2" ~ "Arable",
      variable == "prop_Freshwater_km2" ~ "Freshwater",
      variable == "prop_Saltwater_km2" ~ "Saltwater",
      
      variable == "prop_Unnamed_km2" ~ "Unnamed",
      TRUE ~ "Other"
    ),
    
    # Keep original names for within-group differentiation
    landcover_detail = case_when(
      variable == "prop_0" ~ "Unnamed",
      variable == "prop_1" ~ "Deciduous Woodland",
      variable == "prop_2" ~ "Coniferous Woodland", 
      variable == "prop_3" ~ "Arable",
      variable == "prop_4" ~ "Improved Grassland",
      variable == "prop_5" ~ "Neutral Grassland",
      variable == "prop_6" ~ "Calcareous Grassland",
      variable == "prop_7" ~ "Acid Grassland",
      variable == "prop_8" ~ "Fen",
      variable == "prop_9" ~ "Heather",
      variable == "prop_10" ~ "Heather Grassland",
      variable == "prop_11" ~ "Bog",
      variable == "prop_12" ~ "Inland Rock",
      variable == "prop_13" ~ "Saltwater",
      variable == "prop_14" ~ "Freshwater",
      variable == "prop_15" ~ "Supralittoral Rock",
      variable == "prop_16" ~ "Supralittoral Sediment",
      variable == "prop_17" ~ "Littoral Rock",
      variable == "prop_18" ~ "Littoral Sediment",
      variable == "prop_19" ~ "Saltmarsh",
      variable == "prop_20" ~ "Urban",
      variable == "prop_21" ~ "Suburban"
    )
  )

# Order the groups logically
group_order <- c("Woodland", "Grassland", "Wetland", "Arable", 
                 "Freshwater", "Coastal", "Rock", "Urban", "Suburban", "Unnamed")

landcover_long_grouped$landcover_group <- factor(
  landcover_long_grouped$landcover_group,
  levels = group_order
)


ggplot(landcover_long_grouped) +
  geom_boxplot(aes(x = landcover_group, y = value, fill = landcover_group)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "Landcover Group", y = "Proportion", title = "Landcover Proportions by Group") +
  scale_fill_brewer(palette = "Set3")


library(ggplot2)
library(dplyr)

# Simply rename the column for clarity and plot it
landcover_long_grouped %>%
  rename(proportion = value) %>% # Change 'value' to whatever your column is
  ggplot(aes(x = landcover_group, y = proportion)) +
  geom_boxplot() + # This will now show the distribution of proportions per cell
  labs(
    title = "Landcover Composition per 1km² Grid Cell",
    x = "Landcover Type",
    y = "Proportion of Cell Area"
  ) +
  theme_minimal()



# landcover_map <- grid_clipped %>%
#   left_join(proportions_with_original, by = "grid_id")
# 
# 
# ggplot(landcover_map) +
#   geom_sf(aes(fill = prop_Urban_km2)) +
#   scale_fill_viridis_c() +
#   labs(title = "Urban land cover proportion") +
#   theme_minimal()


library(ggplot2)
library(tidyr)
# join the landcover proportion to the survey df

# landcover_long <- pivot_longer(survey_df, 
#                                cols = 15:24, 
#                                names_to = "variable", 
#                                values_to = "value")


landcover_survey <- left_join(landcover_long_grouped[,c("grid_id", "variable","value", "landcover_group")],
                              survey_df[,c("grid_id")], by = "grid_id")



ggplot(landcover_survey, aes(x = landcover_group, y = value)) +
  geom_boxplot() + 
  labs(x = "UKCEH Landcover type", 
       y = "Sum of area per 1km grid cell (km2)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(filename = "descriptivePlots/landcover_area_plot.png",
       plot = ggplot(landcover_long, aes(x = landcover_group, y = value)) +
         geom_boxplot() + 
         labs(x = "UKCEH Landcover type", 
              y = "Sum of area per 5km grid cell (km2)") +
         theme_minimal() +
         theme(axis.text.x = element_text(angle = 45, hjust = 1)),
       device = "png",
       width = 7,
       height = 8,
       units = "in",
       dpi = 300
)


##TEMPERATURE
hist(full_grid_dates$tasmean_avg)
ggsave(filename = "descriptivePlots/tasmean.png",
       plot = ggplot(full_grid_dates, aes(x = tasmean_avg)) +
         geom_bar() + 
         labs(x = "1-week prior average temperature", 
              y = "Frequency") +
         theme_minimal() +
         theme(axis.text.x = element_text(angle = 0, hjust = 1)),
       device = "png",
       width = 7,
       height = 8,
       units = "in",
       dpi = 300
)

ggsave(filename = "descriptivePlots/rainfall_28d.png",
       plot = ggplot(full_grid_dates, aes(x = rainfall_28d)) +
         geom_bar() + 
         labs(x = "4-week cumulative precipitation", 
              y = "Frequency") +
         theme_minimal() +
         theme(axis.text.x = element_text(angle = 0, hjust = 1)),
       device = "png",
       width = 7,
       height = 8,
       units = "in",
       dpi = 300
)

ggplot(full_grid_dates, aes(as.numeric(rainfall_28d))) +
  geom_bar() + 
  labs(x = "4-week cumulative precipitation", 
       y = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

hist(full_grid_dates$rainfall_28d)


hist(landcover_matrices$z_woodland)


par(mfrow = c(2, 2)) 
for(i in seq_along(landcover_matrices)) {
  matrix_name <- if(!is.null(names(landcover_matrices))) {
    names(landcover_matrices)[i]
  } else {
    paste("Matrix", i)
  }
  hist(landcover_matrices[[i]], 
       main = paste("Histogram of", matrix_name),
       xlab = "Value",
       col = "grey90",
       border = "black")
}

hist(z_livestock)
hist(z_elevation)



### other landcover plots
par(mfrow = c(3, 2)) 

columns_to_plot <- c(3:10, 21, 23, 25)
# for  each column .... create histogram
for(col in columns_to_plot) {
  hist(lambda_grid[[col]], 
       main = paste("Histogram of", names(lambda_grid)[col]),
       xlab = names(lambda_grid)[col],
       col = "lightblue")
}
par(mfrow = c(1, 1))



library(ggplot2)
# prior samples
set.seed(123)
n_samples <- 10000
beta_samples <- rnorm(n_samples, 0, 1.5)
odds_ratio_samples <- exp(beta_samples)

plot_data <- data.frame(
  beta = beta_samples,
  odds_ratio = odds_ratio_samples
)

# plot 1: on coefficient scale (log-odds)
p1 <- ggplot(plot_data, aes(x = beta)) +
  geom_density(fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = c(-2.94, 2.94), linetype = "dashed", color = "red") +
  labs(title = "Prior: Coefficient (Log-Odds Scale)",
       subtitle = "Normal(0, 1.5)",
       x = "β coefficient (log-odds)",
       y = "Density") +
  theme_minimal()

# plot 2: on odds ratio scale (multiplicative)
p2 <- ggplot(plot_data, aes(x = odds_ratio)) +
  geom_density(fill = "forestgreen", alpha = 0.7) +
  geom_vline(xintercept = c(0.053, 18.9), linetype = "dashed", color = "red") +
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                labels = c("0.01x", "0.1x", "1x", "10x", "100x")) +
  labs(title = "Prior: Effect on Abundance (Fold Change)",
       subtitle = "exp(β) ~ Log-Normal(0, 1.5)",
       x = "Fold Change in Abundance",
       y = "Density") +
  theme_minimal()

library(patchwork)
p1 / p2 # combine plots



# What your prior expects vs what literature shows
library(ggplot2)

# thee prior expectation
prior_traps <- rnorm(10000, 0, 0.8)  # trap effects around mean
prior_probs <- plogis(prior_traps)    # must ! convert to probabilities

# literature values
literature <- data.frame(
  trap = c("Lower end prob.", "Higher end prob."),
  logit = c(-0.33, 1.52),
  prob = c(0.42, 0.82)
)

ggplot() +
  geom_density(aes(x = prior_probs), fill = "blue", alpha = 0.5) +
  geom_vline(data = literature, aes(xintercept = prob, color = trap), 
             size = 2, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 1)) +
  labs(title = "Prior Expectation vs Literature Reality",
       x = "Detection Probability",
       y = "Density") +
  theme_minimal()


# --- larvae adult association ------

ggsave(filename = "descriptivePlots/adult-larvae-counts.png",
  plot = ggplot(combined, aes(x = week_start, y = count, color = source)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Weekly Count Comparison",
    x = "Week Start Date",
    y = "Count",
    color = "Data Source"
  ) +
  facet_wrap(~Site_name, scale = "free") +
  theme_minimal() +
  scale_x_date(date_labels = "%b %d", date_breaks = "4 week") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)),
  device = "png",
  width = 12,
  height = 10,
  units = "in",
  dpi = 300
)

### 2023 only
combined %>%
  filter(year(week_start) == 2023) %>%
  ggplot(aes(x = week_start, y = count, color = source)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Weekly Count Comparison (2023)",
    x = "Week Start Date",
    y = "Count",
    color = "Data Source"
  ) +
  scale_x_date(
    date_labels = "%b %d",
    date_breaks = "4 week",
    limits = c(as.Date("2023-01-01"), as.Date("2023-12-31"))
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.margin = margin(10, 10, 30, 10)
  )


# --- counts per week/site: spatial autocorrelation check ---

ggplot(combined, aes(x = Site_name, y = count, color = source)) +
  geom_line() +
  geom_point()  +
  facet_wrap(~week_start, scale = "free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


combined %>%
  filter(year(week_start) == 2023) %>%
  group_by(Site_name, source, week_start) %>%
  summarise(total = n()) %>%
  ggplot(aes(x = Site_name, y = total, color = source)) +
  geom_line() +
  geom_point() +
  facet_wrap(~week_start, scale = "free") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.margin = margin(10, 10, 30, 10)
  )


# --- prior distributions -----
hist(rnorm(1000, 0 , 0.5))
