library(dplyr)
library(readr)

tasmin <- read_csv("csvs/tasmin_grid_new2.csv")
tasmax <- read_csv("csvs/tasmax_grid_new2.csv")

tasmean <- tasmax[, c("grid_id", "date", "mean_max_temp_7d_celsius")] %>%
  inner_join(tasmin[, c("grid_id", "date", "mean_min_temp_7d_celsius")], by = c("grid_id", "date")) %>%
  mutate(mean_temp_7d_celsius = (mean_max_temp_7d_celsius + mean_min_temp_7d_celsius) / 2)

# write.csv(tasmean,
#           file = "csvs/tasmean_grid_new.csv",
#           row.names = FALSE)


# PLOT

grid_avg_temp <- tasmean %>%
  group_by(grid_id) %>%
  summarize(avg_temperature = mean(mean_temp_7d_celsius, na.rm = TRUE)) %>%
  ungroup()

grid_clipped <- st_read("grids/grid_clipped_1km.gpkg")

grid_avg_sf <- grid_clipped %>%
  left_join(grid_avg_temp, by = "grid_id")

ggsave(filename = "maps/mean_temp_study_period.png",
       plot = ggplot(grid_avg_sf) +
         geom_sf(aes(fill = avg_temperature), color = NA) +
         scale_fill_viridis_c(option = "plasma", name = "Avg Temp\n(°C)") +
         labs(title = "Average Temperature by Grid Cell",
              subtitle = "Across all dates during study period") +
  theme_void(),
  device = "png",
  dpi =  300,
  units = "in",
  width = 7,
  height = 8
)


map.temp<-grid_avg_sf

map.temp<-st_make_valid(map.temp)

map.temp<-st_transform(map.temp, 27700)

res <- 1000  # Grid cell size: 10km

map.2 <- st_make_grid(map.temp, cellsize = res, what = "centers")  # Generate grid
map.2 <- st_sf(geometry = map.2)  # Convert to sf object
map.2<-map.2[st_transform(grid_avg_sf, 27700),]


map.2$avg_temperature <- st_join(map.2, map.temp)$avg_temperature



map.r<-rast(data.frame(x=st_coordinates(map.2)[,1],
                       y=st_coordinates(map.2)[,2],
                       z=map.2$avg_temperature),
            crs = st_crs(map.2)$wkt)

library(ggplot2)
library(cowplot)
library(mapview)
library(tidyterra)
avg_temp_map <- ggplot() +
  geom_spatraster(data = map.r, 
                  mask_projection = TRUE) +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid = element_blank()
  )   +
  scale_fill_viridis_c(
    option = "plasma",
    name = "",  
    na.value = "transparent"
  ) +
  theme_void()



ggsave(filename = "maps/new_mean_temp_study_period.png",
       plot = avg_temp_map,
       device = "png",
       dpi =  300,
       units = "in",
       width = 7,
       height = 8
)


########### ISLE OF ARRAN 


mask<-st_read("~/OneDrive - University of Glasgow/PhD/modelling/dataIntegration/MARS/arran_mask.shp")
plot(mask)
# mask<-st_transform(mask, st_crs(avg_rain_map))
# ARRAN<-avg_rain_map[mask,]

# Transform mask to match the raster CRS
mask <- st_transform(mask, st_crs(map.r))

# Now create the plot with the mask
avg_temp_map <- ggplot() +
  geom_spatraster(data = map.r) +
  geom_sf(data = mask, fill = NA, color = "black", linewidth = 1) +  # Add mask outline
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid = element_blank()
  ) +
  scale_fill_viridis_c(
    option = "plasma",
    name = "",  
    na.value = "transparent"
  ) +
  theme_void()
ggsave(filename = "maps/ARRAN_circled_mean_temp_study_period.png",
       plot = avg_temp_map,
       device = "png",
       dpi =  300,
       units = "in",
       width = 7,
       height = 8
)



# If you want to mask the raster data to only show within the polygon:
# Create a masked version of the raster
# Crop to the mask extent first, then mask
map.r_cropped <- crop(map.r, mask)
map.r_masked <- mask(map.r_cropped, mask)


# Then plot the masked version
avg_temp_map_masked <- ggplot() +
  geom_spatraster(data = map.r_masked) +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid = element_blank()
  ) +
  scale_fill_viridis_c(
    option = "plasma",
    name = "",  
    na.value = "transparent"
  ) +
  theme_void()

ggsave(filename = "maps/ARRAN_mean_temp_study_period.png",
       plot = avg_temp_map_masked,
       device = "png",
       dpi =  300,
       units = "in",
       width = 7,
       height = 8
)





# Simple histogram using terra for masked temperature data
hist(map.r_masked, 
     main = "Distribution of Temperature Values\nin Isle of Arran",
     xlab = "Temperature (°C)")

# Transform mask to match the CRS of map.2
mask <- st_transform(mask, st_crs(map.2))
map.2_masked <- map.2[mask, ]



# Create histogram from masked points
histogram_masked_arran <- ggplot(map.2_masked, aes(x = avg_temperature)) +
  geom_histogram(bins = 30, fill = "gray", colour = "white") +
  labs(title = "Distribution of Temperature Values\nin Isle of Arran",
       x = "Temperature (°C)", 
       y = "Frequency") +
  theme_minimal()

print(histogram_masked_arran)

ggsave(filename = "maps/ARRAN_mean_temp_histogram_masked_arran.png",
       plot = histogram_masked_arran,
       device = "png",
       dpi =  300,
       units = "in",
       width = 7,
       height = 8
)

# For comparison - all of Scotland
hist(map.r, 
     main = "Distribution of Temperature Values\nacross Scotland",
     xlab = "Temperature (°C)")

histogram_from_points <- ggplot(map.2, aes(x = avg_temperature)) +
  geom_histogram(bins = 30, fill = "gray", colour = "white") +
  labs(title = "Distribution of Temperature Values\nacross Scotland",
       x = "Temperature (°C)", 
       y = "Frequency") +
  theme_minimal()

print(histogram_from_points)

ggsave(filename = "maps/SCOTLAND_mean_temp_histogram.png",
       plot = histogram_from_points,
       device = "png",
       dpi =  300,
       units = "in",
       width = 7,
       height = 8
)



