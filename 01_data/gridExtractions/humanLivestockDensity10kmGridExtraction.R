library(raster)
library(terra) 
library(sf)
library(mapview)
library(readr)

rm(list = ls())
gc()

scotland <- st_read("gridExtractions/scotlandBoundary/scotland_boundary.shp") %>%
  vect()
plot(scotland)
st_crs(scotland)

grid_clipped <- st_read("grids/grid_clipped_1km.gpkg")
grid_clipped <- st_transform(grid_clipped, st_crs(scotland))
st_crs(grid_clipped)

grid_terra <- vect(grid_clipped) # convrt to spatvector
st_crs(grid_terra)

### --- human density ----- ###
#> 2020 human population density across the UK
#> also performed in QGIS for the 5km extraction 
#> R extraction:
################################
humanDensity <- rast("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/GHS_POP_E2020_GLOBE_R2023A_54009_100_V1_0_R3_C18/GHS_POP_E2020_GLOBE_R2023A_54009_100_V1_0_R3_C18.tif")
humanDensity <- project(humanDensity, scotland)
st_crs(humanDensity)

# crop and mask to scotland's shapefile
humanDensity <- crop(humanDensity,  ext(scotland))
human_masked <- mask(humanDensity, scotland)
humanDensity <- human_masked
plot(humanDensity)

# verify alignment since that are both at 1km
grid_sf <- st_as_sf(grid_terra)
# mapview(humanDensity, alpha.regions = 0.5, layer.name = "Population Density") + 
  mapview(humanDensity, maxpixels = ncell(humanDensity))   + 
  mapview(grid_sf, col.regions = "red", alpha.regions = 0.2, layer.name = "Grid Cells")

### --- EXTRACTION --- ###
grid_1km_human <- grid_clipped
grid_1km_human$human_density <- terra::extract(
  humanDensity,
  grid_terra,
  fun = mean, # mean for density, sum for population counts
  na.rm = TRUE
)[, 2]

summary(grid_1km_human$human_density)
hist(grid_1km_human$human_density, breaks = 50)

library(tmap)
tmap_mode("view")  # Interactive mode (like mapview)
tm_shape(grid_1km_human) +
  tm_fill(
    col = "human_density",
    palette = "YlOrRd",  # Color scale
    style = "quantile",  # Classification method ("equal", "pretty", "quantile")
    n = 5,              # Number of classes
    alpha = 0.7,        # Transparency
    title = "Population Density (per km²)"
  ) +
  tm_borders(col = "gray30", lwd = 0.5)  # Grid cell borders

tmap_mode("plot")  # Static mode (for saving as PNG/PDF)
tm_shape(grid_1km_human) +
  tm_fill(
    col = "human_density",
    palette = "viridis", 
    style = "quantile",  # "log10_pretty", "cont", "quantile",
    title = "Population Density (per km²)"
  ) +
  tm_borders() +
  tm_layout(
    legend.position = c("right", "top"),
    frame = FALSE
  )


write_csv(grid_1km_human, "csvs/humanDensity_1km.csv")

### --- livestock density--- ###
#> gridded livestock density (GLW4) global-2020 10km 
#> spatial resolution: 
#> https://data.amerigeoss.org/dataset/9d1e149b-d63f-4213-978b-317a8eb42d02
#> https://www.nature.com/articles/sdata2018227
#>  buffalo, cattle, sheep, goats, pigs and chicken. 
#>  individual species datasets are available at global extent and 
#>  5 minutes of arc resolution (approx. 10 km at the equator)
#>  EPSG:4326 - WGS84 - Geographic Coordinate System (lat/long)
###################################


livestockDensity <- rast("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/livestockDensity/GLW4-2020.D-DA.GLEAM3-ALL-LU.tif")
crs(livestockDensity) <- "EPSG:4326" # or whatever its original CRS is
st_crs(livestockDensity)
plot(livestockDensity)

# crop and mask to scotland's shpaefile
livestockDensity <- crop(livestockDensity,  ext(scotland))
livestock_masked <- mask(livestockDensity, scotland)
livestockDensity <- raster(livestock_masked)

mapview(livestockDensity, 
        col.regions = viridisLite::viridis,
        layer.name = "Livestock Density")

plot(livestockDensity)
plot(scotland, add = TRUE, alpha = 0.7)
# plot(grid_clipped$geom, add = TRUE, alpha = 0.7)

#---- disaggregate ----
# terra format for extraction
livestock_terra <- rast(livestockDensity)  # convert to spatraster

# resample to 5km : disaggregate
livestock_1km <- disagg(livestock_terra, fact = 10, method = "bilinear")

# replace NA with 0 just temporarily for focal to work
livestock_1km_tmp <- classify(livestock_1km, rcl = matrix(c(NA, 0), ncol = 2), right = FALSE)

# define focal window (3x3 matrix of 1s)
window <- matrix(1, nrow = 3, ncol = 3)

# apply focal with custom function
livestock_1km_filled <- focal(
  livestock_1km_tmp,
  w = window,
  fun = function(x) {
    center <- x[5]  # Middle cell
    neighbors <- x[-5]  # All neighbors
    # if center is 0, calculate mean of non-zero neighbors
    if (center == 0) {
      nz <- neighbors[neighbors != 0]
      if (length(nz) > 0) mean(nz, na.rm = TRUE) else NA
    } else {
      center
    }
  },
  na.policy = "omit",
  filename = "",  # avoid writing to disk
  overwrite = TRUE
)

plot(livestock_1km, main = "Original")
plot(livestock_1km_filled, main = "NA-filled using focal")


# ----- extract mean values to grid -----
grid_1km_livestock <- grid_clipped

grid_1km_livestock$livestock_density <- terra::extract(
  livestock_1km, #livestock_1km
  grid_terra, 
  fun = mean,
  na.rm = TRUE
)[, 2]

grid_1km_livestock$livestock_density_filled <- terra::extract(
  livestock_1km_filled,
  grid_terra, 
  fun = mean,
  na.rm = TRUE
)[, 2]

head(grid_1km_livestock$livestock_density)

which_na <- which(is.na(grid_1km_livestock$livestock_density_filled))
plot(grid_terra[which_na], border = "red")

# mean_value <- mean(values(livestock_5km), na.rm = TRUE)
# grid_5km_livestock$livestock_density[is.na(grid_5km_livestock$livestock_density)] <- mean_value


# --- visualise ---
grid_1km_livestock_plot <- grid_1km_livestock %>%
  mutate(livestock_density = ifelse(is.na(livestock_density), 0, livestock_density))

ggsave(filename = "maps/livestockDensity.png",
       plot = ggplot(grid_1km_livestock_plot) +
         geom_sf(aes(fill = livestock_density), color = NA) +
         scale_fill_viridis_c(option = "plasma", 
                              name = "Livestock\nDensity") +
         labs(title = "Livestock Density",
              subtitle = "per 1km Grid Cell") +
         theme_void() +  # Removes axes, grid lines, and background
         theme(
           axis.text = element_blank(),        # Remove axis text
           axis.title = element_blank(),       # Remove axis titles
           axis.ticks = element_blank(),       # Remove axis ticks
           panel.grid = element_blank(),       # Remove grid lines
           plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center title
           plot.subtitle = element_text(hjust = 0.5, size = 12),  # Center subtitle
           legend.position = "right"           # Position legend
         ),
       device = "png",
       dpi = 300,
       width = 7,
       height = 9,
       units = "in",
       background = "transparent"
)

map.temp<-grid_1km_livestock_plot

map.temp<-st_make_valid(map.temp)

map.temp<-st_transform(map.temp, 27700)

res <- 1000  # Grid cell size: 10km

map.2 <- st_make_grid(map.temp, cellsize = res, what = "centers")  # Generate grid
map.2 <- st_sf(geometry = map.2)  # Convert to sf object
map.2<-map.2[st_transform(grid_1km_livestock_plot, 27700),]


map.2$livestock_density <- st_join(map.2, map.temp)$livestock_density

map.r<-rast(data.frame(x=st_coordinates(map.2)[,1],
                       y=st_coordinates(map.2)[,2],
                       z=map.2$livestock_density),
            crs = st_crs(map.2)$wkt)

library(ggplot2)
library(cowplot)
library(mapview)
library(tidyterra)


png(filename = "livestock_density_raster_map.png")
ggplot() +
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
dev.off()


library(mapview)
mapview(grid_1km_livestock, 
        zcol = "livestock_density",
        col.regions = viridisLite::viridis,
        layer.name = "Livestock Density")

library(tmap)
par(mfrow=c(1,1))
tm_shape(grid_1km_livestock) +
  tm_polygons("livestock_density",
              style = "fixed", breaks = c(0, 50, 100, 150, 200, 250),
              textNA = "No data",                 # Label for NA values
              colorNA = "lightgray",              # Color for NA values
              palette = "RdYlBu",
              title = "Livestock Density") +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(legend.outside = TRUE)

tm_shape(grid_1km_livestock) +
  tm_polygons("livestock_density",
              style = "quantile", 
              palette = "RdYlBu",
              title = "Livestock Density") +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(legend.outside = TRUE)

library(RColorBrewer)
par(mfrow = c(1,2), mar = c(2,2,3,2))

# 10km values-original raster
plot(livestock_10km, 
     main = "Original 10km Raster",
     col = rev(terrain.colors(100)), # High-contrast palette
     range = c(0, max(values(livestock_1km), na.rm = TRUE)), # Full value range
     axes = FALSE)

# 1km grid values
livestock_vect <- vect(grid_1km_livestock)
plot(livestock_vect, "livestock_density",
     main = "Grid Values",
     col = brewer.pal(9, "YlOrRd"), # ColorBrewer palette
     # breaks = "quantile", # Auto-breaks for better distribution
     axes = FALSE)

# plot(livestock_5km, main = "Original 5km Raster")
# plot(livestock_vect, "livestock_density", main = "Grid Values")

length(grid_1km_livestock$grid_id[is.na(grid_1km_livestock$livestock_density)])

# Create a map highlighting NA grids
ggplot(grid_1km_livestock) +
  geom_sf(aes(fill = is.na(livestock_density))) +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray70")) +
  labs(title = "Grids with Missing Livestock Data") +
  theme_bw()


write_csv(grid_1km_livestock, "csvs/livestockDensity_1km.csv")

###########################################

livestock <- st_read("~/Downloads/livestock2020.gpkg")
table(livestock$country)
livestock_scotland <- livestock %>%
  filter(country %in% "Scotland")

plot(livestock_scotland)


