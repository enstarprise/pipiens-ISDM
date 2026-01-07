library(raster)
library(terra) 
library(sf)
library(tmap)
library(readr)

setwd("~/OneDrive - University of Glasgow/PhD/modelling/dataIntegration/ISDM")

scotland <- st_read("gridExtractions/scotlandBoundary/scotland_boundary.shp") %>%
  vect()

#--- LIVESTOCK DATA-----
livestockDensity <- rast("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/livestockDensity/GLW4-2020.D-DA.GLEAM3-ALL-LU.tif")
crs(livestockDensity) <- "EPSG:4326" # or whatever its original CRS is
st_crs(livestockDensity)
# plot(livestockDensity)

# crop and mask to scotland's shpaefile
livestockDensity <- crop(livestockDensity,  ext(scotland))
livestock_masked <- mask(livestockDensity, scotland)
livestockDensity <- raster(livestock_masked)




cattleDensity <- rast("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/livestockDensity/GLW4-2020.D-DA.CTL.tif")
crs(cattleDensity) <- "EPSG:4326"

cattleDensity <- crop(cattleDensity,  ext(scotland))
cattle_masked <- mask(cattleDensity, scotland)
cattleDensity <- raster(cattle_masked)

cattle_terra <- rast(cattleDensity)  # convert to spatraster

# resample to 5km : disaggregate
cattle <- disagg(cattle_terra, fact = 2, method = "bilinear")


grid_5km_cattle <- grid_clipped
grid_5km_cattle$cattle_density <- terra::extract(
  cattle , 
  grid_terra, 
  fun = mean,
  na.rm = TRUE
)[, 2]


par(mfrow=c(1,1))
png("cattle_tmap.png", width = 8, height = 6, units = "in", res = 300)
tm_shape(grid_5km_cattle) +
  tm_polygons("cattle_density",
              style = "quantile",
              palette = "RdYlBu",
              title = "Cattle Density") +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(legend.outside = TRUE)
dev.off()

sheepDensity <- rast("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/livestockDensity/GLW4-2020.D-DA.SHP.tif")

sheepDensity <- crop(sheepDensity,  ext(scotland))
sheep_masked <- mask(sheepDensity, scotland)
sheepDensity <- raster(sheep_masked)

sheep_terra <- rast(sheepDensity)  # convert to spatraster

# resample to 5km : disaggregate
sheep <- disagg(sheep_terra, fact = 2, method = "bilinear")


grid_5km_sheep <- grid_clipped
grid_5km_sheep$sheep_density <- terra::extract(
  sheep , 
  grid_terra, 
  fun = mean,
  na.rm = TRUE
)[, 2]

png("sheep_tmap.png", width = 8, height = 6, units = "in", res = 300)
tm_shape(grid_5km_sheep) +
  tm_polygons("sheep_density",
              style = "quantile",
              palette = "RdYlBu",
              title = "Sheep Density") +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(legend.outside = TRUE)
dev.off()

png("all_livestock_tmap.png", width = 8, height = 6, units = "in", res = 300)
tm_shape(grid_5km_livestock) +
  tm_polygons("livestock_density",
              style = "quantile",
              palette = "RdYlBu",
              title = "Livestock Density") +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(legend.outside = TRUE)
dev.off()

mapview(cattleDensity, 
        col.regions = viridisLite::viridis,
        layer.name = "Cattle Density") +
  mapview(sheepDensity,
          col.regions = viridisLite::viridis,
          layer.name = "Sheep Density")
m <-  mapview(livestockDensity, 
              col.regions = viridisLite::viridis,
              layer.name = "Livestock Density") +
  mapview(cattleDensity, 
             col.regions = viridisLite::viridis,
             layer.name = "Cattle Density") +
  mapview(sheepDensity,
          col.regions = viridisLite::viridis,
          layer.name = "Sheep Density")

# Save as HTML file
mapshot(m, url = "livestock_map.html")



# --- LANDCOVER CLASSIFICATIONS
landcover_map <- landcover_5km %>%
  mutate(
    wetland_km2 = Bog_km2 + Fen_km2 + Saltmarsh_km2,
    woodland_km2 = Deciduous_woodland_km2 + Coniferous_woodland_km2,
    saltwater_km2 = Saltwater_km2,
    freshwater_km2 = Freshwater_km2,
    grassland_heather_km2 = Heather_grassland_km2 + Heather_km2 +
      Improved_grassland_km2 + Acid_grassland_km2 +
      Neutral_grassland_km2 + Calcareous_grassland_km2,
    urban_km2 = Urban_km2, 
    suburban_km2 = Suburban_km2,
    arable_km2 = Arable_km2 
  )

landcover <- rast("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/landCover/UK_CEH/data/7727ce7d-531e-4d77-b756-5cc59ff016bd/ukregion-scotland.tif")
plot(landcover)

library(terra)

# Create a reclassification matrix
reclass_matrix <- matrix(c(
  0, 0, 0,   # Unnamed
  1, 1, 1,   # Deciduous_woodland
  2, 2, 2,   # Coniferous_woodland
  3, 3, 3,   # Arable
  4, 4, 4,   # Improved_grassland
  5, 5, 5,   # Neutral_grassland
  6, 6, 6,   # Calcareous_grassland
  7, 7, 7,   # Acid_grassland
  8, 8, 8,   # Fen
  9, 9, 9,   # Heather
  10, 10, 10, # Heather_grassland
  11, 11, 11, # Bog
  12, 12, 12, # Inland_rock
  13, 13, 13, # Saltwater
  14, 14, 14, # Freshwater
  15, 15, 15, # Supralittoral_rock
  16, 16, 16, # Supralittoral_sediment
  17, 17, 17, # Littoral_rock
  18, 18, 18, # Littoral_sediment
  19, 19, 19, # Saltmarsh
  20, 20, 20, # Urban
  21, 21, 21  # Suburban
), ncol = 3, byrow = TRUE)

landcover_categorical <- landcover
landcover_categorical <- classify(landcover, reclass_matrix)
landcover_categorical[landcover_categorical == 0] <- NA

levels(landcover_categorical) <- data.frame(
  ID = 0:21,
  landcover = c("Unnamed", "Deciduous_woodland", "Coniferous_woodland", "Arable",
                "Improved_grassland", "Neutral_grassland", "Calcareous_grassland",
                "Acid_grassland", "Fen", "Heather", "Heather_grassland", "Bog",
                "Inland_rock", "Saltwater", "Freshwater", "Supralittoral_rock",
                "Supralittoral_sediment", "Littoral_rock", "Littoral_sediment",
                "Saltmarsh", "Urban", "Suburban")
)
is.factor(landcover_categorical)


# Create a custom color palette with more distinct colors
custom_colors <- c(
  "transparent", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
  "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
  "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3",
  "#BEBADA", "#FB8072"
)

levels(landcover_categorical) <- data.frame(
  ID = 0:21,
  landcover = c("Unnamed", 
                "Deciduous_wood", "Coniferous_wood", "Arable",
                "Improved_grass", "Neutral_grass", "Calcareous_grass",
                "Acid_grass", "Heather", "Heather_grass", 
                "Rock", "Saltwater", "Freshwater", "Supralittoral_rock",
                "Supralittoral_sed", "Littoral_rock", "Littoral_sed",
                "Fen", "Bog","Saltmarsh", "Urban", "Suburban")
)

png(filename = "landcover_raster.png", width = 1000, height = 1052, bg = "transparent")
plot(landcover_categorical$landcover,
     col = custom_colors,  # custom colors
     axes = FALSE,
     box = FALSE,
     legend = TRUE,
     colNA = "transparent",
     mar = c(1, 1, 1, 15),
     plg = list(title = "Landcover\nTypes",
                title.cex = 1.5,
                size = 1.5,
                cex = 1.5))
dev.off()


library(RColorBrewer)
# Create a more distinct color palette
n_colors <- 22
qual_colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_colors)

plot(landcover_categorical$landcover,
     col = qual_colors,
     axes = FALSE,
     box = FALSE,
     legend = TRUE,
     colNA = "transparent",
     mar = c(1, 1, 1, 15),
     plg = list(title = "Landcover\nTypes",
                title.cex = 1.5,
                size = 1.5,
                cex = 1.5))

# group similar landcover types
levels(landcover_categorical) <- data.frame(
  ID = 0:21,
  landcover = c("Unnamed", 
                "Deciduous_wood", "Coniferous_wood", "Arable",
                "Improved_grass", "Neutral_grass", "Calcareous_grass",
                "Acid_grass", "Heather", "Heather_grass", 
                "Rock", "Saltwater", "Freshwater", "Supralittoral_rock",
                "Supralittoral_sed", "Littoral_rock", "Littoral_sed",
                "Fen", "Bog","Saltmarsh", "Urban", "Suburban")  # Shorter names
)

png(filename = "landcover_raster.png", width = 1000, height = 1052, bg = "transparent")
plot(landcover_categorical$landcover,
     axes = FALSE,
     box = FALSE,
     legend = TRUE,
     colNA = "transparent",
     mar = c(1, 1, 1, 15),
     plg = list(title = "Landcover\nTypes",
                title.cex = 1.5,
                size = 1.5,
                cex = 1.5))
dev.off()







###### group by modelling grouping
library(tidyterra)
library(ggplot2)
landcover_categorical <- landcover


# Create a reclassification matrix that groups land cover types
reclass_matrix <- matrix(c(
  0, 0, 9,    # Unnamed -> other
  1, 1, 1,    # Deciduous_woodland -> woodland
  2, 2, 1,    # Coniferous_woodland -> woodland
  3, 3, 2,    # Arable -> arable
  4, 4, 3,    # Improved_grassland -> grassland_heather
  5, 5, 3,    # Neutral_grassland -> grassland_heather
  6, 6, 3,    # Calcareous_grassland -> grassland_heather
  7, 7, 3,    # Acid_grassland -> grassland_heather
  8, 8, 4,    # Fen -> wetland
  9, 9, 3,    # Heather -> grassland_heather
  10, 10, 3,  # Heather_grassland -> grassland_heather
  11, 11, 4,  # Bog -> wetland
  12, 12, 9,  # Inland_rock -> other
  13, 13, 5,  # Saltwater -> saltwater
  14, 14, 6,  # Freshwater -> freshwater
  15, 15, 9,  # Supralittoral_rock -> other
  16, 16, 9,  # Supralittoral_sediment -> other
  17, 17, 9,  # Littoral_rock -> other
  18, 18, 9,  # Littoral_sediment -> other
  19, 19, 4,  # Saltmarsh -> wetland
  20, 20, 7,  # Urban -> urban
  21, 21, 8   # Suburban -> suburban
), ncol = 3, byrow = TRUE)

landcover_categorical <- classify(landcover, reclass_matrix)
landcover_categorical <- landcover_categorical[[1]]
names(landcover_categorical) <- "landcover"

# Convert to categorical/factor raster
levels(landcover_categorical) <- data.frame(
  value = 1:9,
  landcover = c("Woodland", "Arable", "Grassland/Heather", "Wetland", 
                "Saltwater", "Freshwater", "Urban", "Suburban", "Other")
)

landcover_colors <- c(
  "Woodland" = "#228B22",
  "Arable" = "#DAA520",
  "Grassland/Heather" = "#90EE90",
  "Wetland" = "#4682B4",
  "Saltwater" = "#000080",
  "Freshwater" = "#87CEEB",
  "Urban" = "#8B0000",
  "Suburban" = "#CD5C5C",
  "Other" = "#D3D3D3"
)

ggplot() +
  geom_spatraster(data = landcover_categorical) +
  scale_fill_manual(
    values = landcover_colors,
    na.translate = FALSE,
    name = "Land Cover"
  ) +
  theme_minimal() +
  labs(title = "Aggregated Land Cover Types") +
  coord_sf()



# --- ELEVATION ---
scotland <- st_read("gridExtractions/scotlandBoundary/scotland_boundary.shp") %>%
  vect()
grid_clipped <- st_read("grids/grid_clipped.gpkg")
grid_clipped <- st_transform(grid_clipped, st_crs(scotland))

elevation <- rast("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/elevation_DEM/scotland_elevation.tif")
st_crs(elevation) == st_crs(grid_clipped)
elevation_5km <- read_csv("gridExtractions/elevation_5km.csv")

grid_with_elevation <- grid_clipped %>%
  left_join(elevation_5km , by = "grid_id")  # replace with your actual ID column name

elevation_map <- tm_shape(grid_with_elevation) +
  tm_polygons("elevation",
              style = "quantile",
              palette = "RdYlBu",
              title = "Elevation (m)") +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(legend.outside = TRUE)

elevation_map

png(filename = "elevation_raster.png", width = 1000, height = 1052, bg = "transparent")
plot(elevation,
     axes = FALSE,    # Remove axis lines
     box = FALSE)
dev.off()

library(tmap)

tm_shape(elevation) +
  tm_raster() +
  tm_layout(frame = FALSE,   
            bg.color = NA) 



#  --- AVERAGE TEMPERATURE -----
#> must be converted to a raster? like the livestock data as above...
#> done in the relevant gridExtrations .R file







# --- AVERAGE RAINFALL -----
#> done in the relevant gridExtrations .R file





