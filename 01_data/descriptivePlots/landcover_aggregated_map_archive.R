

landcover <- rast("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/landCover/UK_CEH/data/7727ce7d-531e-4d77-b756-5cc59ff016bd/ukregion-scotland.tif")


# Create a reclassification matrix that groups land cover types
# reclass_matrix <- matrix(c(
#   0, 0, 0,   # Unnamed
#   1, 1, 1,   # Deciduous_woodland
#   2, 2, 2,   # Coniferous_woodland
#   3, 3, 3,   # Arable
#   4, 4, 4,   # Improved_grassland
#   5, 5, 5,   # Neutral_grassland
#   6, 6, 6,   # Calcareous_grassland
#   7, 7, 7,   # Acid_grassland
#   8, 8, 8,   # Fen
#   9, 9, 9,   # Heather
#   10, 10, 10, # Heather_grassland
#   11, 11, 11, # Bog
#   12, 12, 12, # Inland_rock
#   13, 13, 13, # Saltwater
#   14, 14, 14, # Freshwater
#   15, 15, 15, # Supralittoral_rock
#   16, 16, 16, # Supralittoral_sediment
#   17, 17, 17, # Littoral_rock
#   18, 18, 18, # Littoral_sediment
#   19, 19, 19, # Saltmarsh
#   20, 20, 20, # Urban
#   21, 21, 21  # Suburban
# ), ncol = 3, byrow = TRUE)

# reclass_matrix <- matrix(c(
#   0, 0, NA,   # Unnamed -> NA
#   1, 1, 1,    # Deciduous_woodland -> woodland
#   2, 2, 1,    # Coniferous_woodland -> woodland
#   3, 3, 2,    # Arable -> arable
#   4, 4, 3,    # Improved_grassland -> grassland_heather
#   5, 5, 3,    # Neutral_grassland -> grassland_heather
#   6, 6, 3,    # Calcareous_grassland -> grassland_heather
#   7, 7, 3,    # Acid_grassland -> grassland_heather
#   8, 8, 4,    # Fen -> wetland
#   9, 9, 3,    # Heather -> grassland_heather
#   10, 10, 3,  # Heather_grassland -> grassland_heather
#   11, 11, 4,  # Bog -> wetland
#   12, 12, 9, # Inland_rock -> Rock/sediment
#   13, 13, 5,  # Saltwater -> saltwater
#   14, 14, 6,  # Freshwater -> freshwater
#   15, 15, 9, # Supralittoral_rock -> Rock/sediment
#   16, 16, 9, # Supralittoral_sediment -> Rock/sediment
#   17, 17, 9, # Littoral_rock -> Rock/sediment
#   18, 18, 9, # Littoral_sediment -> Rock/sediment
#   19, 19, 4,  # Saltmarsh -> wetland
#   20, 20, 7,  # Urban -> urban
#   21, 21, 8   # Suburban -> suburban
# ), ncol = 3, byrow = TRUE)

reclass_matrix <- matrix(c(
  0, 0, 9,   # Unnamed -> Rock/sediment (CHANGED from NA to 9)
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
  12, 12, 9,  # Inland_rock -> Rock/sediment
  13, 13, 5,  # Saltwater -> saltwater
  14, 14, 6,  # Freshwater -> freshwater
  15, 15, 9,  # Supralittoral_rock -> Rock/sediment
  16, 16, 9,  # Supralittoral_sediment -> Rock/sediment
  17, 17, 9,  # Littoral_rock -> Rock/sediment
  18, 18, 9,  # Littoral_sediment -> Rock/sediment
  19, 19, 4,  # Saltmarsh -> wetland
  20, 20, 7,  # Urban -> urban
  21, 21, 8   # Suburban -> suburban
), ncol = 3, byrow = TRUE)

# CORRECTED: Use landcover_categorical consistently
landcover_categorical <- landcover
landcover_categorical <- landcover_categorical[[1]]
landcover_categorical <- classify(landcover_categorical, reclass_matrix)  # Fixed this line

landcover_categorical[landcover_categorical == 0] <- NA

# Create a nice color palette and labels for plotting
# landcover_colours <- c(
#   "1" = "#228B22",  # woodland - forest green
#   "2" = "#DAA520",  # arable - goldenrod
#   "3" = "#90EE90",  # grassland_heather - light green
#   "4" = "#4682B4",  # wetland - steel blue
#   "5" = "#000080",  # saltwater - navy
#   "6" = "#87CEEB",  # freshwater - sky blue
#   "7" = "#8B0000",  # urban - dark red
#   "8" = "#CD5C5C",  # suburban - indian red
#   "9" = "#A9A9A9"   # Rock/sediment - dark gray
# )

# landcover_colours <- c(
#   "1" = "#228B22",  # woodland - forest green
#   "2" = "#DAA520",  # arable - goldenrod
#   "3" = "#90EE90",  # grassland_heather - light green
#   "4" = "#8B008B",  # wetland - dark magenta
#   "5" = "#191970",  # saltwater - midnight blue
#   "6" = "#1E90FF",  # freshwater - dodger blue
#   "7" = "#8B0000",  # urban - dark red
#   "8" = "#CD5C5C",  # suburban - indian red
#   "9" = "#696969"   # Rock/sediment - dim gray
# )

landcover_labels <- c(
  "1" = "Woodland",
  "2" = "Arable",
  "3" = "Grassland/Heather",
  "4" = "Wetland (bog,\nsaltmarsh, fen)",
  "5" = "Saltwater",
  "6" = "Freshwater",
  "7" = "Urban",
  "8" = "Suburban",
  "9" = "Rock/sediment"
)


custom_colors <- c("#228B22",  # woodland - forest green
                    "#DAA520",  # arable - goldenrod
                    "#90EE90",  # grassland_heather - light green
                    "#8B008B",  # wetland - dark magenta
                    "#191970",  # saltwater - midnight blue
                    "#1E90FF",  # freshwater - dodger blue
                    "#8B0000",  # urban - dark red
                    "#CD5C5C",  # suburban - indian red
                    "#696969")   # Rock/sediment - dim gray)
                    
levels(landcover_categorical) <- data.frame(
  ID = 0:9,
  landcover = c("Unnamed", "Woodland", "Arable",
                "Grassland", "Wetland",
                 "Saltwater", "Freshwater", "Urban", "Suburban",
                "Rock/sediment")
)

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


















# Create a single plot with proper margins for legend
par(mar = c(5, 4, 4, 8))  # Extra space on right for legend
par(mar = c(5, 4, 4, 2) + 0.1)  # Default margins

png("map_landcover.png", width = 10, height = 12, units = "in", res = 300)
par(mar = c(5, 4, 4, 8))  # Set margins for legend space
plot(landcover_categorical, 
     col = landcover_colours,
     main = "Aggregated Land Cover Types",
     legend = FALSE)

legend("right", 
       legend = landcover_labels, 
       fill = landcover_colours,
       title = "Land Cover",
       inset = c(-0.2, 0),  # Move legend left
       xpd = TRUE)  # Allow plotting outside main area

dev.off()


png("map_landcover_legend.png", width = 4, height = 6, units = "in", res = 300)
plot.new()
legend("right", 
       legend = landcover_labels, 
       fill = landcover_colours,
       title = "Land Cover",
       inset = c(-0.2, 0),  # Move legend left
       xpd = TRUE)  # Allow plotting outside main area
dev.off()


# CROP TO ISLE OF ARRAN

mask <- st_read("~/OneDrive - University of Glasgow/PhD/modelling/dataIntegration/MARS/arran_mask.shp")
mask <- st_transform(mask, st_crs(landcover_categorical))

# Crop and mask the landcover data
landcover_cropped <- crop(landcover_categorical, mask)
landcover_masked <- mask(landcover_cropped, mask)

# Create a single plot with proper margins for legend
par(mar = c(1, 1, 1, 0))  # Extra space on right for legend

# Save the main map
png("map.ARRAN_landcover.png", width = 10, height = 12, units = "in", res = 300)
par(mar = c(5, 4, 4, 8))  # Set margins for legend space
plot(landcover_masked, 
     col = landcover_colours,
     main = "Aggregated Land Cover Types - Arran",
     legend = FALSE)
legend("right", 
       legend = landcover_labels, 
       fill = landcover_colours,
       title = "Land Cover",
       inset = c(-0.2, 0),
       xpd = TRUE)
dev.off()

# Save just the legend separately
png("map.ARRAN_landcover_legend.png", width = 4, height = 6, units = "in", res = 300)
plot.new()
legend("center",
       legend = landcover_labels,
       fill = landcover_colours,
       title = "Land Cover")
dev.off()





# Apply the reclassification
landcover_categorical <- landcover
landcover_categorical <- classify(landcover, reclass_matrix)
landcover_categorical[landcover_categorical == 0] <- NA

# Set the factor levels for your AGGREGATED categories
levels(landcover_categorical) <- data.frame(
  ID = 1:9,
  landcover = c("Woodland", "Arable", "Grassland/Heather", "Wetland",
                "Saltwater", "Freshwater", "Urban", "Suburban", "Rock/sediment")
)

# is.factor(landcover_categorical)

# Create a custom color palette for your 9 aggregated categories
custom_colors <- c(
  "#228B22",  # Woodland - forest green
  "#DAA520",  # Arable - goldenrod
  "#90EE90",  # Grassland/Heather - light green
  "#8B4513",  # Wetland - saddle brown
  "#000080",  # Saltwater - navy blue
  "#00BFFF",  # Freshwater - deep sky blue
  "#8B0000",  # Urban - dark red
  "#CD5C5C",  # Suburban - indian red
  "#A9A9A9"   # Rock/sediment - dark gray
)

# Crop and mask to Arran
landcover_masked <- mask(crop(landcover_categorical, mask), mask)

# Plot with proper settings
png(filename = "ARRAN_landcover.png", width = 1000, height = 1052, bg = "transparent")
plot(landcover_masked$landcover,
     col = custom_colors,
     axes = FALSE,
     box = FALSE,
     legend = TRUE,
     colNA = "transparent",
     mar = c(1, 1, 1, 15),  # Extra space on right for legend
     plg = list(title = "Land Cover\nTypes",
                title.cex = 1.5,
                size = 1.5,
                cex = 1.5,
                x = "right",
                inset = -0.2))  # Adjust legend position if needed
dev.off()






























