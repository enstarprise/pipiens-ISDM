## combine new data with the old

setwd("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/modelling/dataIntegration/ISDM")

library(readr)
library(readxl)
library(dplyr)
library(sf)
library(writexl)

####### LOAD THE DATASHEETS

arran_data <- read_csv("01_data/csvs/Arran-sampling-2025.csv")
cs_new <- read_excel("01_data/csvs/Citizen-science-reports-new.xlsx")
survey_2025summer <- read_excel("01_data/csvs/Surveillance for modelling work.xlsx")


survey_2025summer$Setup_date <- as.Date(survey_2025summer$Setup_date, 
                                 format = "%d/%m/%Y")

survey_2025summer$Collection_date <- as.Date(survey_2025summer$Collection_date, 
                                      format = "%d/%m/%Y")

survey_2025summer$Year <- format(survey_2025summer$Collection_date, "%Y")

survey_2025summer$Month_collection <- format(survey_2025summer$Collection_date, "%m")  # Extracts year as character

namechange <- c("Possil marsh")
survey_2025summer$Site_name<- replace(survey_2025summer$Site_name, survey_2025summer$Site_name
                               %in% namechange, "Possil Marsh")
survey_2025summer <- survey_2025summer %>%
  mutate(Culex = Cx.pipiens.F + Cx.pipiens.M)

# --- citizen science ---
# cs_new <- cs_new %>% filter(Verified_mosquito == "Yes" & Species == "pipiens" & Stage == "Adult")

cs_new$Date_found <- as.Date(cs_new$Date_found, 
                                          format = "%d/%m/%Y")

cs_new$Year <- format(cs_new$Date_found, "%Y")
cs_new$Month <- format(cs_new$Date_found, "%m")  # Extracts year as character



#### COMBINE SURVEY DATA WITH 2025 SURVEYING
survey_data <- read_csv("01_data/csvs/merged_data_sets_28Feb2025.csv")

combined_data <- bind_rows(survey_data, survey_2025summer)


### SAVE THE DATASHEETS
write_csv(cs_new, "cs_df.csv")
write_csv(combined_data, "survey_df.csv")
