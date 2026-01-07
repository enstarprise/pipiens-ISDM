library(ggplot2)
library(tidyr)
landcover_long <- pivot_longer(survey_df, 
                               cols = 15:24, 
                               names_to = "variable", 
                               values_to = "value")
