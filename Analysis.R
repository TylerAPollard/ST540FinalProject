#### ST540 Final Project
## Programmed by: Tyler Pollard, Hanan Ali, Rachel Hardy
## Date Created: 4 April 2024
## Date Modified: 

# Load Libraries ----
library(data.table)
library(MASS)
library(rjags)
library(plyr)
library(tidyverse)

# Read in Global Bleaching data ----
bleaching_data <- fread("global_bleaching_environmental.csv", 
                        na.strings = c("", "NA", "nd"))

missing_df <- ddply(bleaching_data, .(Date_Year), summarise,
                    Total = length(Percent_Bleaching),
                    Non_Missing = sum(!is.na(Percent_Bleaching))) |>
  arrange(Date_Year)


df_2002 <- bleaching_data |> filter(Date_Year == 2002)
sum(!is.na(df_2002$Percent_Bleaching))

df_2002_comments <- bleaching_data |> 
  filter(Date_Year == 2002) |>
  select(Percent_Bleaching, Bleaching_Comments)

filtered_df <- bleaching_data |>
  filter(complete.cases(Percent_Bleaching)) |>
  filter(Date_Year >= 2000)

ggplot(data = filtered_df) +
  geom_histogram(aes(x = Percent_Bleaching))

sum(filtered_df$Percent_Bleaching == 0)

small_filtered_df <- head(filtered_df, 100)

sum(!complete.cases(filtered_df$Distance_to_Shore))
sum(!complete.cases(filtered_df$Exposure))
sum(!complete.cases(filtered_df$Turbidity))
sum(!complete.cases(filtered_df$Cyclone_Frequency))

### Aded comment

# Analysis ####################################################################################################
## Load data ====
### Regression ----









