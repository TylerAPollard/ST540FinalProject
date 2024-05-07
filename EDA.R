#### ST540 Final Project EDA 
## Programmed by: Tyler Pollard, Hanan Ali, Rachel Hardy
## Date Created: 4 April 2024
## Date Modified: 

# Load Libraries ----
library(data.table)
library(MASS)
library(rjags)
library(plyr)
library(GGally)
library(tidyverse)
library(caret)
library(tictoc)
library(brms)
library(BayesFactor)

# Read in Global Bleaching data ----
bleaching_data <- fread("global_bleaching_environmental.csv", 
                        na.strings = c("", "NA", "nd"))

## Check sample sizes from paper
sum(!is.na(bleaching_data$Percent_Bleaching))

# Initial EDA ==============================================================
## Check missing data values for response ----
nonmissing_byYear <- ddply(bleaching_data, .(Date_Year), summarise,
                           Total = length(Percent_Bleaching),
                           Non_Missing = sum(!is.na(Percent_Bleaching)),
                           Percent_NotMissing = round(Non_Missing/Total*100, 2)) |>
  arrange(Date_Year)

nonmissing_bySource <- ddply(bleaching_data, .(Data_Source), summarise,
                             Total = length(Percent_Bleaching),
                             Non_Missing = sum(!is.na(Percent_Bleaching)),
                             Percent_NotMissing = round(Non_Missing/Total*100, 2)) |>
  arrange(Data_Source)

nonmissing_byYearSource <- ddply(bleaching_data, .(Data_Source, Date_Year), summarise,
                                 Total = length(Percent_Bleaching),
                                 Non_Missing = sum(!is.na(Percent_Bleaching)),
                                 Percent_NotMissing = round(Non_Missing/Total*100, 2)) |>
  arrange(Data_Source, Date_Year)

ggplot(data = Reef_Check_df2) +
  geom_col(aes(x = Date_Year, y = Percent_Bleaching)) +
  theme_bw()
## Remove Nuryana and Setiawan due to little data

## Examine each data source ----
### AGRRA ----
AGRRA_df <- bleaching_data |>
  filter(Date_Year >= 1998) |>
  filter(Data_Source == "AGRRA") |>
  arrange(Site_ID, Sample_ID, Date_Year, Date_Month, Date_Day)
unique(AGRRA_df$Percent_Bleaching)
AGRRA_depth_df <- ddply(AGRRA_df, .(Site_ID, Date_Year, Date_Month, Date_Day), summarize,
                        Site_Variance = var(Percent_Bleaching),
                        Observations_N = length(Percent_Bleaching))
AGRRA_depth_vary_df <- AGRRA_depth_df |> filter(Site_Variance != 0)

AGRRA_df2 <- AGRRA_df |>
  filter(Site_ID %in% AGRRA_depth_vary_df$Site_ID)
AGRRA_depth_df2 <- ddply(AGRRA_df2, .(Site_ID, Date_Day, Date_Month, Date_Year), summarize,
                         Site_Variance = var(Percent_Bleaching),
                         Observations_N = length(Percent_Bleaching))
AGRRA_depth_vary_df2 <- AGRRA_depth_df2 |> filter(Site_Variance != 0)

AGRRA_df3 <- AGRRA_df |>
  filter(Site_ID %in% AGRRA_depth_vary_df2$Site_ID) |>
  arrange(Depth_m)

plot(AGRRA_df3$Depth_m, AGRRA_df3$Percent_Bleaching)
cor(AGRRA_df3$Depth_m, AGRRA_df3$Percent_Bleaching, use = "complete.obs")
## Continuous Percent Bleaching values
## Multiple samples were taken at same location same day for various depths
## Only variables that varies across observations from same site and day is depth
## Possibly aggreate?

### Donner ----
Donner_df <- bleaching_data |>
  filter(Date_Year >= 1998) |>
  filter(Data_Source == "Donner") |>
  arrange(Date_Year, Date_Month, Date_Day, Site_ID, Sample_ID)
unique(Donner_df$Percent_Bleaching)
hist(Donner_df$Percent_Bleaching, breaks = 20)
## Large number of category average values 0, 5.5, 30.5, 75

sum(is.na(Donner_df$Depth_m))
# About 1/4 of data has missing depths

Donner_df2 <- Donner_df |>
  filter(complete.cases(Percent_Bleaching)) |>
  arrange(Site_ID, Sample_ID,Date_Year, Date_Month, Date_Day)
Donner_df2 <- ddply(Donner_df2, .(Site_ID), summarize, .drop = FALSE,
                    Avg_Percent_Bleaching = mean(Percent_Bleaching))

Donner_depth_df <- ddply(Donner_df, .(Site_ID), summarize,
                         Site_Variance = var(Percent_Bleaching),
                         Observations_N = length(Percent_Bleaching))
Donner_depth_vary_df <- Donner_depth_df |> filter(Site_Variance != 0) |> filter(!is.na(Site_Variance))
Donner_df2 <- Donner_df |>
  filter(Site_ID %in% Donner_depth_vary_df$Site_ID)
## Tentatively keep Donner_df as is, but need to address missing depth values

### FRRP ----
FRRP_df <- bleaching_data |>
  filter(Date_Year >= 1998) |>
  filter(Data_Source == "FRRP") |>
  arrange(Site_ID, Sample_ID, Date_Year, Date_Month, Date_Day)
unique(FRRP_df$Percent_Bleaching)
unique(FRRP_df$City_Town_Name)
## 5 different counties with varying locations within county
length(unique(FRRP_df$Site_ID))
hist(FRRP_df$Percent_Bleaching, breaks = 20)
## Keep FRRP data as is. Good data

### Kumagai ----
Kumagai_df <- bleaching_data |>
  filter(Date_Year >= 1998) |>
  filter(Data_Source == "Kumagai") |>
  arrange(Site_ID, Sample_ID, Date_Year, Date_Month, Date_Day)
unique(Kumagai_df$Percent_Bleaching)
## Notice that percent bleaching is just the middle of each rating
### ie 0 = 0, 0-10 = 5.5, 11-50 = 30.5, 50-100 = 75
## May need to exclude due to categorical responses


### McClanahan ----
McClanahan_df <- bleaching_data |>
  filter(Date_Year >= 1998) |>
  filter(Data_Source == "McClanahan") |>
  arrange(Site_ID, Sample_ID, Date_Year, Date_Month, Date_Day)
unique(McClanahan_df$Percent_Bleaching)
unique(McClanahan_df$City_Town_Name)
unique(McClanahan_df$Site_ID)
## Data is for only 2016
## Data is complete and appears complete
## Keep McClanahan as is
hist(McClanahan_df$Percent_Bleaching, breaks = 25)

### Reef_Check ----
Reef_Check_df <- bleaching_data |>
  filter(Date_Year >= 1998) |>
  filter(Data_Source == "Reef_Check") |>
  arrange(Site_ID, Date_Year, Date_Month, Date_Day, Sample_ID)
unique(Reef_Check_df$Percent_Bleaching)
## Same Sample ID has multiple readings by substrate Name
Reef_Check_Substrate_df <- ddply(Reef_Check_df, .(Site_ID, Sample_ID), summarize,
                                 Substrate_Var = var(Percent_Bleaching),
                                 Observations_N = length(Percent_Bleaching))
Reef_Check_Substrate_vary_df <- Reef_Check_Substrate_df |>
  filter(Substrate_Var != 0) |>
  filter(!is.na(Substrate_Var))
## Percent_Bleaching does not vary by sample ID but it does for Percent_Cover
Reef_Check_df2 <- Reef_Check_df |>
  distinct(Site_ID, Sample_ID, .keep_all = TRUE) |>
  filter(complete.cases(Percent_Bleaching))
hist(Reef_Check_df2$Percent_Bleaching, breaks = 100)
sum(Reef_Check_df2$Percent_Bleaching == 0)
unique(Reef_Check_df2$Percent_Bleaching)

# Split by 0 or not
Reef_Check_df3_zero <- Reef_Check_df2 |> filter(Percent_Bleaching == 0)

Reef_Check_df3_nonzero <- Reef_Check_df2 |> filter(Percent_Bleaching != 0)
hist(Reef_Check_df3_nonzero$Percent_Bleaching, breaks = 20)
## Agregated data appears to be good as is once NAs are removed

### Safaie ----
Safaie_df <- bleaching_data |>
  filter(Date_Year >= 1998) |>
  filter(Data_Source == "Safaie") |>
  arrange(Site_ID, Date_Year, Date_Month, Date_Day, Sample_ID)
unique(Safaie_df$Percent_Bleaching)
## Percent_Bleaching appears to be categorically values
## Recommend remove data source

## Plot Map ----
map_complete_data <- bleaching_data |>
  filter(!is.na(Percent_Bleaching)) |>
  filter(Data_Source %in% c("Reef_Check")) |>
  arrange(Data_Source, Site_ID, Date_Year, Date_Month, Date_Day, Sample_ID)
world_coordinates <- map_data("world")
ggplot() +
  # geom_map() function takes world coordinates
  # as input to plot world map
  geom_map(
    data = world_coordinates, map = world_coordinates,
    aes(x = long, y = lat, map_id = region)
  ) +
  geom_point(
    data = map_complete_data,
    aes(x = Longitude_Degrees, y = Latitude_Degrees,
        color = Percent_Bleaching)
  ) +
  scale_color_continuous(low = "green", high = "red") +
  theme_bw()

## Notice that percent bleaching is just the middle of each rating
### ie 0 = 0, 0-10 = 5.5, 11-50 = 30.5, 50-100 = 75

## TAKEAWAYS ----
## ID Variable Hierarchy 
## Realm_Name > Country_Name > Ecoregion_Name > State_Island_Province_Name > 
## City_Town_Name > Site _Name (if applicable)
## SiteID has same lat/lon combinations
## Only depth, Sample_ID, and Month/Day/Year varies with siteID
## Same Distance_to_shore, exposure,...,
## Temperatures may change for Site_ID by date
## We will choose to only look at data from 2003 and on because 
## 1. Few data before 1998
## 2. Data from 1998-2002 had a lot of missing data between(42-62% non-missing)
##    - Don't feel comfortbale trusting these data sources
## 3. Data >= 2003 had at least 1000 observations with >82% non-missing


## FINAL DATA SET ----
### Data Sources to keep:
### Reef_Check
final_data1 <- bleaching_data |>
  filter(!is.na(Percent_Bleaching)) |>
  filter(Data_Source == "FRRP") |>
  distinct(Site_ID, Sample_ID, .keep_all = TRUE)

### Filter down to variables of interest ----
final_data2 <- final_data1 |> 
  select(
    # For ordering
    Date,
    # Covariates
    Latitude_Degrees,
    Longitude_Degrees,
    Distance_to_Shore,
    Exposure,
    Turbidity,
    Cyclone_Frequency,
    Date_Year,
    Depth_m,
    ClimSST,
    SSTA,
    SSTA_DHW,
    TSA,
    TSA_DHW,
    Windspeed,
    # Response
    Percent_Bleaching
  ) |>
  filter(Date_Year >= 2003) |>
  arrange(Date)

# FINAL DATA ===============================================================
## Clean data for all Models----
### Original Data ----
bleaching_data <- fread("global_bleaching_environmental.csv", 
                        na.strings = c("", "NA", "nd"))

### Filter Data 1 ----
## Filter to only complete Percent Bleaching and FRRP data set which
## is the 
final_data1 <- bleaching_data |>
  filter(!is.na(Percent_Bleaching)) |>
  filter(Data_Source == "FRRP") |>
  distinct(Site_ID, Sample_ID, .keep_all = TRUE)

### Filter Data 2 ----
## Remove unwanted variables like temperature statistic columns
## and arrange by date for viewing purposes
final_data2 <- final_data1 |> 
  select(
    # For ordering
    Date,
    City_Town_Name,
    # Covariates
    Latitude_Degrees,
    Longitude_Degrees,
    Distance_to_Shore,
    Exposure,
    Turbidity,
    Cyclone_Frequency,
    Date_Year,
    Date_Month,
    Depth_m,
    ClimSST,
    SSTA,
    SSTA_DHW,
    TSA,
    TSA_DHW,
    Windspeed,
    # Response
    Percent_Bleaching
  ) |>
  arrange(Date)

### Filter Data 3 ----
# Remove rows with missing predictors values
final_data3 <- final_data2[complete.cases(final_data2)]
final_data3$Percent_Bleaching_Open <- 
  ifelse(final_data3$Percent_Bleaching == 0, 0.001, 
         ifelse(final_data3$Percent_Bleaching == 100, 99.999, final_data3$Percent_Bleaching))
final_data3$Percent_Bleaching_log <- log(final_data3$Percent_Bleaching_Open)

### Create training/test index vector ----
## for CV and testing model performance and fit
set.seed(52)
trainIndex <- createDataPartition(final_data3$Percent_Bleaching,
                                  p = 0.75,
                                  list = FALSE)
trainIndex <- as.vector(trainIndex)

## EDA for Cleaned Data ----
byYearfun <- function(x){
  mean(x)
}
byYear_EDA <- ddply(final_data3, .(Date_Year), summarize,
                    ClimSST = byYearfun(ClimSST),
                    MedClimSST = median(ClimSST),
                    SSTA = byYearfun(SSTA),
                    SSTA_DHW = byYearfun(SSTA_DHW),
                    TSA = byYearfun(TSA),
                    TSA_DHW = byYearfun(TSA_DHW),
                    Percent_Bleaching = byYearfun(Percent_Bleaching))

temp_data <- final_data3 |> select(
  Distance_to_Shore,
  Cyclone_Frequency,
  Date_Year,
  TSA,
  TSA_DHW,
  Percent_Bleaching
)
ggpairs(temp_data)

temp_data <- final_data3 |> select(
  ClimSST,
  SSTA,
  SSTA_DHW,
  TSA,
  TSA_DHW,
  Percent_Bleaching
)
ggpairs(temp_data)

temp_data2U <- final_data3 |>
  filter(TSA_DHW >= 20)
temp_data2L <- final_data3 |>
  filter(TSA_DHW < 20)


ggplot() +
  geom_density(data = final_data3,
               aes(x = Percent_Bleaching, color = factor(Date_Year)), linewidth = 1) +
  theme_bw()

ggplot() +
  geom_density(data = final_data3,
               aes(x = log(Percent_Bleaching), color = factor(Date_Year)), linewidth = 1) +
  theme_bw()

unique(final_data3$City_Town_Name)
ddply(final_data3, .(Date_Year), summarize,
      Obs = length(Percent_Bleaching),
      Mean = mean(Percent_Bleaching))

trainIndex
testIndex <- 1:2394
testIndex <- testIndex[-trainIndex]

final_data_pred2 <- final_data3 |>
  select(Latitude_Degrees, Longitude_Degrees) |>
  slice(testIndex)
final_data_pred2$Percent_Bleaching_Mean <- YppdMean2
final_data_pred2$Percent_Bleaching_Meed <- YppdMedian2

world_coordinates <- map_data("county") 
ggplot() + 
  # geom_map() function takes world coordinates  
  # as input to plot world map 
  geom_map( 
    data = world_coordinates, map = world_coordinates, 
    aes(x = long, y = lat, map_id = region) 
  ) + 
  geom_point(
    data = final_data_pred2,
    aes(x = Longitude_Degrees, y = Latitude_Degrees, 
        color = Percent_Bleaching_Mean)
  ) +
  xlim(c(-85,-77.5)) +
  ylim(c(23,32.5)) +
  scale_color_continuous(low = "green", high = "red") +
  #facet_wrap(vars(City_Town_Name)) +
  theme_bw()






