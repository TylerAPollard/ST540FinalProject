#### ST540 Final Project
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
library(ggplot2)

# Read in Global Bleaching data ----
bleaching_data <- fread("global_bleaching_environmental.csv", 
                        na.strings = c("", "NA", "nd"))

## Check sample sizes from paper
sum(!is.na(bleaching_data$Percent_Bleaching))

# EDA ===============================================================================
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
  filter(Data_Source == "Reef_Check") |>
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


ggpairs(data = final_data1 |> 
          select(ClimSST,
                 Temperature_Kelvin,
                 Temperature_Mean,
                 Temperature_Minimum,
                 Temperature_Maximum))

# Analysis ==========================================================================
## Load data ----
bleaching_data <- fread("global_bleaching_environmental.csv", 
                        na.strings = c("", "NA", "nd"))

final_data1 <- bleaching_data |>
  filter(!is.na(Percent_Bleaching)) |>
  filter(Data_Source == "Reef_Check") |>
  distinct(Site_ID, Sample_ID, .keep_all = TRUE)

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

# The thingy Tyler told me to add:
missingValues <- c()
for(i in 1:length(colnames(final_data2))){
  missingValues[i] <- sum(is.na(final_data2[[i]]))
}
names(missingValues) <- colnames(final_data2)
missingValues

# Remove rows with missing predictors values
final_data3 <- final_data2[complete.cases(final_data2)]

## Model 1: Simple Linear Model (Rachel) ----
## Modeled with Uninformative Gaussian Priors
final_rachel <- final_data3
final_rachel <- na.omit(final_rachel)
Y1 <- final_rachel$Percent_Bleaching
X1 <- subset(final_rachel, select = -c(Date, Date_Year, Exposure, Percent_Bleaching))
X1 <- subset(X1, select = -c(Turbidity, SSTA, TSA, Windspeed))
X1 <- as.matrix(X1)
X1 <- scale(X1)

n1 <- length(Y1)
p1 <- ncol(X1)

data1   <- list(Y=Y1,X=X1,n=n1,p=p1)
params1 <- c("alpha", "beta")

burn     <- 5000
n.iter   <- 20000
thin     <- 5

# After running the model below, I removed the following insignificant predictors:
# Turbidity, SSTA, TSA, and Windspeed.
# After removal, the model was run again.

model_string1 <- textConnection("model{

   # Likelihood
    for(i in 1:n){
      Y[i]   ~ dnorm(mu[i],tau) 
      mu[i] <- alpha + inprod(X[i,],beta[]) 
    } 
    
    # Priors  
    for(j in 1:p){ 
      beta[j] ~ dnorm(0,0.01) 
    } 
      
    alpha ~  dnorm(0, 0.01) 
    tau   ~  dgamma(0.1, 0.1) 
      
}")

model1 <- jags.model(model_string1, data=data1, n.chains=2, quiet=TRUE)
update(model1, burn, progress.bar="none")
samples1 <- coda.samples(model1, variable.names=params1, n.iter=n.iter, n.thin=thin, progress.bar="none")

summary1 <- summary(samples1)
summary1

par(mar=c(1,1,1,1)) 
# ^This is because of the error: "figure margins too large"
# Unfortunately it makes the plots kinda wonky but at least we can see them.
plot(samples1)

stats1 <- summary1$statistics
rownames(stats1) <- c("Intercept", "Latitude", "Longitude",
                      "Distance to Shore", 
                      "Cyclone Frequency", "Depth", "ClimSST",
                      "SSTA_DHW", "TSA_DHW")
stats1

quantiles1 <- summary1$quantiles
rownames(quantiles1) <- c("Intercept", "Latitude", "Longitude",
                          "Distance to Shore", 
                          "Cyclone Frequency", "Depth", "ClimSST",
                          "SSTA_DHW", "TSA_DHW")
quantiles1

# All of the predictors used above were deemed significant.

### Goodness of Fit Checks for Model 1 ----
### Simple Linear Model with Uninformative Gaussian Priors
# Checking the effective sample sizes.
# All effective sample sizes are very high, well over 10,000!
effectiveSize(samples1)

# R less than 1.1 indicates convergence.
gelman.diag(samples1)

# abs(z) less than 2 indicates convergence.
geweke.diag(samples1[[1]])



## Model 2: Beta regression model (Hanan) following most of Rachel's code but with slight modifications to accommodate for the new model----
## Modeled with Uninformative Gaussian Priors
final_hanan <- final_data3
final_hanan <- na.omit(final_hanan)
Y_2 <- final_hanan$Percent_Bleaching 
Y2 <- Y_2 / 100 # Changing response variable to decimal to fit criteria

X2 <- subset(final_hanan, select = -c(Date, Date_Year, Exposure, Percent_Bleaching))
X2 <- subset(X2, select = -c(Turbidity, SSTA, TSA, Windspeed))
X2 <- as.matrix(X2)
X2 <- scale(X2)

n2 <- length(Y2)
p2 <- ncol(X2)

data2   <- list(Y=Y2,X=X2,n=n2,p=p2)
params1 <- c("alpha", "beta")

burn     <- 10000
n.iter   <- 20000
thin     <- 5

model_string2 <- textConnection("model{
    # Likelihood
    for(i in 1:n){
      Y[i] ~ dbeta(mu[i]*phi, (1-mu[i])*phi) 
      logit(mu[i]) <- alpha + inprod(X[i,], beta[])
    } 
    
    # Priors  
    for(j in 1:p){ 
      beta[j] ~ dnorm(0, 0.01) 
    } 
      
    alpha ~ dnorm(0, 0.01) 
    phi   ~ dgamma(0.1, 0.1) # Shape parameter for beta distribution
      
}")

model2 <- jags.model(model_string2, data=data2, n.chains=2, quiet=TRUE)
#error persists. i've tried increasing burn in. next objective, is to review priors or perhaps try using a different sampler
update(model2, burn, progress.bar="none")
samples2 <- coda.samples(model2, variable.names=params2, n.iter=n.iter, n.thin=thin, progress.bar="none")

summary2 <- summary(samples2)
summary2

