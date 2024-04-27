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
library(tictoc)

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
# missingValues <- c()
# for(i in 1:length(colnames(final_data2))){
#   missingValues[i] <- sum(is.na(final_data2[[i]]))
# }
# names(missingValues) <- colnames(final_data2)
# missingValues

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

burn     <- 200
n.iter   <- 400
thin     <- 5

# After running the model below, I removed the following insignificant predictors:
# Turbidity, SSTA, TSA, and Windspeed.
# After removal, the model was run again.

model_string1 <- textConnection("model{

   # Likelihood
    for(i in 1:n){
      Y[i]   ~ dnorm(mu[i],tau) 
      mu[i] <- alpha + inprod(X[i,],beta[])
      
      # For WAIC
      like[i]    <- dnorm(Y[i],mu[i],tau)
    } 
    
    # Priors  
    for(j in 1:p){ 
      beta[j] ~ dnorm(0,0.01) 
    } 
      
    alpha ~  dnorm(0, 0.01) 
    tau   ~  dgamma(0.1, 0.1) 
    
    # Posterior Predicitve Checks
    for(i in 1:n){
      Yppc[i] ~ dnorm(mu[i], tau)
    }
    D1[1] <- min(Yppc[])
    D1[2] <- max(Yppc[])
    D1[3] <- max(Yppc[]) - min(Yppc[])
    D1[4] <- mean(Yppc[])
    D1[5] <- sd(Yppc[])
}")

data1   <- list(Y=Y1,X=X1,n=n1,p=p1)
params1 <- c("alpha", "beta", "D1")

inits <- list(
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 52),
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 8)
)

## Run this all at once to get how long it took
# Start here
tic()
model1 <- jags.model(model_string1, data=data1, inits = inits,
                     n.chains=2, quiet=TRUE)
update(model1, burn, progress.bar="none")
samples1 <- coda.samples(model1, variable.names=params1, n.iter=n.iter, n.thin=thin, progress.bar="none")
toc()
# Stop here
## with burn = 200 and iters = 400 it takes about 1 minute to run

# Save RData in case code aborts and causes termination
# Make sure to reload libraries after termination if it occurs
save(
  final_data3,
  final_rachel,
  model1,
  samples1,
  X1,
  Y1,
  n1,
  p1,
  file = "Model 1 Data.RData"
)

# Use this if R session terminates
# load(file = "Model 1 Data.RData")

summary1 <- summary(samples1)
summary1

#par(mar=c(1,1,1,1)) 
# ^This is because of the error: "figure margins too large"
# Unfortunately it makes the plots kinda wonky but at least we can see them.
plot(samples1)

# This reduces chance of crashing
dev.off()

stats1 <- summary1$statistics[-c(1:5),]
rownames(stats1) <- c("Intercept", "Latitude", "Longitude",
                      "Distance to Shore", 
                      "Cyclone Frequency", "Depth", "ClimSST",
                      "SSTA_DHW", "TSA_DHW")
stats1

quantiles1 <- summary1$quantiles[-c(1:5),]
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

#### Posterior Predictive Checks ----
DPrintnames <- c(
  "Min of Y",
  "Max of Y",
  "Range of Y",
  "Mean of Y",
  "SD of Y"
)

D0 <- c(
  min(Y1),
  max(Y1),
  max(Y1) - min(Y1),
  mean(Y1),
  sd(Y1)
)
names(D0) <- DPrintnames

D1A <- samples1[[1]][,1:5]
colnames(D1A) <- DPrintnames

D1B <- samples1[[2]][,1:5]
colnames(D1B) <- DPrintnames

pval1A <- rep(0, 5)
names(pval1A) <- DPrintnames
pval1B <- rep(0, 5)
names(pval1B) <- DPrintnames

for(j in 1:5){
  plot(density(D1B[,j]), xlab = "D", ylab = "Posterior Probability", 
       xlim = c(min(D1B[,j], D1A[,j], D0[j]), 
                max(D1B[,j], D1A[,j], D0[j])), 
       main = DPrintnames[j])
  lines(density(D1A[,j]), col = "blue")
  abline(v = D0[j], col = "green", lwd = 2)
  legend("topleft", c("D1B", "D1A", "Observed"), 
         col = c("black", "blue", "green"), lwd = 2)
  
  pval1A[j] <- mean(D1A[,j] > D0[j])
  pval1B[j] <- mean(D1B[,j] > D0[j])
}
pval1A
pval1B

dev.off()

### Model Comparison ----
# Compute DIC
dic1   <- dic.samples(model1, n.iter=n.iter, progress.bar="none")

# Compute WAIC
waic1   <- coda.samples(model1, 
                        variable.names=c("like"), 
                        n.iter=n.iter, progress.bar="none")

save(dic1, waic1, file = "Model 1 Comps.RData")

like1   <- waic1[[1]]
fbar1   <- colMeans(like1)
P1      <- sum(apply(log(like1),2,var))
WAIC1   <- -2*sum(log(fbar1))+2*P1

dic1
WAIC1
P1

### Save Model 1 data ----
save(list = setdiff(ls(.GlobalEnv), c("waic1", "like1")),
     file = "Model 1 All Data.Rdata")

# Delete the other previous RData files now to push to GitHub
file.remove("Model 1 Data.Rdata")
file.remove("Model 1 Comps.Rdata")

# Remove data to free up space in environment (we can load it later to compare models)
rm(list=setdiff(ls(), c("bleaching_data", "final_data3")))


## Model 2: Beta regression model (Hanan) ----
## following most of Rachel's code but with slight modifications to accommodate for the new model
## Modeled with Uninformative Gaussian Priors
final_hanan <- final_data3
final_hanan <- na.omit(final_hanan)
Y_2 <- final_hanan$Percent_Bleaching 
Y2 <- Y_2 / 100 # Changing response variable to decimal to fit criteria
mean(Y2)
sd(Y2)
# Change 0's to 0.00001 and 1's to 0.99999 to avoid infinite density
Y2 <- ifelse(Y2 == 0, 0.00001, 
             ifelse(Y2 == 1, 0.99999, Y2))
mean(Y2)
sd(Y2)

X2 <- subset(final_hanan, select = -c(Date, Date_Year, Exposure, Percent_Bleaching))
X2 <- subset(X2, select = -c(Turbidity, SSTA, TSA, Windspeed))
X2 <- as.matrix(X2)
X2 <- scale(X2)

n2 <- length(Y2)
p2 <- ncol(X2)

burn     <- 200
n.iter   <- 400
thin     <- 5

model_string2 <- textConnection("model{
    # Likelihood
    for(i in 1:n){
      Y[i] ~ dbeta(mu[i]*phi, (1-mu[i])*phi) 
      logit(mu[i]) <- alpha + inprod(X[i,], beta[])
      
      # For WAIC
      like[i]    <- dbeta(Y[i], mu[i]*phi, (1-mu[i])*phi)
    } 
    
    # Priors  
    for(j in 1:p){ 
      beta[j] ~ dnorm(0, 0.01) 
    } 
      
    alpha ~ dnorm(0, 0.01) 
    phi   ~ dgamma(0.1, 0.1) # Shape parameter for beta distribution
    
    # Posterior Predicitve Checks
    for(i in 1:n){
      Yppc[i] ~ dbeta(mu[i]*phi, (1-mu[i])*phi) 
    }
    D2[1] <- min(Yppc[])
    D2[2] <- max(Yppc[])
    D2[3] <- max(Yppc[]) - min(Yppc[])
    D2[4] <- mean(Yppc[])
    D2[5] <- sd(Yppc[])
      
}")

data2   <- list(Y=Y2,X=X2,n=n2,p=p2)
params2 <- c("alpha", "beta", "D2")

inits <- list(
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 52),
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 8)
)

## Run this all at once to get how long it took
# Start here
tic()
model2 <- jags.model(model_string2, data=data2, inits = inits,
                     n.chains=2, quiet=TRUE)
#error persists. i've tried increasing burn in. next objective, is to review priors or perhaps 
#try using a different sampler.
# (tyler) did not receive any errors after changing min and max values but model took forever and I
# didn't have the patience to wait 
update(model2, burn, progress.bar="none")
samples2 <- coda.samples(model2, variable.names=params2, n.iter=n.iter, n.thin=thin, progress.bar="none")
toc()
# Stop here

# Save RData in case code aborts and causes termination
save(
  final_data3,
  final_hanan,
  model2,
  samples2,
  X2,
  Y2,
  n2,
  p2,
  file = "Model 2 Data.RData"
)

# Use this if R session terminates
# load(file = "Model 1 Data.RData")

summary2 <- summary(samples2)
summary2

# ^This is because of the error: "figure margins too large"
# Unfortunately it makes the plots kinda wonky but at least we can see them.
plot(samples2)

# This reduces chance of crashing
dev.off()

stats2 <- summary2$statistics[-c(1:5),]
rownames(stats2) <- c("Intercept", "Latitude", "Longitude",
                      "Distance to Shore", 
                      "Cyclone Frequency", "Depth", "ClimSST",
                      "SSTA_DHW", "TSA_DHW")
stats2

quantiles2 <- summary2$quantiles[-c(1:5),]
rownames(quantiles2) <- c("Intercept", "Latitude", "Longitude",
                          "Distance to Shore", 
                          "Cyclone Frequency", "Depth", "ClimSST",
                          "SSTA_DHW", "TSA_DHW")
quantiles2

# All of the predictors used above were deemed significant.

### Goodness of Fit Checks for Model 2 ----
# Checking the effective sample sizes.
effectiveSize(samples2)

# R less than 1.1 indicates convergence.
gelman.diag(samples2)

# abs(z) less than 2 indicates convergence.
geweke.diag(samples2[[1]])

#### Posterior Predictive Checks ----
DPrintnames <- c(
  "Min of Y",
  "Max of Y",
  "Range of Y",
  "Mean of Y",
  "SD of Y"
)

D0B <- c(
  min(Y2),
  max(Y2),
  max(Y2) - min(Y2),
  mean(Y2),
  sd(Y2)
)
names(D0B) <- DPrintnames

D2A <- samples2[[1]][,1:5]
colnames(D2A) <- DPrintnames

D2B <- samples2[[2]][,1:5]
colnames(D2B) <- DPrintnames

pval2A <- rep(0, 5)
names(pval2A) <- DPrintnames
pval2B <- rep(0, 5)
names(pval2B) <- DPrintnames

for(j in 1:5){
  plot(density(D2B[,j]), xlab = "D", ylab = "Posterior Probability", 
       xlim = c(min(D2B[,j], D2A[,j], D0B[j]), 
                max(D2B[,j], D2A[,j], D0B[j])), 
       main = DPrintnames[j])
  lines(density(D2A[,j]), col = "blue")
  abline(v = D0B[j], col = "green", lwd = 2)
  legend("topleft", c("D2B", "D2A", "Observed"), 
         col = c("black", "blue", "green"), lwd = 2)
  
  pval2A[j] <- mean(D2A[,j] > D0B[j])
  pval2B[j] <- mean(D2B[,j] > D0B[j])
}
pval2A
pval2B

### Model Comparison ----
# Compute DIC
dic2   <- dic.samples(model2, n.iter=n.iter, progress.bar="none")

# Compute WAIC
waic2   <- coda.samples(model2, 
                        variable.names=c("like"), 
                        n.iter=n.iter, progress.bar="none")

save(dic2, waic2, file = "Model 2 Comps.RData")

like2   <- waic2[[1]]
fbar2   <- colMeans(like2)
P2      <- sum(apply(log(like2),2,var))
WAIC2   <- -2*sum(log(fbar2))+2*P2

dic2
WAIC2
P2

### Save Model 2 data ----
save(list = setdiff(ls(.GlobalEnv), c("waic2", "like2")),
     file = "Model 2 All Data.Rdata")

# Delete the other previous RData files now to push to GitHub
file.remove("Model 2 Data.Rdata")
file.remove("Model 2 Comps.Rdata")

# Remove data to free up space in environment (we can load it later to compare models)
rm(list=setdiff(ls(), c("bleaching_data", "final_data3")))


## Model 3: Exponential regression model (Tyler) ----
## Modeled with Uninformative Gaussian Priors
final_tyler <- final_data3
final_tyler <- na.omit(final_tyler)
Y3 <- final_tyler$Percent_Bleaching 

X3 <- subset(final_hanan, select = -c(Date, Date_Year, Exposure, Percent_Bleaching))
X3 <- subset(X3, select = -c(Turbidity, SSTA, TSA, Windspeed))
X3 <- as.matrix(X3)
X3 <- scale(X3)

n3 <- length(Y3)
p3 <- ncol(X3)

burn     <- 200
n.iter   <- 400
thin     <- 5

model_string3 <- textConnection("model{
    # Likelihood
    for(i in 1:n){
      Y[i] ~ dbeta(mu[i]*phi, (1-mu[i])*phi) 
      logit(mu[i]) <- alpha + inprod(X[i,], beta[])
      
      # For WAIC
      like[i]    <- dbeta(Y[i], mu[i]*phi, (1-mu[i])*phi)
    } 
    
    # Priors  
    for(j in 1:p){ 
      beta[j] ~ dnorm(0, 0.01) 
    } 
      
    alpha ~ dnorm(0, 0.01) 
    phi   ~ dgamma(0.1, 0.1) # Shape parameter for beta distribution
    
    # Posterior Predicitve Checks
    for(i in 1:n){
      Yppc[i] ~ dbeta(mu[i]*phi, (1-mu[i])*phi) 
    }
    D3[1] <- min(Yppc[])
    D3[2] <- max(Yppc[])
    D3[3] <- max(Yppc[]) - min(Yppc[])
    D3[4] <- mean(Yppc[])
    D3[5] <- sd(Yppc[])
      
}")

data3   <- list(Y=Y3,X=X3,n=n3,p=p3)
params3 <- c("alpha", "beta", "D3")

inits <- list(
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 52),
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 8)
)

## Run this all at once to get how long it took
# Start here
tic()
model3 <- jags.model(model_string3, data=data3, inits = inits,
                     n.chains=2, quiet=TRUE)
#error persists. i've tried increasing burn in. next objective, is to review priors or perhaps try using a different sampler
update(model3, burn, progress.bar="none")
samples3 <- coda.samples(model3, variable.names=params3, n.iter=n.iter, n.thin=thin, progress.bar="none")
toc()
# Stop here

# Save RData in case code aborts and causes termination
save(
  final_data3,
  final_tyler,
  model3,
  samples3,
  X3,
  Y3,
  n3,
  p3,
  file = "Model 3 Data.RData"
)

# Use this if R session terminates
# load(file = "Model 3 Data.RData")

summary3 <- summary(samples3)
summary3

# ^This is because of the error: "figure margins too large"
# Unfortunately it makes the plots kinda wonky but at least we can see them.
plot(samples3)

# This reduces chance of crashing
dev.off()

stats3 <- summary3$statistics[-c(1:5),]
rownames(stats3) <- c("Intercept", "Latitude", "Longitude",
                      "Distance to Shore", 
                      "Cyclone Frequency", "Depth", "ClimSST",
                      "SSTA_DHW", "TSA_DHW")
stats3

quantiles3 <- summary3$quantiles[-c(1:5),]
rownames(quantiles3) <- c("Intercept", "Latitude", "Longitude",
                          "Distance to Shore", 
                          "Cyclone Frequency", "Depth", "ClimSST",
                          "SSTA_DHW", "TSA_DHW")
quantiles3

# All of the predictors used above were deemed significant.

### Goodness of Fit Checks for Model 3 ----
# Checking the effective sample sizes.
effectiveSize(samples3)

# R less than 1.1 indicates convergence.
gelman.diag(samples3)

# abs(z) less than 2 indicates convergence.
geweke.diag(samples3[[1]])

#### Posterior Predictive Checks ----
DPrintnames <- c(
  "Min of Y",
  "Max of Y",
  "Range of Y",
  "Mean of Y",
  "SD of Y"
)

D0C <- c(
  min(Y3),
  max(Y3),
  max(Y3) - min(Y3),
  mean(Y3),
  sd(Y3)
)
names(D0C) <- DPrintnames

D3A <- samples3[[1]][,1:5]
colnames(D3A) <- DPrintnames

D3B <- samples3[[2]][,1:5]
colnames(D3B) <- DPrintnames

pval3A <- rep(0, 5)
names(pval3A) <- DPrintnames
pval3B <- rep(0, 5)
names(pval3B) <- DPrintnames

for(j in 1:5){
  plot(density(D3B[,j]), xlab = "D", ylab = "Posterior Probability", 
       xlim = c(min(D3B[,j], D3A[,j], D0C[j]), 
                max(D3B[,j], D3A[,j], D0C[j])), 
       main = DPrintnames[j])
  lines(density(D3A[,j]), col = "blue")
  abline(v = D0C[j], col = "green", lwd = 2)
  legend("topleft", c("D3B", "D3A", "Observed"), 
         col = c("black", "blue", "green"), lwd = 2)
  
  pval3A[j] <- mean(D3A[,j] > D0C[j])
  pval3B[j] <- mean(D3B[,j] > D0C[j])
}
pval3A
pval3B

### Model Comparison ----
# Compute DIC
dic3   <- dic.samples(model3, n.iter=n.iter, progress.bar="none")

# Compute WAIC
waic3   <- coda.samples(model3, 
                        variable.names=c("like"), 
                        n.iter=n.iter, progress.bar="none")

save(dic3, waic3, file = "Model 3 Comps.RData")

like3   <- waic3[[1]]
fbar3   <- colMeans(like3)
P3      <- sum(apply(log(like3),2,var))
WAIC3   <- -2*sum(log(fbar3))+2*P3

dic3
WAIC3
P3

### Save Model 3 data ----
save(list = setdiff(ls(.GlobalEnv), c("waic3", "like3")),
     file = "Model 3 All Data.Rdata")

# Delete the other previous RData files now to push to GitHub
file.remove("Model 3 Data.Rdata")
file.remove("Model 3 Comps.Rdata")

# Remove data to free up space in environment (we can load it later to compare models)
rm(list=setdiff(ls(), c("bleaching_data", "final_data3")))

