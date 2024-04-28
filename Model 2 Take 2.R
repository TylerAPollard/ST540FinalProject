#### ST540 Final Project Analysis
## Programmed by: Tyler Pollard, Hanan Ali, Rachel Hardy
## Date Created: 4 April 2024
## Date Modified: 

# Load Libraries ----
library(data.table)
library(MASS)
library(bayesplot)
library(rstanarm)
library(rjags)
library(plyr)
library(GGally)
library(tidyverse)
library(caret)
library(tictoc)

# Clean data for all Models----
## Original Data ----
bleaching_data <- fread("Desktop/ST540FinalProject/global_bleaching_environmental.csv", 
                        na.strings = c("", "NA", "nd"))

## Filter Data 1 ----
## Filter to only complete Percent Bleaching and FRRP data set which
## is the 
final_data1 <- bleaching_data |>
  filter(!is.na(Percent_Bleaching)) |>
  filter(Data_Source == "FRRP") |>
  distinct(Site_ID, Sample_ID, .keep_all = TRUE)

## Filter Data 2 ----
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

## Filter Data 3 ----
# Remove rows with missing predictors values
final_data3 <- final_data2[complete.cases(final_data2)]

## Create training/test index vector ----
## for CV and testing model performance and fit
set.seed(52)
trainIndex <- createDataPartition(final_data3$Percent_Bleaching,
                                  p = 0.75,
                                  list = FALSE)
trainIndex <- as.vector(trainIndex)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Model 1: Beta Regression Model =====
## Programmer: Hanan
## Modeled with Uninformative Gaussian Priors

## Load Data ----
### Closed Support Data ----
## Y2 is original data with closed support [0,100]
Y2 <- final_data3$Percent_Bleaching
Y2 <- Y2 / 100 # Changing response variable to decimal to fit criteria
mean(Y2)
sd(Y2)
zeroindex <- Y2 == 0

## We may not need to use these but could be interesting to model data
## separately like in a hierarchical model
Y2zero <- Y2[zeroindex]
Y2nonzero <- Y2[!zeroindex]

### Open Support Data ----
## Change 0's to 0.00001 and 1's to 0.99999 to avoid infinite density
## Y2B is original data to have open support (0,100) 
Y2B <- ifelse(Y2 == 0, 0.001, 
              ifelse(Y2 == 100, 99.999, Y2))
mean(Y2B)
sd(Y2B)
Y2Bzero <- Y2B[zeroindex]
Y2Bnonzero <- Y2B[!zeroindex]

### Log Transformed Data ----
## Y2log is the log of the open supported Y
Y2log <- log(Y2B)
mean(Y2log)
sd(Y2log)
Y2logzero <- Y2log[zeroindex]
Y2lognonzero <- Y2log[!zeroindex]

### Covariate Matrix ----
## Variables removed were from multiple linear regression with no random
## effects. Feel Free to modify this X1 for your specific model if needed
## but rename it X2 for model and so on. Same with Y's just to avoid
## writing over someone 
X2 <- final_data3 |> 
  select(-c(
    Date,
    City_Town_Name,
    Latitude_Degrees,
    Longitude_Degrees,
    Exposure,
    ClimSST,
    SSTA,
    SSTA_DHW,
    Percent_Bleaching
  ))
X2$Exposure <- ifelse(X2$Exposure == "Sheltered", 0, 1)
X2 <- scale(X2)

#### Split Data ----
Y2train <- Y2[trainIndex]
Y2test <- Y2[-trainIndex]
X2train <- X2[trainIndex,]
X2test <- X2[-trainIndex,]

## Simulation Variables ----
n2train <- length(Y2train)
p2 <- ncol(X2train)

## Test models with this these simulation variables:
burn <- 500
n.iter <- 1000
thin <- 5
## We will increase to final model
# burn     <- 5000
# n.iter   <- 10000
# thin     <- 5

## Define Model ----
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


### Compile Model Inputs ----
data2   <- list(Y = Y2train,
                X = X2train,
                n = n2train,
                p = p2)
params2 <- c("alpha", "beta", "D2")

inits <- list(
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 52),
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 8)
)

## Fit Model ----
## Run this all at once to get how long it took
# Start here
tic()
model2 <- jags.model(model_string2, data=data2, inits = inits,
                     n.chains=2, quiet=FALSE)
update(model2, burn, progress.bar="text")
samples2 <- coda.samples(model2, variable.names=params2, n.iter=n.iter, n.thin=thin, 
                         progress.bar="text")
toc()