#### ST540 Final Project Analysis
## Programmed by: Tyler Pollard, Hanan Ali, Rachel Hardy
## Date Created: 4 April 2024
## Date Modified: 28 April 2024

# Load Libraries ----
library(data.table)
library(MASS)
library(caret)
library(posterior)
library(bayesplot)
library(rstanarm)
library(rjags)
library(plyr)
library(GGally)
library(tidyverse)
library(tictoc)
library(brms)
library(BayesFactor)
library(tidybayes)
library(bayesanova)


# Clean data for all Models----
## Original Data ----
bleaching_data <- fread("global_bleaching_environmental.csv", 
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
final_data3$Percent_Bleaching_Open <- 
  ifelse(final_data3$Percent_Bleaching == 0, 0.01,
         ifelse(final_data3$Percent_Bleaching == 100, 99.99, 
                final_data3$Percent_Bleaching))
final_data3$Percent_Bleaching_Log <- log(final_data3$Percent_Bleaching_Open)

## Create training/test index vector ----
## for CV and testing model performance and fit
set.seed(52)
trainIndex <- createDataPartition(final_data3$Percent_Bleaching,
                                  p = 0.75,
                                  list = FALSE)
trainIndex <- as.vector(trainIndex)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Model 2: Beta Regression Model =====
## Modeled with Uninformative Gaussian Priors

## Load Data ----
## Y2 is original data with closed support [0,100]
Y2 <- final_data3$Percent_Bleaching
Y2 <- Y2 / 100 # Changing response variable to decimal to fit criteria
mean(Y2)
sd(Y2)

## Use this for Beta regression
Y2open <- final_data3$Percent_Bleaching_Open/100

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
    Date_Month,
    Turbidity,
    Exposure,
    ClimSST,
    SSTA,
    SSTA_DHW,
    Windspeed,
    Depth_m,
    Percent_Bleaching,
    Percent_Bleaching_Open,
    Percent_Bleaching_Log
  ))
X2unscale <- X2
#X2$Exposure <- ifelse(X2$Exposure == "Sheltered", 0, 1)
X2 <- scale(X2)

#### Split Data ----
Y2train <- Y2open[trainIndex]
Y2test <- Y2open[-trainIndex]
X2train <- X2[trainIndex,]
X2test <- X2[-trainIndex,]

## Simulation Variables ----
n2train <- length(Y2train)
n2test <- length(Y2test)
p2 <- ncol(X2train)

## Test models with this these simulation variables:
burn2 <- 2000
n.iters2 <- 4000
thin2 <- 1


## Fit Model ----
model_data <- cbind(
  X2train,
  Percent_Bleaching = Y2train
)
model_data <- data.frame(model_data)

beta_fit <- brm(
  Percent_Bleaching ~ .,
  data = model_data,
  family = Beta(),
  save_pars = save_pars(all = TRUE),
  chains = 2,
  iter = n.iters2,
  warmup = burn2,
  thin = thin2,
  cores = 2,
  seed = 52
)

summary(beta_fit)

beta_fit_draws <- as_draws(beta_fit)

Yppc <- post
Yppd <- posterior_predict(beta_fit, newdata = X2test)

ppc_density_plot2 <- 
  ppc_dens_overlay(Y2train, Yppd) +
  labs(title = "Posterior Predictive Checks of Beta Regression on Training Data",
       subtitle = "Simulated Data Sets Compared to Training Data") +
  theme_bw() +
  legend_none() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.5)),
    plot.subtitle = element_text(size = rel(1)))

YppdL2 <- apply(Yppd, 2, function(x) quantile(x, 0.025))
YppdU2 <- apply(Yppd, 2, function(x) quantile(x, 0.975))
YppdMedian2 <- apply(Yppd, 2, median)
YppdMean2 <- apply(Yppd, 2, mean)

BIAS2B  <- mean(YppdMean2-Y2test)
MSE2B   <- mean((YppdMean2-Y2test)^2)
MAE2B   <- mean(abs(YppdMean2 -Y2test))
MAD2B   <- mean(abs(YppdMedian2-Y2test))
COV2B   <- mean( (YppdL2 <= Y2test) & (Y2test <= YppdU2))
WIDTH2B <- mean(YppdU2-YppdL2)

predStats2B <- c(
  "BIAS" = BIAS2B,
  "MSE" = MSE2B, 
  "MAE" = MAE2B,
  "MAD" = MAD2B,
  "COV" = COV2B,
  "WIDTH" = WIDTH2B
)
predStats2B


# Model 2: Zero Infoated Beta Regression Model =====
## Modeled with Uninformative Gaussian Priors

## Load Data ----
## Y2open is original data with closed support (0,100)
Y2zero <- ifelse(Y2 == 1, 0.9999, Y2)
mean(Y2)
sd(Y2)

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
    Date_Month,
    Turbidity,
    Exposure,
    ClimSST,
    SSTA,
    SSTA_DHW,
    Windspeed,
    Depth_m,
    Percent_Bleaching,
    Percent_Bleaching_Open,
    Percent_Bleaching_Log
  ))
X2unscale <- X2
#X2$Exposure <- ifelse(X2$Exposure == "Sheltered", 0, 1)
X2 <- scale(X2)

#### Split Data ----
Y2train <- Y2zero[trainIndex]
Y2test <- Y2zero[-trainIndex]
X2train <- X2[trainIndex,]
X2test <- X2[-trainIndex,]

## Simulation Variables ----
n2train <- length(Y2train)
n2test <- length(Y2test)
p2 <- ncol(X2train)

## Test models with this these simulation variables:
burn2 <- 2000
n.iters2 <- 4000
thin2 <- 1


## Fit Model ----
model_data <- cbind(
  X2train,
  Percent_Bleaching = Y2train
)
model_data <- data.frame(model_data)

zero_inflated_beta_fit <- brm(
  Percent_Bleaching ~ .,
  data = model_data,
  family = zero_inflated_beta(),
  save_pars = save_pars(all = TRUE),
  chains = 2,
  iter = n.iters2,
  warmup = burn2,
  thin = thin2,
  cores = 2,
  seed = 52
)

summary(zero_beta_fit)

zero_beta_fit_draws <- as_draws(zero_beta_fit)

Yppd <- posterior_predict(zero_beta_fit, newdata = X2test)

ppc_density_plot2 <- 
  ppc_dens_overlay(Y2test, Yppd) +
  labs(title = "Posterior Predictive Checks of Beta Regression on Training Data",
       subtitle = "Simulated Data Sets Compared to Training Data") +
  theme_bw() +
  legend_none() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.5)),
    plot.subtitle = element_text(size = rel(1)))

YppdL2 <- apply(Yppd, 2, function(x) quantile(x, 0.025))
YppdU2 <- apply(Yppd, 2, function(x) quantile(x, 0.975))
YppdMedian2 <- apply(Yppd, 2, median)
YppdMean2 <- apply(Yppd, 2, mean)

BIAS2B  <- mean(YppdMean2-Y2test)
MSE2B   <- mean((YppdMean2-Y2test)^2)
MAE2B   <- mean(abs(YppdMean2 -Y2test))
MAD2B   <- mean(abs(YppdMedian2-Y2test))
COV2B   <- mean( (YppdL2 <= Y2test) & (Y2test <= YppdU2))
WIDTH2B <- mean(YppdU2-YppdL2)

predStats2B <- c(
  "BIAS" = BIAS2B,
  "MSE" = MSE2B, 
  "MAE" = MAE2B,
  "MAD" = MAD2B,
  "COV" = COV2B,
  "WIDTH" = WIDTH2B
)
predStats2B





