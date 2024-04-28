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

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Model 1: Beta Regression Model =====
## Programmer: Hanan
## Modeled with Uninformative Gaussian Priors

## Load Data ----
### Closed Support Data ----
## Y2 is original data with closed support [0,100]
Y2 <- final_data3$Percent_Bleaching_Open
Y2 <- Y2 / 100 # Changing response variable to decimal to fit criteria
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
    Percent_Bleaching,
    Percent_Bleaching_Open,
    Percent_Bleaching_Log
  ))
#X2$Exposure <- ifelse(X2$Exposure == "Sheltered", 0, 1)
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
n.iters <- 1000
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
samples2 <- coda.samples(model2, variable.names=params2, n.iter=n.iters,
                         n.thin=thin, progress.bar="text")
toc()
# Stop here

### SAVE CHECKPOINT ----
# Save RData in case code aborts and causes termination
# Make sure to reload libraries after termination if it occurs
filename2 <- paste0(
  "Model 2 Fit Data ", 
  "(B", burn, "-",
  "I", n.iters, "-",
  "T", thin, 
  ").RData"
)
save(
  list = setdiff(
    ls(.GlobalEnv), 
    c("bleaching_data", "final_data1", "final_data2")
  ),
  file = filename2
)
# Use this if R session terminates
# load(file = "Model 2 Fit Data.RData")

### Convergence Diagnostics ----
#### Trace Plots ----
# Trace and Density plots
# color_scheme_set("mix-blue-red")
# mcmc_trace(samples1)
# mcmc_rank_overlay(samples1)
# dimnames(samples1)
plot(samples2)

# Remove plots if it is bogging down environment
dev.off()

#### Effective Sample Size ----
# Checking the effective sample sizes.
# All effective sample sizes are very high, well over 10,000!
effectiveSize(samples2)

#### Gelman-Rubin Diagnostics ----
# R less than 1.1 indicates convergence.
gelman.diag(samples2)

#### Geweke Diagnostics ----
# abs(z) less than 2 indicates convergence.
geweke.diag(samples2)

### Model Summaries ----
summary2 <- summary(samples2)

#### Parameter Estimates ----
stats2 <- summary2$statistics[-c(1:5),]
rownames(stats2) <- c("Intercept", colnames(X2))
stats2

#### Parameter 95% CIs ----
# CI = Credible Interval
quantiles2 <- summary2$quantiles[-c(1:5),]
rownames(quantiles2) <- c("Intercept", colnames(X2))
quantiles2

## Goodness of Fit Checks ----
### Posterior Predictive Checks ----
# Create naming vector to be used throughout checks and plots
DPrintnames <- c(
  "Min of Y",
  "Max of Y",
  "Range of Y",
  "Mean of Y",
  "SD of Y"
)

# Calculate observed values for checks
D0 <- c(
  min(Y2train),
  max(Y2train),
  max(Y2train) - min(Y2train),
  mean(Y2train),
  sd(Y2train)
)
names(D0) <- DPrintnames

# Chain 1 PPCs
D2A <- samples2[[1]][,1:5]
colnames(D2A) <- DPrintnames

# Chain 2 PPCs
D2B <- samples2[[2]][,1:5]
colnames(D2B) <- DPrintnames

# Create empty vectors to store Bayesian p-values for each chain
pval2A <- rep(0, 5)
names(pval2A) <- DPrintnames
pval2B <- rep(0, 5)
names(pval2B) <- DPrintnames

#### Plot all PPCs ----
#### and calculate p-values
# May end up converting this to ggplot for report

# For all plots in one
par(mfrow = c(3,2))

# For individual plots
# dev.off()

for(j in 1:5){
  plot(density(D2B[,j]), xlab = "D", ylab = "Posterior Probability", 
       xlim = c(min(D2B[,j], D2A[,j], D0[j]), 
                max(D2B[,j], D2A[,j], D0[j])), 
       main = DPrintnames[j])
  lines(density(D2A[,j]), col = "blue")
  abline(v = D0[j], col = "green", lwd = 2)
  legend("topleft", c("D1B", "D2A", "Observed"), 
         col = c("black", "blue", "green"), lwd = 2)
  
  pval2A[j] <- mean(D2A[,j] > D0[j])
  pval2B[j] <- mean(D2B[,j] > D0[j])
}
pval2A
pval2B

# Remove plots if it is bogging down environment
dev.off()

## Model Comparison ----
### DIC ----
dic2   <- dic.samples(model2, n.iter=n.iters, progress.bar="text")

#### SAVE CHECKPOINT ----
save(dic2, file = "Model 2 DIC.RData")

### WAIC ----
waic2   <- coda.samples(model2, 
                        variable.names=c("like"), 
                        n.iter=n.iters, progress.bar="text")

#### SAVE CHECKPOINT ----
save(waic2, file = "Model 2 WAIC.RData")

#### Manually Compute WAIC and P ----
like2   <- waic2[[1]]
fbar2   <- colMeans(like2)
P2      <- sum(apply(log(like2),2,var))
WAIC2   <- -2*sum(log(fbar2))+2*P2

### Output Results ----
dic2
WAIC2
P2

## FINAL SAVE ----
save(list = setdiff(
  ls(.GlobalEnv), 
  c("waic2", "like2", "bleaching_data", "final_data1", "final_data2")),
  file = filename2
)

# Delete the other previous RData files now to push to GitHub
file.remove("Model 2 DIC.RData")
file.remove("Model 2 WAIC.RData")

# Remove data to free up space in environment (we can load it later to compare models)
rm(list=setdiff(ls(), c("bleaching_data", "final_data3")))

