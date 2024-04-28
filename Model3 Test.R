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
  ifelse(final_data3$Percent_Bleaching == 0, 0.001, 
         ifelse(final_data3$Percent_Bleaching == 100, 99.999, final_data3$Percent_Bleaching))
final_data3$Percent_Bleaching_perc <- (final_data3$Percent_Bleaching_Open)/100
final_data3$Percent_Bleaching_log <- log(final_data3$Percent_Bleaching_Open)


## Create training/test index vector ----
## for CV and testing model performance and fit
set.seed(52)
trainIndex <- createDataPartition(final_data3$Percent_Bleaching,
                                  p = 0.75,
                                  list = FALSE)
trainIndex <- as.vector(trainIndex)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Model 3: Zero-Inflated Regression ----
## Programmed by Tyler
## Modeled Heirarchically

## Load Data ----
### Closed Support Data ----
## Y3 is original data with closed support [0,100]
Y3 <- final_data3$Percent_Bleaching_log
Z3 <- ifelse(final_data3$Percent_Bleaching == 0, 0, 1)

### Covariate Matrix ----
## Variables removed were from multiple linear regression with no random
## effects. Feel Free to modify this X1 for your specific model if needed
## but rename it X2 for model and so on. Same with Y's just to avoid
## writing over someone 
X3 <- final_data3 |> 
  select(-c(
    Date,
    City_Town_Name,
    Latitude_Degrees,
    Longitude_Degrees,
    #Exposure,
    Depth_m,
    ClimSST,
    SSTA,
    SSTA_DHW,
    Percent_Bleaching,
    Percent_Bleaching_Open,
    Percent_Bleaching_log
  ))
X3$Exposure <- ifelse(X3$Exposure == "Sheltered", 0, 1)
X3 <- scale(X3)

#### Split Data ----
Y3train <- Y3[trainIndex]
Y3test <- Y3[-trainIndex]
X3train <- X3[trainIndex,]
X3test <- X3[-trainIndex,]
Z3train <- Z3[trainIndex]
Z3test <- Z3[-trainIndex]

## Simulation Variables ----
n3train <- length(Y3train)
p3 <- ncol(X3train)

## Test models with this these simulation variables:
# burn <- 500
# n.iters <- 3000
# thin <- 5
## We will increase to final model
burn     <- 1000
n.iter   <- 2000
thin     <- 5

## Define Model ----
model_string3 <- textConnection("model{
    # Likelihood
    for(i in 1:n){
      Y[i] ~ dnorm(mu[i], taue)
      mu[i] <- exp(alphaY + inprod(X[i,], betaY[]))
      
      #Z[i] ~ dbern(q[i])
      #logit(q[i]) <- alphaZ + inprod(X[i,], betaZ[])
    } 
    
    # Priors  
    for(j in 1:p){ 
      betaY[j] ~ dnorm(0, 0.01)
      #betaZ[j] ~ dnorm(0, 0.01)
    } 
    alphaY ~ dnorm(0, 0.01)
    #alphaZ ~ dnorm(0, 0.01)
    
    taue ~ dgamma(0.1,0.1)
    sigma <- 1/sqrt(taue)
    
    # Posterior Predicitve Checks
    for(i in 1:n){
      Yppc[i] ~ dnorm(mu[i], taue)
      #Yppc[i] <- exp(YppcLog[i])
    }
    D3[1] <- min(Yppc[])
    D3[2] <- max(Yppc[])
    D3[3] <- max(Yppc[]) - min(Yppc[])
    D3[4] <- mean(Yppc[])
    D3[5] <- sd(Yppc[])
    D3[6] <- 
}")

### Compile Model Inputs ----
data3   <- list(Y = Y3train,
                #Z = Z3train,
                X = X3train,
                n = n3train,
                p = p3)
params3 <- c("alphaY", "betaY", 
             #"alphaZ", "betaZ",
             "sigma",
             "D3"
             )

inits <- list(
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 52),
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 8)
)

## Fit Model ----
## Run this all at once to get how long it took
# Start here
tic()
model3 <- jags.model(model_string3, data=data3, inits = inits,
                     n.chains=2, quiet=FALSE)
update(model3, burn, progress.bar="text")
samples3 <- coda.samples(model3, variable.names=params3, n.iter=n.iter, 
                         n.thin=thin, progress.bar="text")
toc()
# Stop here

### SAVE CHECKPOINT ----
# Save RData in case code aborts and causes termination
# Make sure to reload libraries after termination if it occurs
filename3 <- paste0(
  "Model 3 Fit Data ", 
  "(B", burn, "-",
  "I", n.iter, "-",
  "T", thin, 
  ").RData"
)
save(
  list = setdiff(
    ls(.GlobalEnv), 
    c("bleaching_data", "final_data3", "final_data2")
  ),
  file = filename3
)
# Use this if R session terminates
# load(file = "Model 3 Fit Data.RData")

### Convergence Diagnostics ----
#### Trace Plots ----
# Trace and Density plots
# color_scheme_set("mix-blue-red")
# mcmc_trace(samples3)
# mcmc_rank_overlay(samples3)
# dimnames(samples3)
plot(samples3)

# Remove plots if it is bogging down environment
dev.off()

#### Effective Sample Size ----
# Checking the effective sample sizes.
# All effective sample sizes are very high, well over 30,000!
effectiveSize(samples3)

#### Gelman-Rubin Diagnostics ----
# R less than 3.3 indicates convergence.
gelman.diag(samples3)

#### Geweke Diagnostics ----
# abs(z) less than 2 indicates convergence.
geweke.diag(samples3)

### Model Summaries ----
summary3 <- summary(samples3)

#### Parameter Estimates ----
paramNames3A <- c("Intercept", colnames(X3))
paramNames3B <- rep(paramNames3A, each = 2)
paramNames3C <- paste0(c("Y_", "Z_"), paramNames3B)
paramNames3D <- c(paramNames3C, "sigma")

stats3 <- summary3$statistics[-c(1:5),]
rownames(stats3) <- c(paramNames, "sigma")
stats3

#### Parameter 95% CIs ----
# CI = Credible Interval
quantiles3 <- summary3$quantiles[-c(1:5),]
rownames(quantiles3) <- c(paramNames, "sigma")
quantiles3

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
Y3Exp <- final_data3$Percent_Bleaching_Open
Y3trainExp <- Y3Exp[trainIndex]
Y3testExp <- Y3Exp[-trainIndex]

D0 <- c(
  min(Y3train),
  max(Y3train),
  max(Y3train) - min(Y3train),
  mean(Y3train),
  sd(Y3train)
)
D0 <- c(
  min(Y3trainExp),
  max(Y3trainExp),
  max(Y3trainExp) - min(Y3trainExp),
  mean(Y3trainExp),
  sd(Y3trainExp)
)
names(D0) <- DPrintnames

# Chain 3 PPCs
D3A <- samples3[[1]][,1:5]
colnames(D3A) <- DPrintnames

# Chain 2 PPCs
D3B <- samples3[[2]][,1:5]
colnames(D3B) <- DPrintnames

# Create empty vectors to store Bayesian p-values for each chain
pval3A <- rep(0, 5)
names(pval3A) <- DPrintnames
pval3B <- rep(0, 5)
names(pval3B) <- DPrintnames

#### Plot all PPCs ----
#### and calculate p-values
# May end up converting this to ggplot for report

# For all plots in one
par(mfrow = c(3,2))

# For individual plots
# dev.off()

for(j in 1:5){
  plot(density(D3B[,j]), xlab = "D", ylab = "Posterior Probability", 
       xlim = c(min(D3B[,j], D3A[,j], D0[j]), 
                max(D3B[,j], D3A[,j], D0[j])), 
       main = DPrintnames[j])
  lines(density(D3A[,j]), col = "blue")
  abline(v = D0[j], col = "green", lwd = 2)
  legend("topleft", c("D3B", "D3A", "Observed"), 
         col = c("black", "blue", "green"), lwd = 2)
  
  pval3A[j] <- mean(D3A[,j] > D0[j])
  pval3B[j] <- mean(D3B[,j] > D0[j])
}
pval3A
pval3B

# Remove plots if it is bogging down environment
dev.off()

## Model Comparison ----
### DIC ----
dic3   <- dic.samples(model3, n.iter=n.iter, progress.bar="text")

#### SAVE CHECKPOINT ----
save(dic3, file = "Model 3 DIC.RData")

### WAIC ----
waic3   <- coda.samples(model3, 
                        variable.names=c("like"), 
                        n.iter=n.iter, progress.bar="text")

#### SAVE CHECKPOINT ----
save(waic3, file = "Model 3 WAIC.RData")

#### Manually Compute WAIC and P ----
like3   <- waic3[[3]]
fbar3   <- colMeans(like3)
P3      <- sum(apply(log(like3),2,var))
WAIC3   <- -2*sum(log(fbar3))+2*P3

### Output Results ----
dic3
WAIC3
P3

## FINAL SAVE ----
save(list = setdiff(ls(.GlobalEnv), c("waic3", "like3")),
     file = filename3)

# Delete the other previous RData files now to push to GitHub
file.remove("Model 3 DIC.Rdata")
file.remove("Model 3 WAIC.Rdata")

# Remove data to free up space in environment (we can load it later to compare models)
rm(list=setdiff(ls(), c("bleaching_data", "final_data3")))



