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
# Model 1: Simple Linear Model =====
## Modeled with Uninformative Gaussian Priors

## Load Data ----
## Y1 is original data with closed support [0,100]
Y1 <- final_data3$Percent_Bleaching

### Covariate Matrix ----
## Variables removed were from multiple linear regression with no random
## effects. Variables that were removed below were removed for 
## convenience or because they were insignificant (like Windspeed).
X1 <- final_data3 |> 
  select(-c(
    Date,
    City_Town_Name,
    Latitude_Degrees,
    Longitude_Degrees,
    Exposure,
    ClimSST,
    SSTA,
    SSTA_DHW,
    Windspeed,
    Percent_Bleaching,
    Percent_Bleaching_Open,
    Percent_Bleaching_Log
  ))
# X1$Exposure <- ifelse(X1$Exposure == "Sheltered", 0, 1)
X1 <- scale(X1)

#### Split Data ----
# Observed data
Y1train <- Y1[trainIndex]
X1train <- X1[trainIndex,]

# Set aside for prediction
Y1test <- Y1[-trainIndex]
X1test <- X1[-trainIndex,]

## Simulation Variables ----
n1 <- length(Y1)
n1train <- length(Y1train)
n1test <- length(Y1test)
p1 <- ncol(X1train)

## Test models with this these simulation variables:
# burn <- 500
# n.iter <- 1000
# thin <- 5
## We will increase to final model
burn1     <- 5000
n.iters1   <- 10000
thin1     <- 5

## Define Model ----
model_string1 <- textConnection("model{
  # Likelihood
  for(i in 1:n){
    Y[i] ~ dnorm(mu[i],tau)
    mu[i] <- alpha + inprod(X[i,],beta[])
    
    # For WAIC
    # like[i] <- dnorm(Y[i],mu[i],tau)
  } 
    
  # Priors  
  for(j in 1:p){ 
    beta[j] ~ dnorm(0,0.01) 
  } 
    
  alpha ~  dnorm(0, 0.01) 
  tau   ~  dgamma(0.1, 0.1) 
  sigma <- 1/sqrt(tau)
    
  # Posterior Predicitve Checks
  for(i in 1:n){
    Yppc[i] ~ dnorm(mu[i], tau)
  }
  D1[1] <- min(Yppc[])
  D1[2] <- max(Yppc[])
  D1[3] <- max(Yppc[]) - min(Yppc[])
  D1[4] <- mean(Yppc[])
  D1[5] <- sd(Yppc[])
  
  # Predictions
  for(i in 1:n_pred){
    Ypred[i] ~ dnorm(mu_pred[i], tau)
    mu_pred[i] <- alpha + inprod(Xpred[i,], beta[])
  }
}")

### Compile Model Inputs ----
data1   <- list(Y = Y1train,
                X = X1train,
                n = n1train,
                p = p1,
                Xpred = X1test,
                n_pred = n1test
)
params1 <- c("alpha", "beta", "sigma", "D1", "Yppc", "Ypred")

inits <- list(
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 52),
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 8)
)

## Fit Model ----
## Run this all at once to get how long it took
# Start here
tic()
model1 <- jags.model(model_string1, data=data1, inits = inits,
                     n.chains=2, quiet=FALSE)
update(model1, burn1, progress.bar="text")
samples1 <- coda.samples(model1, variable.names=params1, n.iter=n.iters1, n.thin=thin1, 
                         progress.bar="text")
toc()
# Stop here

### SAVE CHECKPOINT ----
# Save RData in case code aborts and causes termination
# Make sure to reload libraries after termination if it occurs
filename1 <- paste0(
  "Model 1 Fit Data ", 
  "(B", burn1, "-",
  "I", n.iters1, "-",
  "T", thin1, 
  ").RData"
)
save(
  list = setdiff(
    ls(.GlobalEnv), 
    c("bleaching_data")
  ),
  file = filename1
)

# Use this if R session terminates
# load(file = "Model 1 Fit Data.RData")

### Save samples by variables ----
samples1A <- samples1[[1]]
samples1B <- samples1[[2]]
colnames(samples1A)

#### Parameters ----
paramSamps1A <- samples1A[,c(602:610)]
paramSamps1B <- samples1B[,c(602:610)]

### Convergence Diagnostics ----
#### Trace Plots ----


# Remove plots if it is bogging down environment
dev.off()

#### Effective Sample Size ----
# Checking the effective sample sizes.
# All effective sample sizes are very high, well over 10,000!
effectiveSize(samples1)

#### Gelman-Rubin Diagnostics ----
# R less than 1.1 indicates convergence.
gelman.diag(samples1)

#### Geweke Diagnostics ----
# abs(z) less than 2 indicates convergence.
geweke.diag(samples1)

### Model Summaries ----
summary1 <- summary(samples1)

#### Parameter Estimates ----
stats1 <- summary1$statistics[-c(1:5),]
rownames(stats1) <- c("Intercept", colnames(X1))
stats1

#### Parameter 95% CIs ----
# CI = Credible Interval
quantiles1 <- summary1$quantiles[-c(1:5),]
rownames(quantiles1) <- c("Intercept", colnames(X1))
quantiles1

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
  min(Y1train),
  max(Y1train),
  max(Y1train) - min(Y1train),
  mean(Y1train),
  sd(Y1train)
)
names(D0) <- DPrintnames

# Chain 1 PPCs
D1A <- samples1[[1]][,1:5]
colnames(D1A) <- DPrintnames

# Chain 2 PPCs
D1B <- samples1[[2]][,1:5]
colnames(D1B) <- DPrintnames

# Create empty vectors to store Bayesian p-values for each chain
pval1A <- rep(0, 5)
names(pval1A) <- DPrintnames
pval1B <- rep(0, 5)
names(pval1B) <- DPrintnames

#### Plot all PPCs ----
#### and calculate p-values
# May end up converting this to ggplot for report

# For all plots in one
par(mfrow = c(3,2))

# For individual plots
# dev.off()

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

# Remove plots if it is bogging down environment
dev.off()

## Model Comparison ----
### DIC ----
dic1   <- dic.samples(model1, n.iter=n.iter, progress.bar="text")

#### SAVE CHECKPOINT ----
save(dic1, file = "Model 1 DIC.RData")

### WAIC ----
waic1   <- coda.samples(model1, 
                        variable.names=c("like"), 
                        n.iter=n.iter, progress.bar="text")

#### SAVE CHECKPOINT ----
save(waic1, file = "Model 1 WAIC.RData")

#### Manually Compute WAIC and P ----
like1   <- waic1[[1]]
fbar1   <- colMeans(like1)
P1      <- sum(apply(log(like1),2,var))
WAIC1   <- -2*sum(log(fbar1))+2*P1

### Output Results ----
dic1
WAIC1
P1

## Cross Validation ----
# Observed data
Y1o <- Y1[trainIndex] # Observed (Y1train)
X1o <- X1[trainIndex,] # Observed (X1train)

# Set aside for prediction
Y1p <- Y1[-trainIndex] # For prediction (Y1test)
X1p <- X1[-trainIndex,] # For prediction (X1test)

n1o    <- length(Y1o)
n1p    <- length(Y1p)
p1     <- ncol(X1o)

# Model below
model_string <- "model{

  # Likelihood
  for(i in 1:n1o){
    Y1o[i]   ~ dnorm(muo[i],tau)
    muo[i] <- alpha + inprod(X1o[i,],beta[])
  }

  # Prediction
  for(i in 1:n1p){
    Y1p[i]  ~ dnorm(mup[i],tau)
    mup[i] <- alpha + inprod(X1p[i,],beta[])
  }

  # Priors
  for(j in 1:p1){
    beta[j] ~ dnorm(0,0.01)
  }
  
  alpha  ~ dnorm(0, 0.01)
  tau    ~  dgamma(0.1, 0.1) 
  sigma  <- 1/sqrt(tau)
  
}"

# NOTE: Y1p is not sent to JAGS!
model1_cv <- jags.model(textConnection(model_string), 
                        data = list(Y1o=Y1o,n1o=n1o,n1p=n1p,p1=p1,X1o=X1o,X1p=X1p))
update(model1_cv, burn1, progress.bar="none")
samp1_cv <- coda.samples(model1_cv, 
                         variable.names=c("alpha", "beta", "Y1p"), 
                         n.iter=n.iters1, progress.bar="text")

#### SAVE CHECKPOINT ----
save(dic1, file = "Model 1 DIC.RData")

# Summary
summary1_cv <- summary(samp1_cv[,-c(1:n1p)])
summary1_cv

# Extract the samples for each parameter
samps1       <- samp1_cv[[1]]
Yp.samps1    <- samps1[,1:n1p] 
alpha.samps1 <- samps1[,n1p+1]
beta.samps1  <- samps1[,n1p+1+1:p1]
sigma.samps1  <- samps1[,ncol(samps1)]

# Beta means
beta.mn1  <- colMeans(beta.samps1)
beta.mn1

# Sigma mean
sigma.mn1 <- mean(sigma.samps1)
sigma.mn1

# Alpha (intercept) mean
alpha.mn1 <- mean(alpha.samps1)
alpha.mn1

# Graphical representation of plug-in vs PPD vs truth
# Only graphing the first few!
for(j in 1:10){
  
  # Plug-in
  mu <- alpha.mn1+sum(X1p[j,]*beta.mn1)
  y  <- rnorm(n.iters1,mu,sigma.mn1)
  plot(density(y),col=2,xlab="Y",main="PPD")
  
  # PPD
  lines(density(Yp.samps1[,j]))
  
  # Truth
  abline(v=Y1p[j],col=3,lwd=2)
  
  legend("topright",c("PPD","Plug-in","Truth"),col=1:3,lty=1,inset=0.05)
}

# Remove plots if it is bogging down environment
dev.off()

# Create empty data frames to store values
df_mu1   <- data.frame()
df_low1  <- data.frame()
df_high1 <- data.frame()

# Looping through values to calculate the mean and credible interval
for(j in 1:n1p){
  mu   <- alpha.mn1 + sum(X1p[j,]*beta.mn1)
  low  <- mu - 1.96*sigma.mn1 
  high <- mu + 1.96*sigma.mn1
  
  df_mu1   <- rbind(df_mu1, mu)
  df_low1  <- rbind(df_low1, low)
  df_high1 <- rbind(df_high1, high)
}

# Predictions baby!!!!!
df_cv1_pred <- cbind(df_mu1, df_low1, df_high1)

# Creating vectors of the predictions and actual values
df_mu1 <- as.matrix(df_mu1)
Y1p <- as.matrix(Y1p)

# Create empty data frame to store values
est1 <- c()

# Looping through and taking the abs difference of the two vectors
for(i in 1:n1p){
  diff <- abs(Y1p[i]-df_mu1[i])
  est1[i]  <- diff
}

# First column of the est1 data frame has wonky name
mean_abs_error1 <- mean(est1)
mean_abs_error1

## FINAL SAVE ----
save(list = setdiff(
  ls(.GlobalEnv), 
  c("waic1", "like1", "bleaching_data", "final_data1", "final_data2")
),
file = filename1
)

# Delete the other previous RData files now to push to GitHub
file.remove("Model 1 DIC.Rdata")
file.remove("Model 1 WAIC.Rdata")

# Remove data to free up space in environment (we can load it later to compare models)
# rm(list=setdiff(ls(), c("bleaching_data", "final_data3")))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Model 2: Beta Regression Model =====
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
    Depth_m,
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
n2test <- length(Y2test)
p2 <- ncol(X2train)

## Test models with this these simulation variables:
burn2 <- 2000
n.iters2 <- 4000
thin2 <- 5
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
      # like[i]    <- dbeta(Y[i], mu[i]*phi, (1-mu[i])*phi)
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
    
    # Predictions
    for(i in 1:n_pred){
      Ypred[i] ~ dbeta(mu_pred[i]*phi, (1-mu_pred[i])*phi) 
      logit(mu_pred[i]) <- alpha + inprod(Xpred[i,], beta[])
    }
      
}")

### Compile Model Inputs ----
data2   <- list(Y = Y2train,
                X = X2train,
                n = n2train,
                p = p2,
                Xpred = X2test,
                n_pred = n2test
)
params2 <- c("alpha", "beta", "phi", "D2", "Yppc", "Ypred")

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
update(model2, burn2, progress.bar="text")
samples2 <- coda.samples(model2, variable.names=params2, n.iter=n.iters2,
                         n.thin=thin2, progress.bar="text")
toc()
# Stop here

### SAVE CHECKPOINT ----
# Save RData in case code aborts and causes termination
# Make sure to reload libraries after termination if it occurs
filename2 <- paste0(
  "Model 2 Fit Data ", 
  "(B", burn2, "-",
  "I", n.iters2, "-",
  "T", thin2, 
  ").RData"
)
save(
  list = setdiff(
    ls(.GlobalEnv), 
    c("bleaching_data")
  ),
  file = filename2
)
# Use this if R session terminates
# load(file = "Model 2 Fit Data.RData")

### Parse Samples by Variable ----
samples2A <- samples2[[1]]
samples2B <- samples2[[2]]

#### Parameters ----
paramSamps2A <- samples2A[ ,2400:2406]
paramSamps2B <- samples2B[ ,2400:2406]
colnames(paramSamps2A) <- c("Intercept", colnames(X2), "phi")
colnames(paramSamps2B) <- c("Intercept", colnames(X2), "phi")
paramSamps2 <- mcmc.list(paramSamps2A,paramSamps2B)
paramSamps2C <- as_draws_array(paramSamps2)

#### PPCs ----
##### Data ----
YppcSamps2A <- samples2A[ ,6:1803]
YppcSamps2B <- samples2B[ ,6:1803]
YppcSamps2 <- mcmc.list(YppcSamps2A,YppcSamps2B)

##### Checks ----
ppcSamps2A <- samples2A[ ,1:5]
ppcSamps2B <- samples2B[ ,1:5]
DPrintnames <- c(
  "Min of Y",
  "Max of Y",
  "Range of Y",
  "Mean of Y",
  "SD of Y"
)
colnames(ppcSamps2A) <- DPrintnames
colnames(ppcSamps2B) <- DPrintnames
ppcSamps2 <- mcmc.list(ppcSamps2A,ppcSamps2B)

#### PPDs ----
ppdSamps2A <- samples2A[ ,1804:2399]
ppdSamps2B <- samples2B[ ,1804:2399]
ppdSamps2 <- mcmc.list(ppdSamps2A, ppdSamps2B)
PPDcomb2 <- rbind(ppdSamps2A,ppdSamps2B)



### Convergence Diagnostics ----
#### Trace Plots ----
# Trace and Density plots
color_scheme_set("teal")
#color_scheme_set("brewer-Dark2")
mcmc_combo(paramSamps2,
           combo = c("trace", "dens_overlay"),
           # pars = c("alpha", 
           #          paste0("beta[", 1:p2, "]"),
           #          "phi"),
           widths = c(2,1))

# mcmc_rank_overlay(samples1)
# dimnames(samples1)
#plot(paramSamps2A)

# Remove plots if it is bogging down environment
dev.off()

#### Effective Sample Size ----
# Checking the effective sample sizes.
# All effective sample sizes are very high, well over 10,000!
effectiveSize(paramSamps2)

#### Gelman-Rubin Diagnostics ----
# R less than 1.1 indicates convergence.
gelman.diag(paramSamps2)

#### Geweke Diagnostics ----
# abs(z) less than 2 indicates convergence.
geweke.diag(paramSamps2)

### Model Summaries ----
summary2 <- summary(paramSamps2)

#### Parameter Estimates ----
stats2 <- summary2$statistics
stats2

#### Parameter 95% CIs ----
# CI = Credible Interval
quantiles2 <- summary2$quantiles
quantiles2

#### Plot parameter distributions
mcmc_areas(
  paramSamps2,
  pars = colnames(X2),
  point_est = "mean",
  prob = 0.95) +
  labs(
    title = "Posterior Distributions of Predictors",
    subtitle = "95% Credible Interval about the Point Estimate"
  ) +
  vline_0() +
  theme_bw()

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

# For all plots in one
par(mfrow = c(3,2))

# For individual plots
# dev.off()

for(j in 1:5){
  plot(density(D2B[,j]), xlab = "D", ylab = "Posterior Probability", 
       xlim = c(min(D2B[,j], D2A[,j], D0[j]), 
                max(D2B[,j], D2A[,j], D0[j])), 
       main = DPrintnames[j])
  lines(density(D2B[,j]), col = "blue")
  abline(v = D0[j], col = "green", lwd = 2)
  legend("topleft", c("D1B", "D2A", "Observed"), 
         col = c("black", "blue", "green"), lwd = 2)
  
  pval2A[j] <- mean(D2B[,j] > D0[j])
  pval2B[j] <- mean(D2B[,j] > D0[j])
}
pval2A
pval2B

# Remove plots if it is bogging down environment
dev.off()

#### Plot all PPCs ----
#### and calculate p-values
# May end up converting this to ggplot for report
ppc_obsN2 <- 30
ppc_obs2 <- sample(1:length(Y2train), ppc_obsN2)
ppc_density_plot2 <- ppc_dens_overlay(
  Y2train, YppcSamps2A[ppc_obs2,]) +
  labs(title = "Posterior Predictive Checks vs Training Set",
       subtitle = "N = 30 Randomly selected posteriors from JAGS") +
  theme_bw() +
  legend_none()

ppc_q2.5_plot2 <- ppc_stat(
  Y2train, YppcSamps2A, 
  stat = function(y) quantile(y, 0.025)) +
  labs(title = "2.5% Quantile") +
  theme_bw() +
  legend_none()
ppc_q97.5_plot2 <- ppc_stat(
  Y2train, YppcSamps2A, 
  stat = function(y) quantile(y, 0.975)) +
  labs(title = "97.5% Quantile") +
  theme_bw() +
  legend_none()
ppc_median_plot2 <- ppc_stat(
  Y2train, YppcSamps2A, 
  stat = "median") +
  labs(title = "Median") +
  theme_bw() +
  legend_none()
ppc_mad_plot2 <- ppc_stat(
  Y2train, YppcSamps2A, 
  stat = "mad") +
  labs(title = "MAD") +
  theme_bw() +
  legend_none()
ppc_mean_plot2 <- ppc_stat(
  Y2train, YppcSamps2A, 
  stat = "mean") +
  labs(title = "Mean") +
  theme_bw() +
  legend_none()
ppc_sd_plot2 <- ppc_stat(
  Y2train, YppcSamps2A, 
  stat = "sd") +
  labs(title = "Standard Deviation") +
  theme_bw() +
  legend_none()

ppc_lay2 <- rbind(c(1,1),
                  c(2,3),
                  c(4,5),
                  c(6,7)
)
bayesplot_grid(
  plots = list(
    ppc_density_plot2,
    ppc_q2.5_plot2,
    ppc_q97.5_plot2,
    ppc_median_plot2,
    ppc_mad_plot2,
    ppc_mean_plot2,
    ppc_sd_plot2),
  grid_args = list(
    layout_matrix = ppc_lay2
  )
)

# Remove saved plots for saving
rm(ppc_density_plot2, ppc_q2.5_plot2, ppc_median_plot2, ppc_q97.5_plot2,
   ppc_mean_plot2, ppc_sd_plot2)

#### Plot PDDs ----
PPDmean2 <- apply(PPDcomb2, 2, mean)
PPDmedian2 <- apply(PPDcomb2, 2, median)
PPDlb2 <- apply(PPDcomb2, 2, function(x) quantile(x, 0.025))
PPDub2 <- apply(PPDcomb2, 2, function(x) quantile(x, 0.975))

PPDcomb2B <- mcmc(PPDcomb2)

deviance2A <- abs(PPDmedian2 - Y2test)
MAD2B <- mean(deviance2A)

yt <- example_y_data()
yrept <- example_yrep_draws()
yreptMed <- apply(yrept, 2, median)
mean(abs(yreptMed - yt))
mean(abs(median(yrept)))

ppc_stat(yt, yrept, stat = "mad") + theme_bw()

ppdComb_density_plot2 <- ppc_dens_overlay(
  Y2test, PPDcomb2) +
  labs(title = "Posterior Predictive Distribution vs Test Set",
       subtitle = "N = 30 Randomly selected posteriors from JAGS") +
  theme_bw() +
  legend_none()
ppdComb_density_plot2


ppdComb_density_plot2B <- ppc_dens_overlay(
  Y2test, PPDcomb2) +
  labs(title = "Posterior Predictive Distribution vs Test Set",
       subtitle = "95% Credible Interval about Mean") +
  theme_bw() +
  legend_none()
ppdComb_density_plot2B

ppdComb_density_plot2B <- ppc_dens_overlay_data(
  Y2test, PPDcomb2)

ppd_ribbon(PPDcomb2, prob = 0.95)

ppd_obsN2 <- 30
ppd_obs2 <- sample(1:length(Y2test), ppc_obsN2)
ppd_density_plot2 <- ppc_dens_overlay(
  Y2test, PPDcomb2[ppc_obs2,]) +
  labs(title = "Posterior Predictive Distribution vs Test Set",
       subtitle = "N = 30 Randomly selected posteriors from JAGS") +
  theme_bw() +
  legend_none()

ppd_q2.5_plot2 <- ppc_stat(
  Y2test, PPDcomb2, 
  stat = function(y) quantile(y, 0.025)) +
  labs(title = "2.5% Quantile") +
  theme_bw() +
  legend_none()

ppd_q97.5_plot2 <- ppc_stat(
  Y2test, PPDcomb2, 
  stat = function(y) quantile(y, 0.975)) +
  labs(title = "97.5% Quantile") +
  theme_bw() +
  legend_none()

ppd_median_plot2 <- ppc_stat(
  Y2test, PPDcomb2, 
  stat = "median") +
  labs(title = "Median") +
  theme_bw() +
  legend_none()

ppd_mad_plot2 <- ppc_stat(
  Y2test, PPDcomb2, 
  stat =  "mad") +
  labs(title = "MAD") +
  theme_bw() +
  legend_none()


ggplot() +
  geom_density(aes(x = deviance2A))


ppd_mse_plot2 <- ppc_stat(
  Y2test, PPDcomb2, 
  stat = "mse") +
  labs(title = "MSE") +
  theme_bw() +
  legend_none()

ppd_mean_plot2 <- ppc_stat(
  Y2test, PPDcomb2, 
  stat = "mean") +
  labs(title = "Mean") +
  theme_bw() +
  legend_none()

ppd_sd_plot2 <- ppc_stat(
  Y2test, PPDcomb2, 
  stat = "sd") +
  labs(title = "Standard Deviation") +
  theme_bw() +
  legend_none()

ppd_lay2 <- rbind(c(1,1),
                  c(2,3),
                  c(4,5),
                  c(6,7)
)
bayesplot_grid(
  plots = list(
    ppd_density_plot2,
    ppd_q2.5_plot2,
    ppd_q97.5_plot2,
    ppd_median_plot2,
    ppd_mad_plot2,
    ppd_mean_plot2,
    ppd_sd_plot2),
  grid_args = list(
    layout_matrix = ppd_lay2
  )
)


## Model Comparison ----
### DIC ----
# dic2   <- dic.samples(model2, n.iter=2000, progress.bar="text")
# 
# #### SAVE CHECKPOINT ----
# save(dic2, file = "Model 2 DIC.RData")
# 
# ### WAIC ----
# waic2   <- coda.samples(model2, 
#                         variable.names=c("like"), 
#                         n.iter=2000, progress.bar="text")
# 
# #### SAVE CHECKPOINT ----
# save(waic2, file = "Model 2 WAIC.RData")
# 
# #### Manually Compute WAIC and P ----
# like2   <- waic2[[1]]
# fbar2   <- colMeans(like2)
# P2      <- sum(apply(log(like2),2,var))
# WAIC2   <- -2*sum(log(fbar2))+2*P2
# 
# ### Output Results ----
# dic2
# WAIC2
# P2

## FINAL SAVE ----
save(list = setdiff(
  ls(.GlobalEnv), 
  c("bleaching_data")),
  file = filename2
)

save(
  final_data3,
  samples2,
  paramSamps2,
  ppcSamps2,
  YppcSamps2,
  ppdSamps2,
  quantiles2,
  stats2,
  Y2,
  Y2test,
  Y2train,
  X2,
  X2train,
  X2test,
  trainIndex,
  file = filename2
)

# Delete the other previous RData files now to push to GitHub
file.remove("Model 2 DIC.RData")
file.remove("Model 2 WAIC.RData")

# Remove data to free up space in environment (we can load it later to compare models)
rm(list=setdiff(ls(), c("bleaching_data", "final_data3")))
