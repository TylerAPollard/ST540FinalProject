#### ST540 Final Project Analysis
## Programmed by: Tyler Pollard, Hanan Ali, Rachel Hardy
## Date Created: 4 April 2024
## Date Modified: 28 April 2024

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

## Create training/test index vector ----
## for CV and testing model performance and fit
set.seed(52)
trainIndex <- createDataPartition(final_data3$Percent_Bleaching,
                                  p = 0.75,
                                  list = FALSE)
trainIndex <- as.vector(trainIndex)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Model 1: Simple Linear Model =====
## Programmer: Rachel
## Modeled with Uninformative Gaussian Priors

## Load Data ----
### Closed Support Data ----
## Y1 is original data with closed support [0,100]
Y1 <- final_data3$Percent_Bleaching
mean(Y1)
sd(Y1)
zeroindex <- Y1 == 0

## We may not need to use these but could be interesting to model data
## separately like in a hierarchical model
Y1zero <- Y1[zeroindex]
Y1nonzero <- Y1[!zeroindex]

### Open Support Data ----
## Change 0's to 0.00001 and 1's to 0.99999 to avoid infinite density
## Y1B is original data to have open support (0,100) 
Y1B <- ifelse(Y1 == 0, 0.001, 
              ifelse(Y1 == 100, 99.999, Y1))
mean(Y1B)
sd(Y1B)
Y1Bzero <- Y1B[zeroindex]
Y1Bnonzero <- Y1B[!zeroindex]

### Log Transformed Data ----
## Y1log is the log of the open supported Y
Y1log <- log(Y1B)
mean(Y1log)
sd(Y1log)
Y1logzero <- Y1log[zeroindex]
Y1lognonzero <- Y1log[!zeroindex]

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
    Percent_Bleaching
  ))
# X1$Exposure <- ifelse(X1$Exposure == "Sheltered", 0, 1)
X1 <- scale(X1)

#### Split Data ----
Y1train <- Y1[trainIndex]
Y1test <- Y1[-trainIndex]
X1train <- X1[trainIndex,]
X1test <- X1[-trainIndex,]

## Simulation Variables ----
n1 <- length(Y1)
n1train <- length(Y1train)
p1 <- ncol(X1train)

## Test models with this these simulation variables:
# burn <- 500
# n.iter <- 1000
# thin <- 5
## We will increase to final model
burn     <- 10000
n.iter   <- 50000
thin     <- 10

## Define Model ----
model_string1 <- textConnection("model{
  # Likelihood
  for(i in 1:n){
    Y[i] ~ dnorm(mu[i],tau)
    mu[i] <- alpha + inprod(X[i,],beta[])
    
    # For WAIC
    like[i] <- dnorm(Y[i],mu[i],tau)
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

### Compile Model Inputs ----
data1   <- list(Y = Y1train,
                X = X1train,
                n = n1train,
                p = p1)
params1 <- c("alpha", "beta", "D1")

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
update(model1, burn, progress.bar="text")
samples1 <- coda.samples(model1, variable.names=params1, n.iter=n.iter, n.thin=thin, 
                         progress.bar="text")
toc()
# Stop here

### SAVE CHECKPOINT ----
# Save RData in case code aborts and causes termination
# Make sure to reload libraries after termination if it occurs
filename1 <- paste0(
  "Model 1 Fit Data ", 
  "(B", burn, "-",
  "I", n.iter, "-",
  "T", thin, 
  ").RData"
)
save(
  list = setdiff(
    ls(.GlobalEnv), 
    c("bleaching_data", "final_data1", "final_data2")
  ),
  file = filename1
)

# Use this if R session terminates
# load(file = "Model 1 Fit Data.RData")

### Convergence Diagnostics ----
#### Trace Plots ----
# Trace and Density plots
# color_scheme_set("mix-blue-red")
# mcmc_trace(samples1)
# mcmc_rank_overlay(samples1)
# dimnames(samples1)
par(mar=c(1,1,1,1))
plot(samples1)

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
  min(Y1),
  max(Y1),
  max(Y1) - min(Y1),
  mean(Y1),
  sd(Y1)
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
update(model, 10000, progress.bar="none")
samp1_cv <- coda.samples(model, 
                     variable.names=c("alpha", "beta", "Y1p"), 
                     n.iter=20000, progress.bar="text")

# Summary
summary1_cv <- summary(samp[,-c(1:n1p)])
summary1_cv

# Extract the samples for each parameter
samps1       <- samp1_cv[[1]]
Yp.samps1    <- samps1[,1:n1p] 
alpha.samps1 <- samps1[,n1p+1]
beta.samps1  <- samps1[,n1p+1+1:p1]
sigma.samps1  <- samps1[,ncol(samps1)]

# Beta means
beta.mn1  <- colMeans(beta.samps)
beta.mn1

# Sigma mean
sigma.mn1 <- mean(sigma.samps)
sigma.mn1

# Alpha (intercept) mean
alpha.mn1 <- mean(alpha.samps)
alpha.mn1

# Graphical representation of plug-in vs PPD vs truth
# Only graphing the first few!
for(j in 1:10){
  
  # Plug-in
  mu <- alpha.mn1+sum(X1p[j,]*beta.mn1)
  y  <- rnorm(20000,mu,sigma.mn1)
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
est1 <- data.frame()

# Looping through and taking the abs difference of the two vectors
for(i in 1:n1p){
  diff <- abs(Y1p[i]-df_mu1[i])
  est1  <- rbind(est1, diff)
}

# First column of the est1 data frame has wonky name
mean_abs_error <- mean(est1$X5.43069500446891)
mean_abs_error

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Model 2: Beta Regression Model (Hanan) ----
## Following most of Rachel's code but with slight modifications to accommodate for the new model
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


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Model 3: Exponential Regression Model (Tyler) ----
## Modeled with Uninformative Gaussian Priors
final_tyler <- final_data3
Y3 <- final_tyler$Percent_Bleaching
mean(Y3)
sd(Y3)
zeroindex <- Y3 == 0
Y3zero <- Y3[zeroindex]
Y3nonzero <- Y3[!zeroindex]
plot(density(Y3zero), xlim = c(0,100))
lines(density(Y3nonzero))
plot(density(1/Y3nonzero))

# Change 0's to 0.00001 and 1's to 0.99999 to avoid infinite density
# Log link
Y3B <- ifelse(Y3 == 0, 0.001, 
              ifelse(Y3 == 100, 99.999, Y3))
mean(Y3B)
sd(Y3B)
Y3B <- log(Y3B)
plot(density(Y3B))

Y3Bzero <- Y3B[zeroindex]
Y3Bnonzero <- Y3B[!zeroindex]
plot(density(Y3Bzero))
plot(density(Y3Bnonzero))

# Inverse link
Y3C <- ifelse(Y3 == 0, 0.001, 
              ifelse(Y3 == 100, 99.999, Y3))
Y3C <- 1/(Y3C)
plot(density(Y3C))
curve(dnorm(x, 0, 66), xlim = c(-400,400))

Y3Czero <- Y3C[zeroindex]
Y3Cnonzero <- Y3C[!zeroindex]
plot(density(Y3Czero))
plot(density(Y3Cnonzero))

Z3 <- ifelse(Y3 == 0, 0, 1)

# X3 <- subset(final_hanan, select = -c(Date, Date_Year, Exposure, Percent_Bleaching))
# X3 <- subset(X3, select = -c(Turbidity, SSTA, TSA, Windspeed))
X3 <- final_tyler |> 
  select(-c(
    Date,
    Exposure,
    Date_Year,
    Turbidity,
    SSTA,
    TSA,
    Windspeed,
    Percent_Bleaching
  ))
X3 <- scale(X3)

n3 <- length(Y3)
p3 <- ncol(X3)

burn     <- 100
n.iter   <- 200
thin     <- 1

model_string3A <- textConnection("model{
    # Likelihood
    for(i in 1:n){
      Z[i] ~ dbinom(q[i], 1)
      logit(q[i]) <- int + inprod(X[i,], alpha[])
    } 
    
    # Priors  
    for(j in 1:p){ 
      alpha[j] ~ dnorm(0, 0.01)
    } 
    int ~ dnorm(0, 0.01)
    
    # Posterior Predicitve Checks
    for(i in 1:n){
      Zppc[i] ~ dbinom(q[i], 1)
    }
    # D3[1] <- min(Yppc[])
    # D3[2] <- max(Yppc[])
    # D3[3] <- max(Yppc[]) - min(Yppc[])
    # D3[4] <- mean(Yppc[])
    # D3[5] <- sd(Yppc[])
      
}")

data3   <- list(Z=Z3,X=X3,n=n3,p=p3)
params3 <- c("int", "alpha", "Zppc") #, "D3")

inits <- list(
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 52),
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 8)
)

## Run this all at once to get how long it took
# Start here
tic()
model3A <- jags.model(model_string3A, data=data3, inits = inits,
                      n.chains=2, quiet=TRUE)
#error persists. i've tried increasing burn in. next objective, is to review priors or perhaps try using a different sampler
update(model3A, burn, progress.bar="none")
samples3A <- coda.samples(model3A, variable.names=params3, n.iter=n.iter, n.thin=thin, 
                          progress.bar="none")
toc()

model_string3B <- textConnection("model{
    # Likelihood
    for(i in 1:n){
      # Y[i] ~ dgamma(shape[i], rate[i])
      # mu[i] <- exp(alpha + inprod(X[i,], beta[]))
      # shape[i] <- pow(mu[i], 2)/pow(taue, 2)
      # rate[i] <- mu[i]/pow(taue, 2)
      Y[i] ~ 
      
      # For WAIC
      #like[i]    <- dgamma(Y[i], shape[i], rate[i])
    } 
    
    # Priors  
    for(j in 1:p){ 
      beta[j] ~ dnorm(0, 0.01)
      #b[j] <- exp(beta[j])
    } 
    alpha ~ dnorm(0, 0.01) 
    #a <- exp(alpha)
    taue ~ dgamma(0.1, 0.1) # Shape parameter for beta distribution
    #sigma2 <- 1/taue
    
    # Posterior Predicitve Checks
    for(i in 1:n){
      Yppc[i] ~ dgamma(shape[i], rate[i])
    }
    # D3[1] <- min(Yppc[])
    # D3[2] <- max(Yppc[])
    # D3[3] <- max(Yppc[]) - min(Yppc[])
    # D3[4] <- mean(Yppc[])
    # D3[5] <- sd(Yppc[])
      
}")

data3   <- list(Y=Y3,X=X3,n=n3,p=p3)
params3 <- c("alpha", "beta") #, "D3")

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

### Model 3 with glm
glm1 <- glm(Y3 ~ X3, )

