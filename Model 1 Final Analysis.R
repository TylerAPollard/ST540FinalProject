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
Y <- final_data3$Percent_Bleaching_Open/100
Y1 <- log(Y) 
Z1 <- ifelse(Y == 0, 0, 1)

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
X1unscale <- X1
# X1$Exposure <- ifelse(X1$Exposure == "Sheltered", 0, 1)
X1 <- scale(X1)

#### Split Data ----
# Observed data
Y1train <- Y1[trainIndex]
Z1train <- Z1[trainIndex]
X1train <- X1[trainIndex,]

# Set aside for prediction
Y1test <- Y1[-trainIndex]
Z1test <- Z1[-trainIndex]
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
burn1     <- 200
n.iters1   <- 400
thin1     <- 5

## Define Model ----
model_string1 <- textConnection("model{
  # Likelihood
  for(i in 1:n){
    Y[i] ~ dnorm(mu[i], tau)
    mu[i] <- ifelse(Z[i], alpha + inprod(X[i,],beta[]) - 0.5/tau, 0)
    
    Z[i] ~ dbern(q[i])
    logit(q[i]) <- alphaZ + inprod(X[i,], betaZ[])
    
    # For WAIC
    # like[i] <- dnorm(Y[i],mu[i],tau)
  } 
    
  # Priors  
  for(j in 1:p){ 
    beta[j] ~ dnorm(0,0.01)
    betaZ[j] ~ dnorm(0,0.01) 
  } 
    
  alpha ~  dnorm(0, 0.01)
  alphaZ ~ dnorm(0, 0.01)
  tau   ~  dgamma(0.1, 0.1) 
  sigma <- 1/sqrt(tau)
    
  # Posterior Predicitve Checks
  for(i in 1:n){
    Yppc[i] ~ dnorm(alpha + inprod(X[i,],beta[]), tau)
    Y2A[i] <- exp(Yppc[i])
  }
  D1[1] <- min(Yppc[])
  D1[2] <- max(Yppc[])
  D1[3] <- max(Yppc[]) - min(Yppc[])
  D1[4] <- mean(Yppc[])
  D1[5] <- sd(Yppc[])
  D2[1] <- min(Y2A[])
  D2[2] <- max(Y2A[])
  D2[3] <- max(Y2A[]) - min(Y2A[])
  D2[4] <- mean(Y2A[])
  D2[5] <- sd(Y2A[])
  
  # Predictions
  # for(i in 1:n_pred){
  #   Ypred[i] ~ dnorm(mu_pred[i], tau)
  #   mu_pred[i] <- ifelse(Zpred[i], alpha + inprod(Xpred[i,],beta[]), 0)
  #   
  #   Zpred[i] ~ dbern(qpred[i])
  #   logit(qpred[i]) <- alphaZ + inprod(Xpred[i,], betaZ[])
  # }
}")

### Compile Model Inputs ----
data1   <- list(Y = Y1train,
                X = X1train,
                n = n1train,
                p = p1,
                #Xpred = X1test,
                #n_pred = n1test,
                Z = Z1train
)
params1 <- c("alpha", "beta", "sigma", "D1", "D2", "Yppc", "Y2A",
             "alphaZ", "betaZ")

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
samples1 <- coda.samples(model1, variable.names=params1, n.iter=n.iters1, 
                         n.thin=thin1, progress.bar="text")
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

save(samples1, file = "Model 1 Samples.RData")

# Use this if R session terminates
# load(file = "Model 1 Fit Data.RData")

### Parse Samples by Variable ----
samples1A <- samples1[[1]]
samples1B <- samples1[[2]]

#### Parameters ----
paramSamps1A <- samples1A[ ,2400:2409]
paramSamps1B <- samples1B[ ,2400:2409]
colnames(paramSamps1A) <- c("Intercept", colnames(X1), "sigma")
colnames(paramSamps1B) <- c("Intercept", colnames(X1), "sigma")
paramSamps1 <- mcmc.list(paramSamps1A,paramSamps1B)
paramSamps1Comb <- rbind(paramSamps1A,paramSamps1B)
class(paramSamps1Comb)
str(paramSamps1Comb)
dim(paramSamps1Comb)
dimnames(paramSamps1Comb)

#### PPCs ----
##### Data ----
YppcSamps1A <- samples1A[ ,6:1803]
YppcSamps1B <- samples1B[ ,6:1803]
YppcSamps1 <- mcmc.list(YppcSamps1A,YppcSamps1B)
YppcSamps1Comb <- rbind(YppcSamps1A,YppcSamps1B)
class(YppcSamps1Comb)
str(YppcSamps1Comb)
dim(YppcSamps1Comb)
dimnames(YppcSamps1Comb)

##### Checks ----
ppcSamps1A <- samples1A[ ,1:5]
ppcSamps1B <- samples1B[ ,1:5]
DPrintnames <- c(
  "Min of Y",
  "Max of Y",
  "Range of Y",
  "Mean of Y",
  "SD of Y"
)
colnames(ppcSamps1A) <- DPrintnames
colnames(ppcSamps1B) <- DPrintnames
ppcSamps1 <- mcmc.list(ppcSamps1A,ppcSamps1B)
ppcSamps1Comb <- rbind(ppcSamps1A,ppcSamps1B)

ppcSamps1C <- samples1A[ ,6:10]
ppcSamps1D <- samples1B[ ,6:10]
DPrintnames <- c(
  "Min of Y",
  "Max of Y",
  "Range of Y",
  "Mean of Y",
  "SD of Y"
)
colnames(ppcSamps1C) <- DPrintnames
colnames(ppcSamps1D) <- DPrintnames
ppcSamps1B <- mcmc.list(ppcSamps1C,ppcSamps1D)
ppcSamps1BComb <- rbind(ppcSamps1C,ppcSamps1D)

#### PPDs ----
ppdSamps1A <- samples1A[ ,1804:2399]
ppdSamps1B <- samples1B[ ,1804:2399]
ppdSamps1 <- mcmc.list(ppdSamps1A, ppdSamps1B)
ppdSamps1Comb <- rbind(ppdSamps1A,ppdSamps1B)

### Convergence Diagnostics ----
#### Trace Plots ----
# Trace and Density plots
color_scheme_set("teal")
mcmc_combo(paramSamps1,
           combo = c("trace", "dens_overlay"),
           # pars = c("alpha", 
           #          paste0("beta[", 1:p2, "]"),
           #          "phi"),
           widths = c(2,1))

# Remove plots if it is bogging down environment
dev.off()

#### Effective Sample Size ----
# Checking the effective sample sizes.
# All effective sample sizes are very high, well over 10,000!
effectiveSize(paramSamps1)

#### Gelman-Rubin Diagnostics ----
# R less than 1.1 indicates convergence.
gelman.diag(paramSamps1)

#### Geweke Diagnostics ----
# abs(z) less than 2 indicates convergence.
geweke.diag(paramSamps1)

### Model Summaries ----
summary1 <- summary(paramSamps1)

#### Parameter Estimates ----
stats1 <- summary1$statistics
stats1

#### Parameter 95% CIs ----
# CI = Credible Interval
quantiles1 <- summary1$quantiles
quantiles1

#### Plot parameter distributions ----
mcmc_areas(
  paramSamps1,
  pars = colnames(X1),
  point_est = "mean",
  prob = 0.95) +
  labs(
    title = "Posterior Distributions of Predictors",
    subtitle = "95% Credible Interval about the Point Estimate"
  ) +
  vline_0(linewidth = 1.5) +
  theme_bw()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Goodness of Fit Checks ----
### Posterior Predictive Checks ----
#### Old Plots ----
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

D1C <- rbind(D1A,D1B)

# Create empty vectors to store Bayesian p-values for each chain
pval1A <- rep(0, 5)
names(pval1A) <- DPrintnames
pval1B <- rep(0, 5)
names(pval1B) <- DPrintnames
pval1C <- rep(0, 5)
names(pval1C) <- DPrintnames

#### Plot all PPCs ----
#### and calculate p-values
# May end up converting this to ggplot for report

# For all plots in one
par(mfrow = c(3,2))

# For individual plots
# dev.off()

for(j in 1:5){
  plot(density(D1A[,j]), xlab = "D", ylab = "Posterior Probability", 
       xlim = c(min(D1A[,j], D1B[,j], D1C[,j], D0[j]), 
                max(D1A[,j], D1B[,j], D1C[,j], D0[j])), 
       main = DPrintnames[j])
  lines(density(D1B[,j]), col = "blue")
  lines(density(D1C[,j]), col = "red")
  abline(v = D0[j], col = "green", lwd = 1)
  legend("topleft", c("D1A", "D1B", "D1C", "Observed"), 
         col = c("black", "blue", "red", "green"), lwd = 1)
  
  pval1A[j] <- mean(D1A[,j] > D0[j])
  pval1B[j] <- mean(D1B[,j] > D0[j])
  pval1C[j] <- mean(D1C[,j] > D0[j])
}
pval1A
pval1B
pval1C

# Remove plots if it is bogging down environment
dev.off()

#### Bayesplot PPCs ----
#### and calculate p-values
# May end up converting this to ggplot for report
ppcL1 <- apply(YppcSamps1Comb, 1, function(x) quantile(x, 0.025))
ppcU1 <- apply(YppcSamps1Comb, 1, function(x) quantile(x, 0.975))
ppcMedian1 <- apply(YppcSamps1Comb, 1, median)
ppcMean1 <- apply(YppcSamps1Comb, 1, mean)
ppcSD1 <- apply(YppcSamps1Comb, 1, sd)

DppcL1 <- quantile(Y1train, 0.025)
DppcU1 <- quantile(Y1train, 0.975)
DppcMedian1 <- median(Y1train)
DppcMean1 <-mean(Y1train)
DppcSD1 <- sd(Y1train)

##### Overall Density ----
draw1 <- sample(1:20000, 8000)
ppc_density_plot1 <- 
  ppc_dens_overlay(Y1train, YppcSamps1Comb[draw1,]) +
  labs(title = "Posterior Predictive Checks of Linear Regression on Training Data",
       subtitle = "Simulated Data Sets Compared to Training Data") +
  theme_bw() +
  legend_none() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.5)),
    plot.subtitle = element_text(size = rel(1)))

##### Q2.5%  ----
ppc_q2.5_plot1 <- 
  ppc_stat(Y1train, YppcSamps1Comb, stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = "2.5% Quantile") +
  theme_bw() +
  legend_none()

ppcL1dens <- density(ppcL1)
ppcL1dens <- cbind(ppcL1dens$x, ppcL1dens$y)
ppcL1densB <- ppcL1dens[between(ppcL1dens[,1], quantile(ppcL1, 0.025), quantile(ppcL1, 0.975)), ] 
ppc_q2.5_plot1B <- ggplot() +
  # geom_histogram(aes(x = ppcL2,  after_stat(density)),
  #                fill = "#bcdcdc", color = "#99c7c7") +
  geom_ribbon(aes(x = ppcL1densB[,1], ymin = 0, ymax = ppcL1densB[,2]),
              fill = "#bcdcdc") +
  geom_density(aes(x = ppcL1), color = "#99c7c7", linewidth = .75) +
  #geom_vline(aes(xintercept = quantile(ppcL2, 0.975)), color = "#007C7C", linewidth = 2) +
  geom_vline(aes(xintercept = DppcL1), color = "#007C7C", linewidth = 1) +
  #geom_text(aes(label = paste0("p-value = ", round(mean(ppcL1 > DppcL1), 3))),
            #x = 0.84*max(ppcL1), y = max(ppcL1dens[,2]), size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  labs(title = "2.5% Quantile",
       x = NULL,
       y = "Posterior Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1)))
ppc_q2.5_plot1B

##### Q97.5%  ----
ppc_q97.5_plot1 <- 
  ppc_stat(Y1train, YppcSamps1Comb, stat = function(y) quantile(y, 0.975)) +
  labs(title = "97.5% Quantile") +
  theme_bw() +
  legend_none()

ppcU1dens <- density(ppcU1)
ppcU1dens <- cbind(ppcU1dens$x, ppcU1dens$y)
ppcU1densB <- ppcU1dens[between(ppcU1dens[,1], quantile(ppcU1, 0.025), quantile(ppcU1, 0.975)), ] 
ppc_q97.5_plot1B <- ggplot() +
  # geom_histogram(aes(x = ppcU1,  after_stat(density)),
  #                fill = "#bcdcdc", color = "#99c7c7") +
  geom_ribbon(aes(x = ppcU1densB[,1], ymin = 0, ymax = ppcU1densB[,2]),
              fill = "#bcdcdc") +
  geom_density(aes(x = ppcU1), color = "#99c7c7", linewidth = .75) +
  #geom_vline(aes(xintercept = quantile(ppcU1, 0.975)), color = "#007C7C", linewidth = 1) +
  geom_vline(aes(xintercept = DppcU1), color = "#007C7C", linewidth = 1) +
  #geom_text(aes(label = paste0("p-value = ", round(mean(ppcU1 > DppcU1),3))),
            #x = 0.82*max(ppcU1), y = max(ppcU1dens[,1]), size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  labs(title = "97.5% Quantile",
       x = NULL,
       y = "Posterior Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1)))
ppc_q97.5_plot1B

##### Median ----
ppc_median_plot1 <- 
  ppc_stat(Y1train, YppcSamps1Comb, stat = "median") +
  labs(title = "Median") +
  theme_bw() +
  legend_none()

ppcMedian1dens <- density(ppcMedian1)
ppcMedian1dens <- cbind(ppcMedian1dens$x, ppcMedian1dens$y)
ppcMedian1densB <- ppcMedian1dens[between(ppcMedian1dens[,1], 
                                          quantile(ppcMedian1, 0.025), 
                                          quantile(ppcMedian1, 0.975)), ] 
ppc_median_plot1B <- ggplot() +
  # geom_histogram(aes(x = ppcMedian1,  after_stat(density)),
  #                fill = "#bcdcdc", color = "#99c7c7") +
  geom_ribbon(aes(x = ppcMedian1densB[,1], ymin = 0, ymax = ppcMedian1densB[,2]),
              fill = "#bcdcdc") +
  geom_density(aes(x = ppcMedian1), color = "#99c7c7", linewidth = .75) +
  #geom_vline(aes(xintercept = quantile(ppcMedian1, 0.975)), color = "#007C7C", linewidth = 1) +
  geom_vline(aes(xintercept = DppcMedian1), color = "#007C7C", linewidth = 1) +
  #geom_text(aes(label = paste0("p-value = ", round(mean(ppcMedian1 > DppcMedian1),3))),
            #x = 0.72*max(ppcMedian1), y = max(ppcMedian1dens[,1]), size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  labs(title = "Median",
       x = NULL,
       y = "Posterior Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1)))
ppc_median_plot1B

##### Mean ----
ppc_mean_plot1 <- 
  ppc_stat(Y1train, YppcSamps1Comb, stat = "mean") +
  labs(title = "Mean") +
  theme_bw() +
  legend_none()

ppcMean1dens <- density(ppcMean1)
ppcMean1dens <- cbind(ppcMean1dens$x, ppcMean1dens$y)
ppcMean1densB <- ppcMean1dens[between(ppcMean1dens[,1], 
                                      quantile(ppcMean1, 0.025), 
                                      quantile(ppcMean1, 0.975)), ] 
ppc_mean_plot1B <- ggplot() +
  # geom_histogram(aes(x = ppcMean1,  after_stat(density)),
  #                fill = "#bcdcdc", color = "#99c7c7") +
  geom_ribbon(aes(x = ppcMean1densB[,1], ymin = 0, ymax = ppcMean1densB[,2]),
              fill = "#bcdcdc") +
  geom_density(aes(x = ppcMean1), color = "#99c7c7", linewidth = .75) +
  #geom_vline(aes(xintercept = quantile(ppcMean1, 0.975)), color = "#007C7C", linewidth = 1) +
  geom_vline(aes(xintercept = DppcMean1), color = "#007C7C", linewidth = 1) +
  #geom_text(aes(label = paste0("p-value = ", round(mean(ppcMean1 > DppcMean1),3))),
            #x = 0.96*max(ppcMean1), y = max(ppcMean1dens[,1]), size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  labs(title = "Mean",
       x = NULL,
       y = "Posterior Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1)))
ppc_mean_plot1B

##### Standard Deviation ----
ppc_sd_plot1 <- 
  ppc_stat(Y1train, YppcSamps1Comb, stat = "sd") +
  labs(title = "Standard Deviation") +
  theme_bw() +
  legend_none()

ppcSD1dens <- density(ppcSD1)
ppcSD1dens <- cbind(ppcSD1dens$x, ppcSD1dens$y)
ppcSD1densB <- ppcSD1dens[between(ppcSD1dens[,1], quantile(ppcSD1, 0.025), quantile(ppcSD1, 0.975)), ] 
ppc_sd_plot1B <- ggplot() +
  # geom_histogram(aes(x = ppcSD1,  after_stat(density)),
  #                fill = "#bcdcdc", color = "#99c7c7") +
  geom_ribbon(aes(x = ppcSD1densB[,1], ymin = 0, ymax = ppcSD1densB[,2], linetype = "A"),
              fill = "#bcdcdc") +
  geom_density(aes(x = ppcSD1), color = "#99c7c7", linewidth = .75) +
  #geom_vline(aes(xintercept = quantile(ppcSD1, 0.975)), color = "#007C7C", linewidth = 1) +
  geom_vline(aes(xintercept = DppcSD1, linetype = "B"), color = "#99c7c7", alpha = 0, linewidth = 0.75) +
  geom_vline(aes(xintercept = DppcSD1, linetype = "C"), color = "#007C7C", linewidth = 1) +
  geom_text(aes(label = paste0("p-value = ", round(mean(ppcSD1 > DppcSD1),3))),
            x = 0.965*max(ppcSD1), y = max(ppcSD1dens[,1]), size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  labs(title = "Standard Deviation",
       x = NULL,
       y = "Posterior Density") +
  scale_linetype_manual(
    name = element_blank(),
    values = c(1,1,1),
    breaks = c("C", "B", "A"),
    labels = c("Observed", "Posterior", "95% Credible Interval")
  ) +
  guides(
    linetype = guide_legend(
      override.aes = list(
        alpha = c(1,1,1),
        shape = c(NA, "--", NA)
      )
    )
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1)))
ppc_sd_plot1B
ppc_legend1b <- ggpubr::get_legend(ppc_sd_plot1B)

#### Combination Plot ----
ppc_lay1 <- rbind(c(NA,NA, rep(1,8),NA, NA),
                  c(rep(2,4),rep(3,4),rep(4,4)),
                  c(NA, NA, rep(5,4), rep(6,6))
)

bayesplot_grid(
  plots = list(
    ppc_density_plot1,
    ppc_q2.5_plot1B,
    ppc_median_plot1B,
    ppc_q97.5_plot1B,
    ppc_mean_plot1B,
    ppc_sd_plot1B),
  grid_args = list(
    layout_matrix = ppc_lay1
  )
)

# Remove saved plots for saving
rm(ppc_density_plot2, ppc_q2.5_plot2, ppc_median_plot2, ppc_q97.5_plot2,
   ppc_mean_plot2, ppc_sd_plot2)

### Plot PDDs ----
ppdL1 <- apply(ppdSamps1Comb, 1, function(x) quantile(x, 0.025))
ppdU1 <- apply(ppdSamps1Comb, 1, function(x) quantile(x, 0.975))
ppdMedian1 <- apply(ppdSamps1Comb, 1, median)
ppdMean1 <- apply(ppdSamps1Comb, 1, mean)
ppdSD1 <- apply(ppdSamps1Comb, 1, sd)

DppdL1 <- quantile(Y1test, 0.025)
DppdU1 <- quantile(Y1test, 0.975)
DppdMedian1 <- median(Y1test)
DppdMean1 <-mean(Y1test)
DppdSD1 <- sd(Y1test)

#### Predictive Statisitcs ----
YppdL1 <- apply(ppdSamps1Comb, 2, function(x) quantile(x, 0.025))
YppdU1 <- apply(ppdSamps1Comb, 2, function(x) quantile(x, 0.975))
YppdMedian1 <- apply(ppdSamps1Comb, 2, median)
YppdMean1 <- apply(ppdSamps1Comb, 2, mean)
YppdSD1 <- apply(ppdSamps1Comb, 2, sd)

Y_mod1 <- rbeta(nrow(paramSamps1Comb), 
                paramSamps1Comb[,1] + X1test)


Y1_mean <- c()
Y1_median <- c()
Y1_low <- c()
Y1_high <- c()
expit <- function(x){1/(1+exp(-x))}
for(i in 1:n1test){
  mu <- expit(paramSamps1Comb[,1] + paramSamps1Comb[,2:6]%*%X1test[i,])
  Y_mod1 <- rbeta(nrow(paramSamps1Comb), 
                  paramSamps1Comb[,7]*mu,
                  paramSamps1Comb[,7]*(1-mu))
  Y1_mean[i] <- mean(Y_mod1)
  Y1_median[i] <- median(Y_mod1)
  Y1_low[i] <- quantile(Y_mod1, 0.025)
  Y1_high[i] <- quantile(Y_mod1, 0.975)
}

BIAS1  <- mean(Y1_mean-Y1test)
MSE1   <- mean((Y1_mean-Y1test)^2)
MAE1   <- mean(abs(Y1_mean-Y1test))
MAD1   <- mean(abs(Y1_median-Y1test))
COV1   <- mean( (Y1_low <= Y1test) & (Y1test <= Y1_high))
WIDTH1 <- mean(Y1_high-Y1_low)

predStats1A <- c(
  "BIAS" = BIAS1,
  "MSE" = MSE1, 
  "MAE" = MAE1,
  "MAD" = MAD1,
  "COV" = COV1,
  "WIDTH" = WIDTH1
)

BIAS1B  <- mean(YppdMean1-Y1test)
MSE1B   <- mean((YppdMean1-Y1test)^2)
MAE1B   <- mean(abs(YppdMean1 -Y1test))
MAD1B   <- mean(abs(YppdMedian1-Y1test))
COV1B   <- mean( (YppdL1 <= Y1test) & (Y1test <= YppdU1))
WIDTH1B <- mean(YppdU1-YppdL1)

predStats1B <- c(
  "BIAS" = BIAS1B,
  "MSE" = MSE1B, 
  "MAE" = MAE1B,
  "MAD" = MAD1B,
  "COV" = COV1B,
  "WIDTH" = WIDTH1B
)
predStats1B
predStats1 <- rbind(predStats1A, predStats1B)

#### Bayesplot PPD ----
##### Overall Density ----
ppd_density_plot1 <- 
  ppc_dens_overlay(Y1test, ppdSamps1Comb[,]) +
  labs(title = "Posterior Predictive Distribution of Beta Regression on Test Data",
       subtitle = "Simulated Data Sets Compared to Test Data",
       y = "Posterior Density") +
  annotate("text", label = paste0("BIAS   = ", round(BIAS1B,4)), x = 0.865, y = 4.75, size = 4.5,) +
  annotate("text", label = paste0("MSE    = ", round(MSE1B,4)), x = 0.865, y = 4.5, size = 4.5) +
  annotate("text", label = paste0("MAE    = ", round(MAE1B,4)), x = 0.865, y = 4.25, size = 4.5) +
  annotate("text", label = paste0("MAD    = ", round(MAD1B,4)), x = 0.865, y = 4, size = 4.5) +
  annotate("text", label = paste0("COV    = ", round(COV1B,4)), x = 0.865, y = 3.75, size = 4.5) +
  annotate("text", label = paste0("WIDTH = ", round(WIDTH1B,4)), x = 0.865, y = 3.5, size = 4.5) +
  theme_bw() +
  legend_none() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.5)),
    plot.subtitle = element_text(size = rel(1)))
ppd_density_plot1

##### Coverage ----
ppdOrder1 <- order(Y1test)
ppdXOrder1 <- sort(Y1test)
ppd_cov1_plot <- ggplot() +
  geom_ribbon(aes(x = 1:n1test, ymin = YppdL1, ymax = YppdU1), 
              fill = "#bcdcdc") +
  geom_path(aes(x = 1:n1test, y = Y1test), color = "#007C7C") +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  labs(
    title = "Model Coverage of Test Data",
    x = "Index",
    y = "Percent Bleaching"
  ) +
  geom_text(aes(label = paste0("Coverage = ", round(predStats1B["COV"], 3))),
            x = 75, y = 0.97, size = 4.5) +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1))
  )
ppd_cov1_plot

ppd_lay1 <- rbind(c(1),
                  c(2))

#### Combination Plot ----
ppdComb_plot1 <- bayesplot_grid(
  plots = list(
    ppd_density_plot1,
    ppd_cov1_plot
  ),
  grid_args = list(
    nrow = 2
  )
)

## Model Comparison ----
### DIC ----
# dic1   <- dic.samples(model1, n.iter=n.iter, progress.bar="text")
# 
#### SAVE CHECKPOINT ----
# save(dic1, file = "Model 1 DIC.RData")
# 
### WAIC ----
# waic1   <- coda.samples(model1, 
#                         variable.names=c("like"), 
#                         n.iter=n.iter, progress.bar="text")
# 
#### SAVE CHECKPOINT ----
# save(waic1, file = "Model 1 WAIC.RData")
# 
#### Manually Compute WAIC and P ----
# like1   <- waic1[[1]]
# fbar1   <- colMeans(like1)
# P1      <- sum(apply(log(like1),2,var))
# WAIC1   <- -2*sum(log(fbar1))+2*P1
# 
### Output Results ----
# dic1
# WAIC1
# P1

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