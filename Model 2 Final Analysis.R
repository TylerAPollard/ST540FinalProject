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
X2unscale <- X2
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
paramSamps2Comb <- rbind(paramSamps2A,paramSamps2B)
class(paramSamps2Comb)
str(paramSamps2Comb)
dim(paramSamps2Comb)
dimnames(paramSamps2Comb)

#### PPCs ----
##### Data ----
YppcSamps2A <- samples2A[ ,6:1803]
YppcSamps2B <- samples2B[ ,6:1803]
YppcSamps2 <- mcmc.list(YppcSamps2A,YppcSamps2B)
YppcSamps2Comb <- rbind(YppcSamps2A,YppcSamps2B)
class(YppcSamps2Comb)
str(YppcSamps2Comb)
dim(YppcSamps2Comb)
dimnames(YppcSamps2Comb)

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
ppcSamps2Comb <- rbind(ppcSamps2A,ppcSamps2B)

#### PPDs ----
ppdSamps2A <- samples2A[ ,1804:2399]
ppdSamps2B <- samples2B[ ,1804:2399]
ppdSamps2 <- mcmc.list(ppdSamps2A, ppdSamps2B)
ppdSamps2Comb <- rbind(ppdSamps2A,ppdSamps2B)

#### Bayesplot Exmaple Data ----
set.seed(52)
# MCMC
ex_mcmc_draws <- example_mcmc_draws(chains = 2, params = 6)
class(ex_mcmc_draws)
str(ex_mcmc_draws)
dim(ex_mcmc_draws)
dimnames(ex_mcmc_draws)

# Observed Y
ex_y_draws <- example_y_data()
class(ex_y_draws)
str(ex_y_draws)

# Pred Y
ex_yrep_draws <- example_yrep_draws()
class(ex_yrep_draws)
str(ex_yrep_draws)
dim(ex_yrep_draws)
dimnames(ex_yrep_draws)

# X
ex_x_draws <- example_x_data()
class(ex_x_draws)
str(ex_x_draws)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Convergence Diagnostics ----
#### Trace Plots ----
# Trace and Density plots
color_scheme_set("teal")
mcmc_combo(paramSamps2,
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

#### Plot parameter distributions ----
color_scheme_set("teal")
mcmc_areas(
  paramSamps2,
  pars = colnames(X2),
  point_est = "mean",
  prob = 0.95) +
  labs(
    title = "Posterior Distributions of Predictors",
    subtitle = "95% Credible Interval about the Point Estimate"
  ) +
  vline_0(linewidth = 1.5) +
  scale_y_discrete(expand = expansion(add = c(0.5,1.25)),
                   breaks = rev(c("Distance_to_Shore",
                              "Cyclone_Frequency",
                              "Date_Year",
                              "TSA",
                              "TSA_DHW"))) +
  annotate("text", x = quantiles2["Distance_to_Shore", 3], y = 2-0.05,
           label = paste0("(", round(quantiles2["Distance_to_Shore", 1], 3), ", ",
                          round(quantiles2["Distance_to_Shore", 5], 3), ")")) +
  annotate("text", x = quantiles2["Cyclone_Frequency", 3], y = 3,
           label = paste0("(", round(quantiles2["Cyclone_Frequency", 1], 3), ", ",
                          round(quantiles2["Cyclone_Frequency", 5], 3), ")")) +
  annotate("text", x = quantiles2["Date_Year", 3], y = 4-0.1,
           label = paste0("(", round(quantiles2["Date_Year", 1], 3), ", ",
                          round(quantiles2["Date_Year", 5], 3), ")")) +
  annotate("text", x = quantiles2["TSA", 3], y = 5-0.1,
           label = paste0("(", round(quantiles2["TSA", 1], 3), ", ",
                          round(quantiles2["TSA", 5], 3), ")")) +
  annotate("text", x = quantiles2["TSA_DHW", 3]-0.04, y = 6-0.1,
           label = paste0("(", round(quantiles2["TSA_DHW", 1], 3), ", ",
                          round(quantiles2["TSA_DHW", 5], 3), ")")) +
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

D2C <- rbind(D2A,D2B)

# Create empty vectors to store Bayesian p-values for each chain
pval2A <- rep(0, 5)
names(pval2A) <- DPrintnames
pval2B <- rep(0, 5)
names(pval2B) <- DPrintnames
pval2C <- rep(0, 5)
names(pval2C) <- DPrintnames

# For all plots in one
par(mfrow = c(3,2))

# For individual plots
#dev.off()

for(j in 1:5){
  plot(density(D2A[,j]), xlab = "D", ylab = "Posterior Probability", 
       xlim = c(min(D2A[,j], D2B[,j], D2C[,j], D0[j]), 
                max(D2A[,j], D2B[,j], D2C[,j], D0[j])), 
       main = DPrintnames[j])
  lines(density(D2B[,j]), col = "blue")
  lines(density(D2C[,j]), col = "red")
  abline(v = D0[j], col = "green", lwd = 2)
  legend("topleft", c("D2A", "D2B", "D2C", "Observed"), 
         col = c("black", "blue", "red", "green"), lwd = 2)
  
  pval2A[j] <- mean(D2A[,j] > D0[j])
  pval2B[j] <- mean(D2B[,j] > D0[j])
  pval2C[j] <- mean(D2C[,j] > D0[j])
}
pval2A
pval2B
pval2C

# Remove plots if it is bogging down environment
dev.off()

#### Bayesplot PPCs ----
#### and calculate p-values
# May end up converting this to ggplot for report
ppcL2 <- apply(YppcSamps2Comb, 1, function(x) quantile(x, 0.025))
ppcU2 <- apply(YppcSamps2Comb, 1, function(x) quantile(x, 0.975))
ppcMedian2 <- apply(YppcSamps2Comb, 1, median)
ppcMean2 <- apply(YppcSamps2Comb, 1, mean)
ppcSD2 <- apply(YppcSamps2Comb, 1, sd)

DppcL2 <- quantile(Y2train, 0.025)
DppcU2 <- quantile(Y2train, 0.975)
DppcMedian2 <- median(Y2train)
DppcMean2 <-mean(Y2train)
DppcSD2 <- sd(Y2train)

##### Overall Density ----
ppc_density_plot2 <- 
  ppc_dens_overlay(Y2train, YppcSamps2Comb) +
  labs(title = "Posterior Predictive Checks of Beta Regression on Training Data",
       subtitle = "Simulated Data Sets Compared to Training Data") +
  theme_bw() +
  legend_none() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.5)),
    plot.subtitle = element_text(size = rel(1)))

##### Q2.5%  ----
ppc_q2.5_plot2 <- 
  ppc_stat(Y2train, YppcSamps2Comb, stat = function(y) quantile(y, 0.025), freq = FALSE) +
  labs(title = "2.5% Quantile") +
  theme_bw() +
  legend_none()

ppcL2dens <- density(ppcL2)
ppcL2dens <- cbind(ppcL2dens$x, ppcL2dens$y)
ppcL2densB <- ppcL2dens[between(ppcL2dens[,1], quantile(ppcL2, 0.025), quantile(ppcL2, 0.975)), ] 
ppc_q2.5_plot2B <- ggplot() +
  # geom_histogram(aes(x = ppcL2,  after_stat(density)),
  #                fill = "#bcdcdc", color = "#99c7c7") +
  geom_ribbon(aes(x = ppcL2densB[,1], ymin = 0, ymax = ppcL2densB[,2]),
              fill = "#bcdcdc") +
  geom_density(aes(x = ppcL2), color = "#99c7c7", linewidth = .75) +
  #geom_vline(aes(xintercept = quantile(ppcL2, 0.975)), color = "#007C7C", linewidth = 2) +
  geom_vline(aes(xintercept = DppcL2), color = "#007C7C", linewidth = 1) +
  geom_text(aes(label = paste0("p-value = ", round(mean(ppcL2 > DppcL2), 3))),
            x = 0.84*max(ppcL2), y = max(ppcL2dens[,2]), size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  labs(title = "2.5% Quantile",
       x = NULL,
       y = "Posterior Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1)))
ppc_q2.5_plot2B

##### Q97.5%  ----
ppc_q97.5_plot2 <- 
  ppc_stat(Y2train, YppcSamps2Comb, stat = function(y) quantile(y, 0.975)) +
  labs(title = "97.5% Quantile") +
  theme_bw() +
  legend_none()

ppcU2dens <- density(ppcU2)
ppcU2dens <- cbind(ppcU2dens$x, ppcU2dens$y)
ppcU2densB <- ppcU2dens[between(ppcU2dens[,1], quantile(ppcU2, 0.025), quantile(ppcU2, 0.975)), ] 
ppc_q97.5_plot2B <- ggplot() +
  # geom_histogram(aes(x = ppcU2,  after_stat(density)),
  #                fill = "#bcdcdc", color = "#99c7c7") +
  geom_ribbon(aes(x = ppcU2densB[,1], ymin = 0, ymax = ppcU2densB[,2]),
              fill = "#bcdcdc") +
  geom_density(aes(x = ppcU2), color = "#99c7c7", linewidth = .75) +
  #geom_vline(aes(xintercept = quantile(ppcU2, 0.975)), color = "#007C7C", linewidth = 2) +
  geom_vline(aes(xintercept = DppcU2), color = "#007C7C", linewidth = 1) +
  geom_text(aes(label = paste0("p-value = ", round(mean(ppcU2 > DppcU2),3))),
            x = 0.82*max(ppcU2), y = max(ppcU2dens[,2]), size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  labs(title = "97.5% Quantile",
       x = NULL,
       y = "Posterior Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1)))
ppc_q97.5_plot2B

##### Median ----
ppc_median_plot2 <- 
  ppc_stat(Y2train, YppcSamps2Comb, stat = "median") +
  labs(title = "Median") +
  theme_bw() +
  legend_none()

ppcMedian2dens <- density(ppcMedian2)
ppcMedian2dens <- cbind(ppcMedian2dens$x, ppcMedian2dens$y)
ppcMedian2densB <- ppcMedian2dens[between(ppcMedian2dens[,1], quantile(ppcMedian2, 0.025), quantile(ppcMedian2, 0.975)), ] 
ppc_median_plot2B <- ggplot() +
  # geom_histogram(aes(x = ppcMedian2,  after_stat(density)),
  #                fill = "#bcdcdc", color = "#99c7c7") +
  geom_ribbon(aes(x = ppcMedian2densB[,1], ymin = 0, ymax = ppcMedian2densB[,2]),
              fill = "#bcdcdc") +
  geom_density(aes(x = ppcMedian2), color = "#99c7c7", linewidth = .75) +
  #geom_vline(aes(xintercept = quantile(ppcMedian2, 0.975)), color = "#007C7C", linewidth = 2) +
  geom_vline(aes(xintercept = DppcMedian2), color = "#007C7C", linewidth = 1) +
  geom_text(aes(label = paste0("p-value = ", round(mean(ppcMedian2 > DppcMedian2),3))),
            x = 0.72*max(ppcMedian2), y = max(ppcMedian2dens[,2]), size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  labs(title = "Median",
       x = NULL,
       y = "Posterior Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1)))
ppc_median_plot2B

##### Mean ----
ppc_mean_plot2 <- 
  ppc_stat(Y2train, YppcSamps2Comb, stat = "mean") +
  labs(title = "Mean") +
  theme_bw() +
  legend_none()

ppcMean2dens <- density(ppcMean2)
ppcMean2dens <- cbind(ppcMean2dens$x, ppcMean2dens$y)
ppcMean2densB <- ppcMean2dens[between(ppcMean2dens[,1], quantile(ppcMean2, 0.025), quantile(ppcMean2, 0.975)), ] 
ppc_mean_plot2B <- ggplot() +
  # geom_histogram(aes(x = ppcMean2,  after_stat(density)),
  #                fill = "#bcdcdc", color = "#99c7c7") +
  geom_ribbon(aes(x = ppcMean2densB[,1], ymin = 0, ymax = ppcMean2densB[,2]),
              fill = "#bcdcdc") +
  geom_density(aes(x = ppcMean2), color = "#99c7c7", linewidth = .75) +
  #geom_vline(aes(xintercept = quantile(ppcMean2, 0.975)), color = "#007C7C", linewidth = 2) +
  geom_vline(aes(xintercept = DppcMean2), color = "#007C7C", linewidth = 1) +
  geom_text(aes(label = paste0("p-value = ", round(mean(ppcMean2 > DppcMean2),3))),
            x = 0.96*max(ppcMean2), y = max(ppcMean2dens[,2]), size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  labs(title = "Mean",
       x = NULL,
       y = "Posterior Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1)))
ppc_mean_plot2B

##### Standard Deviation ----
ppc_sd_plot2 <- 
  ppc_stat(Y2train, YppcSamps2Comb, stat = "sd") +
  labs(title = "Standard Deviation") +
  theme_bw() +
  legend_none()

ppcSD2dens <- density(ppcSD2)
ppcSD2dens <- cbind(ppcSD2dens$x, ppcSD2dens$y)
ppcSD2densB <- ppcSD2dens[between(ppcSD2dens[,1], quantile(ppcSD2, 0.025), quantile(ppcSD2, 0.975)), ] 
ppc_sd_plot2B <- ggplot() +
  # geom_histogram(aes(x = ppcSD2,  after_stat(density)),
  #                fill = "#bcdcdc", color = "#99c7c7") +
  geom_ribbon(aes(x = ppcSD2densB[,1], ymin = 0, ymax = ppcSD2densB[,2], linetype = "A"),
              fill = "#bcdcdc") +
  geom_density(aes(x = ppcSD2), color = "#99c7c7", linewidth = .75) +
  #geom_vline(aes(xintercept = quantile(ppcSD2, 0.975)), color = "#007C7C", linewidth = 2) +
  geom_vline(aes(xintercept = DppcSD2, linetype = "B"), color = "#99c7c7", alpha = 0, linewidth = 0.75) +
  geom_vline(aes(xintercept = DppcSD2, linetype = "C"), color = "#007C7C", linewidth = 1) +
  geom_text(aes(label = paste0("p-value = ", round(mean(ppcSD2 > DppcSD2),3))),
            x = 0.965*max(ppcSD2), y = max(ppcSD2dens[,2]), size = 3) +
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
ppc_sd_plot2B
ppc_legend2b <- ggpubr::get_legend(ppc_sd_plot2B)

#### Combination Plot ----
ppc_lay2 <- rbind(c(NA,NA, rep(1,8),NA, NA),
                  c(rep(2,4),rep(3,4),rep(4,4)),
                  c(NA, NA, rep(5,4), rep(6,6))
)

ppcComb_plot <- bayesplot_grid(
  plots = list(
    ppc_density_plot2,
    ppc_q2.5_plot2B,
    ppc_median_plot2B,
    ppc_q97.5_plot2B,
    ppc_mean_plot2B,
    ppc_sd_plot2B),
  grid_args = list(
    layout_matrix = ppc_lay2
  )
)

# Remove saved plots for saving
rm(ppc_density_plot2, ppc_q2.5_plot2, ppc_median_plot2, ppc_q97.5_plot2,
   ppc_mean_plot2, ppc_sd_plot2)

### Plot PDDs ----
ppdL2 <- apply(ppdSamps2Comb, 1, function(x) quantile(x, 0.025))
ppdU2 <- apply(ppdSamps2Comb, 1, function(x) quantile(x, 0.975))
ppdMedian2 <- apply(ppdSamps2Comb, 1, median)
ppdMean2 <- apply(ppdSamps2Comb, 1, mean)
ppdSD2 <- apply(ppdSamps2Comb, 1, sd)

DppdL2 <- quantile(Y2test, 0.025)
DppdU2 <- quantile(Y2test, 0.975)
DppdMedian2 <- median(Y2test)
DppdMean2 <-mean(Y2test)
DppdSD2 <- sd(Y2test)

#### Predictive Statisitcs ----
YppdL2 <- apply(ppdSamps2Comb, 2, function(x) quantile(x, 0.025))
YppdU2 <- apply(ppdSamps2Comb, 2, function(x) quantile(x, 0.975))
YppdMedian2 <- apply(ppdSamps2Comb, 2, median)
YppdMean2 <- apply(ppdSamps2Comb, 2, mean)
YppdSD2 <- apply(ppdSamps2Comb, 2, sd)

Y_mod2 <- rbeta(nrow(paramSamps2Comb), 
                paramSamps2Comb[,1] + X2test)




Y2_mean <- c()
Y2_median <- c()
Y2_low <- c()
Y2_high <- c()
expit <- function(x){1/(1+exp(-x))}
for(i in 1:n2test){
  mu <- expit(paramSamps2Comb[,1] + paramSamps2Comb[,2:6]%*%X2test[i,])
  Y_mod2 <- rbeta(nrow(paramSamps2Comb), 
                  paramSamps2Comb[,7]*mu,
                  paramSamps2Comb[,7]*(1-mu))
  Y2_mean[i] <- mean(Y_mod2)
  Y2_median[i] <- median(Y_mod2)
  Y2_low[i] <- quantile(Y_mod2, 0.025)
  Y2_high[i] <- quantile(Y_mod2, 0.975)
}

BIAS2  <- mean(Y2_mean-Y2test)
MSE2   <- mean((Y2_mean-Y2test)^2)
MAE2   <- mean(abs(Y2_mean-Y2test))
MAD2   <- mean(abs(Y2_median-Y2test))
COV2   <- mean( (Y2_low <= Y2test) & (Y2test <= Y2_high))
WIDTH2 <- mean(Y2_high-Y2_low)

predStats2A <- c(
  "BIAS" = BIAS2,
  "MSE" = MSE2, 
  "MAE" = MAE2,
  "MAD" = MAD2,
  "COV" = COV2,
  "WIDTH" = WIDTH2
)

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
predStats2 <- rbind(predStats2A, predStats2B)

#### Bayesplot PPD ----
##### Overall Density ----
ppd_density_plot2 <- 
  ppc_dens_overlay(Y2test, ppdSamps2Comb[,]) +
  labs(title = "Posterior Predictive Distribution of Beta Regression on Test Data",
       subtitle = "Simulated Data Sets Compared to Test Data",
       y = "Posterior Density") +
  annotate("text", label = paste0("BIAS   = ", round(BIAS2B,4)), x = 0.865, y = 4.75, size = 4.5,) +
  annotate("text", label = paste0("MSE    = ", round(MSE2B,4)), x = 0.865, y = 4.5, size = 4.5) +
  annotate("text", label = paste0("MAE    = ", round(MAE2B,4)), x = 0.865, y = 4.25, size = 4.5) +
  annotate("text", label = paste0("MAD    = ", round(MAD2B,4)), x = 0.865, y = 4, size = 4.5) +
  annotate("text", label = paste0("COV    = ", round(COV2B,4)), x = 0.865, y = 3.75, size = 4.5) +
  annotate("text", label = paste0("WIDTH = ", round(WIDTH2B,4)), x = 0.865, y = 3.5, size = 4.5) +
  theme_bw() +
  legend_none() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.5)),
    plot.subtitle = element_text(size = rel(1)))
ppd_density_plot2

##### Coverage ----
ppdOrder2 <- order(Y2test)
ppdXOrder2 <- sort(Y2test)
ppd_cov2_plot <- ggplot() +
  geom_ribbon(aes(x = 1:n2test, ymin = YppdL2, ymax = YppdU2), 
              fill = "#bcdcdc") +
  geom_path(aes(x = 1:n2test, y = Y2test), color = "#007C7C") +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  labs(
    title = "Model Coverage of Test Data",
    x = "Index",
    y = "Percent Bleaching"
  ) +
  geom_text(aes(label = paste0("Coverage = ", round(predStats2B["COV"], 3))),
            x = 75, y = 0.97, size = 4.5) +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1))
  )
ppd_cov2_plot

ppd_lay2 <- rbind(c(1),
                  c(2))

#### Combination Plot ----
ppdComb_plot2 <- bayesplot_grid(
  plots = list(
    ppd_density_plot2,
    ppd_cov2_plot
  ),
  grid_args = list(
    nrow = 2
  )
)


trainIndex
testIndex <- 1:2394
testIndex <- testIndex[-trainIndex]

ppdRand <- sample(1:596, 6)
ppdRandDraw <- t(ppdSamps2Comb[ppdRand, ])

final_data_pred2 <- final_data3 |>
  #select(Latitude_Degrees, Longitude_Degrees) |>
  slice(testIndex)
final_data_pred2 <- cbind(final_data_pred2, ppdRandDraw)
colnames(final_data_pred2)[5:10] <- paste0("Random Draw ", 1:6)

long_final_data_pred2 <- rbind(
  final_data_pred2,
  final_data_pred2,
  final_data_pred2,
  final_data_pred2,
  final_data_pred2,
  final_data_pred2
)
long_final_data_pred2$Random_Draw <- rep(1:6, each = n2test)
long_final_data_pred2$Percent_Bleaching_Pred <- c(
  ppdRandDraw[,1],
  ppdRandDraw[,2],
  ppdRandDraw[,3],
  ppdRandDraw[,4],
  ppdRandDraw[,5],
  ppdRandDraw[,6]
)

world_coordinates <- map_data("county") 
ggplot() + 
  # geom_map() function takes world coordinates  
  # as input to plot world map 
  geom_map( 
    data = world_coordinates, map = world_coordinates, 
    aes(x = long, y = lat, map_id = region) 
  ) + 
  geom_point(
    data = long_final_data_pred2 |> slice(1:596) |> filter(Date_Year %in% c(2006, 2011, 2016)),
    aes(x = Longitude_Degrees, y = Latitude_Degrees, 
        color = Percent_Bleaching_Pred)
  ) +
  xlim(c(-83.5,-79.5)) +
  ylim(c(24.25,27.5)) +
  scale_color_continuous(low = "green", high = "red", limits = c(0,1)) +
  labs(
    title = "Posterior Prediction of Percent Bleaching Test Set Locations",
    subtitle = "Comparison of Beginning, Middle, and End Years from Data") +
  facet_wrap(vars(Date_Year), ncol = 3) +
  theme_bw()

##### Prediction vs Year ----
X2unscale <- final_data3 |> 
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
scaleMeans <- apply(X2unscale, 2, mean)
scaleSDs <- apply(X2unscale, 2, sd)
ppd_byCovariates <- cbind(X2unscale[-trainIndex,], YppdMean2, YppdL2, YppdU2)

ppd_byYear2_plot <- ggplot(data = ppd_byCovariates) +
  geom_line(aes(x = Date_Year, y = YppdMean2))
ppd_byYear2_plot

ppd_intervals(ypred = ppdSamps2Comb, x = ppd_byCovariates$Date_Year)

ppc_violin_grouped(y = Y2test, yrep = ppdSamps2Comb, group = ppd_byCovariates$Date_Year) +
  theme_bw() +
  coord_flip()

# ppdL2dens <- density(ppdL2)
# ppdL2dens <- cbind(ppdL2dens$x, ppdL2dens$y)
# 
# ppdU2dens <- density(ppdU2)
# ppdU2dens <- cbind(ppdU2dens$x, ppdU2dens$y)
# 
# ppdMean2dens <- density(ppdMean2)
# ppdMean2dens <- cbind(ppdMean2dens$x, ppdMean2dens$y)
# 
# ppdmu2curve <- 
# 
# ppd_density_plot2B <- ggplot() +
#   geom_density(aes(x = Y2test)) +
#   geom_density(aes(x = ppdSamps2Comb))

##### Q2.5%  ----
# ppd_q2.5_plot2 <- 
#   ppc_stat(Y2test, ppdSamps2Comb, stat = function(y) quantile(y, 0.025), freq = FALSE) +
#   labs(title = "2.5% Quantile") +
#   theme_bw() +
#   legend_none()
# 
# ppdL2dens <- density(ppdL2)
# ppdL2dens <- cbind(ppdL2dens$x, ppdL2dens$y)
# ppdL2densB <- ppdL2dens[between(ppdL2dens[,1], quantile(ppdL2, 0.025), quantile(ppdL2, 0.975)), ] 
# ppd_q2.5_plot2B <- ggplot() +
#   # geom_histogram(aes(x = ppdL2,  after_stat(density)),
#   #                fill = "#bcdcdc", color = "#99c7c7") +
#   geom_ribbon(aes(x = ppdL2densB[,1], ymin = 0, ymax = ppdL2densB[,2]),
#               fill = "#bcdcdc") +
#   geom_density(aes(x = ppdL2), color = "#99c7c7", linewidth = .75) +
#   #geom_vline(aes(xintercept = quantile(ppdL2, 0.975)), color = "#007C7C", linewidth = 2) +
#   geom_vline(aes(xintercept = DppdL2), color = "#007C7C", linewidth = 1) +
#   geom_text(aes(label = paste0("p-value = ", round(mean(ppdL2 > DppdL2), 3))),
#             x = 0.84*max(ppdL2), y = max(ppdL2dens[,2]), size = 3) +
#   scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
#   labs(title = "2.5% Quantile",
#        x = NULL,
#        y = "Posterior Density") +
#   theme_bw() +
#   theme(
#     plot.title = element_text(size = rel(1)))
# ppd_q2.5_plot2B
# 
##### Q97.5%  ----
# ppd_q97.5_plot2 <- 
#   ppc_stat(Y2test, ppdSamps2Comb, stat = function(y) quantile(y, 0.975)) +
#   labs(title = "97.5% Quantile") +
#   theme_bw() +
#   legend_none()
# 
# ppdU2dens <- density(ppdU2)
# ppdU2dens <- cbind(ppdU2dens$x, ppdU2dens$y)
# ppdU2densB <- ppdU2dens[between(ppdU2dens[,1], quantile(ppdU2, 0.025), quantile(ppdU2, 0.975)), ] 
# ppd_q97.5_plot2B <- ggplot() +
#   # geom_histogram(aes(x = ppdU2,  after_stat(density)),
#   #                fill = "#bcdcdc", color = "#99c7c7") +
#   geom_ribbon(aes(x = ppdU2densB[,1], ymin = 0, ymax = ppdU2densB[,2]),
#               fill = "#bcdcdc") +
#   geom_density(aes(x = ppdU2), color = "#99c7c7", linewidth = .75) +
#   #geom_vline(aes(xintercept = quantile(ppdU2, 0.975)), color = "#007C7C", linewidth = 2) +
#   geom_vline(aes(xintercept = DppdU2), color = "#007C7C", linewidth = 1) +
#   geom_text(aes(label = paste0("p-value = ", round(mean(ppdU2 > DppdU2),3))),
#             x = 0.82*max(ppdU2), y = max(ppdU2dens[,2]), size = 3) +
#   scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
#   labs(title = "97.5% Quantile",
#        x = NULL,
#        y = "Posterior Density") +
#   theme_bw() +
#   theme(
#     plot.title = element_text(size = rel(1)))
# ppd_q97.5_plot2B
# 
##### Median ----
# ppd_median_plot2 <- 
#   ppc_stat(Y2test, ppdSamps2Comb, stat = "median") +
#   labs(title = "Median") +
#   theme_bw() +
#   legend_none()
# 
# ppdMedian2dens <- density(ppdMedian2)
# ppdMedian2dens <- cbind(ppdMedian2dens$x, ppdMedian2dens$y)
# ppdMedian2densB <- ppdMedian2dens[between(ppdMedian2dens[,1], quantile(ppdMedian2, 0.025), quantile(ppdMedian2, 0.975)), ] 
# ppd_median_plot2B <- ggplot() +
#   # geom_histogram(aes(x = ppdMedian2,  after_stat(density)),
#   #                fill = "#bcdcdc", color = "#99c7c7") +
#   geom_ribbon(aes(x = ppdMedian2densB[,1], ymin = 0, ymax = ppdMedian2densB[,2]),
#               fill = "#bcdcdc") +
#   geom_density(aes(x = ppdMedian2), color = "#99c7c7", linewidth = .75) +
#   #geom_vline(aes(xintercept = quantile(ppdMedian2, 0.975)), color = "#007C7C", linewidth = 2) +
#   geom_vline(aes(xintercept = DppdMedian2), color = "#007C7C", linewidth = 1) +
#   geom_text(aes(label = paste0("p-value = ", round(mean(ppdMedian2 > DppdMedian2),3))),
#             x = 0.72*max(ppdMedian2), y = max(ppdMedian2dens[,2]), size = 3) +
#   scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
#   labs(title = "Median",
#        x = NULL,
#        y = "Posterior Density") +
#   theme_bw() +
#   theme(
#     plot.title = element_text(size = rel(1)))
# ppd_median_plot2B
# 
##### Mean ----
# ppd_mean_plot2 <- 
#   ppc_stat(Y2test, ppdSamps2Comb, stat = "mean") +
#   labs(title = "Mean") +
#   theme_bw() +
#   legend_none()
# 
# ppdMean2dens <- density(ppdMean2)
# ppdMean2dens <- cbind(ppdMean2dens$x, ppdMean2dens$y)
# ppdMean2densB <- ppdMean2dens[between(ppdMean2dens[,1], quantile(ppdMean2, 0.025), quantile(ppdMean2, 0.975)), ] 
# ppd_mean_plot2B <- ggplot() +
#   # geom_histogram(aes(x = ppdMean2,  after_stat(density)),
#   #                fill = "#bcdcdc", color = "#99c7c7") +
#   geom_ribbon(aes(x = ppdMean2densB[,1], ymin = 0, ymax = ppdMean2densB[,2]),
#               fill = "#bcdcdc") +
#   geom_density(aes(x = ppdMean2), color = "#99c7c7", linewidth = .75) +
#   #geom_vline(aes(xintercept = quantile(ppdMean2, 0.975)), color = "#007C7C", linewidth = 2) +
#   geom_vline(aes(xintercept = DppdMean2), color = "#007C7C", linewidth = 1) +
#   geom_text(aes(label = paste0("p-value = ", round(mean(ppdMean2 > DppdMean2),3))),
#             x = 0.96*max(ppdMean2), y = max(ppdMean2dens[,2]), size = 3) +
#   scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
#   labs(title = "Mean",
#        x = NULL,
#        y = "Posterior Density") +
#   theme_bw() +
#   theme(
#     plot.title = element_text(size = rel(1)))
# ppd_mean_plot2B
# 
##### Standard Deviation ----
# ppd_sd_plot2 <- 
#   ppc_stat(Y2test, ppdSamps2Comb, stat = "sd") +
#   labs(title = "Standard Deviation") +
#   theme_bw() +
#   legend_none()
# 
# ppdSD2dens <- density(ppdSD2)
# ppdSD2dens <- cbind(ppdSD2dens$x, ppdSD2dens$y)
# ppdSD2densB <- ppdSD2dens[between(ppdSD2dens[,1], quantile(ppdSD2, 0.025), quantile(ppdSD2, 0.975)), ] 
# ppd_sd_plot2B <- ggplot() +
#   # geom_histogram(aes(x = ppdSD2,  after_stat(density)),
#   #                fill = "#bcdcdc", color = "#99c7c7") +
#   geom_ribbon(aes(x = ppdSD2densB[,1], ymin = 0, ymax = ppdSD2densB[,2], linetype = "A"),
#               fill = "#bcdcdc") +
#   geom_density(aes(x = ppdSD2), color = "#99c7c7", linewidth = .75) +
#   #geom_vline(aes(xintercept = quantile(ppdSD2, 0.975)), color = "#007C7C", linewidth = 2) +
#   geom_vline(aes(xintercept = DppdSD2, linetype = "B"), color = "#99c7c7", alpha = 0, linewidth = 0.75) +
#   geom_vline(aes(xintercept = DppdSD2, linetype = "C"), color = "#007C7C", linewidth = 1) +
#   geom_text(aes(label = paste0("p-value = ", round(mean(ppdSD2 > DppdSD2),3))),
#             x = 0.965*max(ppdSD2), y = max(ppdSD2dens[,2]), size = 3) +
#   scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
#   labs(title = "Standard Deviation",
#        x = NULL,
#        y = "Posterior Density") +
#   scale_linetype_manual(
#     name = element_blank(),
#     values = c(1,1,1),
#     breaks = c("C", "B", "A"),
#     labels = c("Observed", "Posterior", "95% Credible Interval")
#   ) +
#   guides(
#     linetype = guide_legend(
#       override.aes = list(
#         alpha = c(1,1,1),
#         shape = c(NA, "--", NA)
#       )
#     )
#   ) +
#   theme_bw() +
#   theme(
#     plot.title = element_text(size = rel(1)))
# ppd_sd_plot2B
# #ppd_legend2b <- ggpubr::get_legend(ppd_sd_plot2B)
# 
#### Combination Plot 
# ppd_lay2 <- rbind(c(NA,NA, rep(1,8),NA, NA),
#                   c(rep(2,4),rep(3,4),rep(4,4)),
#                   c(NA, NA, rep(5,4), rep(6,6))
# )
# 
# bayesplot_grid(
#   plots = list(
#     ppd_density_plot2,
#     ppd_q2.5_plot2B,
#     ppd_median_plot2B,
#     ppd_q97.5_plot2B,
#     ppd_mean_plot2B,
#     ppd_sd_plot2B),
#   grid_args = list(
#     layout_matrix = ppd_lay2
#   )
# )



## Model Comparison ----
### DIC ----
# dic2   <- dic.samples(model2, n.iter=2000, progress.bar="text")
# 
#### SAVE CHECKPOINT ----
# save(dic2, file = "Model 2 DIC.RData")
# 
### WAIC ----
# waic2   <- coda.samples(model2, 
#                         variable.names=c("like"), 
#                         n.iter=2000, progress.bar="text")
# 
#### SAVE CHECKPOINT ----
# save(waic2, file = "Model 2 WAIC.RData")
# 
#### Manually Compute WAIC and P ----
# like2   <- waic2[[1]]
# fbar2   <- colMeans(like2)
# P2      <- sum(apply(log(like2),2,var))
# WAIC2   <- -2*sum(log(fbar2))+2*P2
# 
### Output Results ----
# dic2
# WAIC2
# P2

## FINAL SAVE ----
save(
  # Data
  final_data3,
  Y2,
  Y2test,
  Y2train,
  X2,
  X2train,
  X2test,
  trainIndex,
  
  # MCMC Samples
  samples2,
  
  # Statisitcs
  quantiles2,
  stats2,
  predStats2B,
  
  # File Name
  file = filename2
)
