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

# Load Data ----
load("Model 1 Fit Data (B5000-I10000-T5).RData")
load("Model 1 Draws.RData")

# Old Plots ----
## Closed Support Densities ----
### Percent_Bleaching ----
ggplot() +
  geom_histogram(data = final_data3,
                 aes(x = Percent_Bleaching, after_stat(density)),
                 color = "black", fill = "blue",
                 bins = 100) +
  geom_density(data = final_data3,
               aes(x = Percent_Bleaching),
               color = "red",
               linewidth = 1) +
  theme_bw()

#### by Date_Year ----
ggplot() +
  geom_density(data = final_data3,
               aes(x = Percent_Bleaching, color = factor(Date_Year)),
               linewidth = 1) +
  theme_bw()

#### by Date_Month ----
ggplot() +
  geom_density(data = final_data3,
               aes(x = Percent_Bleaching, color = factor(Date_Month)),
               linewidth = 1) +
  theme_bw()

ggplot() +
  geom_density(data = final_data3,
               aes(x = log(Percent_Bleaching), color = factor(Date_Year)), linewidth = 1) +
  theme_bw()

## Open Support Densities ----
final_data3$Percent_Bleaching_Open <- 
  ifelse(final_data3$Percent_Bleaching == 0, 0.001, 
    ifelse(final_data3$Percent_Bleaching == 100, 99.999, final_data3$Percent_Bleaching))
  
### Percent_Bleaching ----
ggplot() +
  geom_histogram(data = final_data3,
                 aes(x = Percent_Bleaching_Open, after_stat(density)),
                 color = "black", fill = "blue",
                 bins = 100) +
  geom_density(data = final_data3,
               aes(x = Percent_Bleaching_Open),
               color = "red",
               linewidth = 1) +
  theme_bw()

#### by Date_Year ----
ggplot() +
  geom_density(data = final_data3,
               aes(x = Percent_Bleaching_Open, color = factor(Date_Year)),
               linewidth = 1) +
  theme_bw()

#### by Date_Month ----
ggplot() +
  geom_density(data = final_data3,
               aes(x = Percent_Bleaching_Open, color = factor(Date_Month)),
               linewidth = 1) +
  theme_bw()

## Log Support Densities ----
final_data3$Percent_Bleaching_log <- log(final_data3$Percent_Bleaching_Open)
  
### Percent_Bleaching ----
ggplot() +
  geom_histogram(data = final_data3,
                 aes(x = Percent_Bleaching_log, after_stat(density)),
                 color = "black", fill = "blue",
                 bins = 100) +
  geom_density(data = final_data3,
               aes(x = Percent_Bleaching_log),
               color = "red",
               linewidth = 1) +
  theme_bw()

#### by Date_Year ----
ggplot() +
  geom_density(data = final_data3,
               aes(x = Percent_Bleaching_log, color = factor(Date_Year)),
               linewidth = 1) +
  theme_bw()

#### by Date_Month ----
ggplot() +
  geom_density(data = final_data3,
               aes(x = Percent_Bleaching_log, color = factor(Date_Month)),
               linewidth = 1) +
  theme_bw()


plot(samples1)

###
# For all plots in one
par(mfrow = c(3,2))
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


# New Plots ----
## Load bayesian plotting library
library(bayesplot)

example <- example_mcmc_draws(chains = 3, params = 6)

## Trace plots only
color_scheme_set("mix-red-blue")
mcmc_trace(example, regex_pars = c("alpha", "beta")) +
  #scale_color_manual(values = c("red","blue")) +
  theme_bw()

## Density plots only
mcmc_dens_overlay(draws2, regex_pars = c("alpha", "beta")) +
  #scale_color_manual(values = c("red","blue")) +
  theme_bw()

## Combo plot
mcmc_combo(draws2, regex_pars = c("alpha", "beta"), 
           combo = c("dens_overlay", "trace")) +
  scale_color_manual(values = c("red","blue")) +
  theme_bw()

## Convert MCMC List to Draws Array
draws1 <- as_draws_list(samples1)
draws2 <- posterior::as_draws_array(samples1)

## Trace plots only
color_scheme_set("mix-red-blue")
mcmc_trace(draws2, regex_pars = c("alpha", "beta")) +
  #scale_color_manual(values = c("red","blue")) +
  theme_bw()

## Density plots only
mcmc_dens_overlay(draws2, regex_pars = c("alpha", "beta")) +
  #scale_color_manual(values = c("red","blue")) +
  theme_bw()

## Combo plot
mcmc_combo(draws2, regex_pars = c("alpha", "beta"), 
           combo = c("dens_overlay", "trace")) +
  scale_color_manual(values = c("red","blue")) +
  theme_bw()


