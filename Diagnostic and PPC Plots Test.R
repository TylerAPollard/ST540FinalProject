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
Y <- final_data3 |> pull(Percent_Bleaching)
Y <- Y/100
Nonzer0 <- final_data3 |> filter(Percent_Bleaching != 0) |> pull(Percent_Bleaching)
Nonzer0 <- Nonzer0/100
mean(Nonzer0/100)
sd(Nonzer0/100)
var(Nonzer0/100)
Nonzer0B <- log(Nonzer0)
Nonzer0C <- log(Nonzer0 + var(Nonzer0)/2)

curve(dnorm(x, mean(Y), sd(Y)))
curve(dnorm(x, mean(Nonzer0), sd(Nonzer0)))
curve(dnorm(x, mean(Nonzer0B), sd(Nonzer0B)))
curve(dnorm(x, mean(Nonzer0C), sd(Nonzer0C)))

curve(dnorm(x, mean(log(Y)), sd(log(Y))))
curve(dnorm(x, mean(log(Nonzer0)), sd(log(Nonzer0))))
curve(dnorm(x, mean(log(Nonzer0B)), sd(log(Nonzer0B))))
curve(dnorm(x, mean(Nonzer0C), sd(Nonzer0C)))

ggplot() +
  geom_histogram(#data = final_data3 |> filter(Percent_Bleaching != 0),
    aes(x = Y, after_stat(density)),
    color = "#99c7c7", fill = "#bcdcdc",
    bins = 100) +
  geom_density(#data = final_data3,
    aes(x = Y),
    color = "#007C7C", 
    linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  #scale_x_continuous(limits = c(-0.5,1)) +
  labs(title = "Density Plot of Percent Bleaching from 2,394 Coral Reef Samples",
       subtitle = "Data was Collected by the Florida Reef Resilience Program from 2006 to 2016",
       x = "Percent Bleaching",
       y = "Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

ggplot() +
  geom_histogram(#data = final_data3 |> filter(Percent_Bleaching != 0),
    aes(x = Nonzer0, after_stat(density)),
    color = "#99c7c7", fill = "#bcdcdc",
    bins = 100) +
  geom_density(#data = final_data3,
    aes(x = Nonzer0),
    color = "#007C7C", 
    linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  #scale_x_continuous(limits = c(-0.5,1)) +
  labs(title = "Density Plot of Percent Bleaching from 2,394 Coral Reef Samples",
       subtitle = "Data was Collected by the Florida Reef Resilience Program from 2006 to 2016",
       x = "Percent Bleaching",
       y = "Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

ggplot() +
  geom_histogram(#data = final_data3 |> filter(Percent_Bleaching != 0),
    aes(x = Nonzer0B, after_stat(density)),
    color = "#99c7c7", fill = "#bcdcdc",
    bins = 100) +
  geom_density(#data = final_data3,
    aes(x = Nonzer0B),
    color = "#007C7C", 
    linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  #scale_x_continuous(limits = c(-0.5,1)) +
  labs(title = "Density Plot of Percent Bleaching from 2,394 Coral Reef Samples",
       subtitle = "Data was Collected by the Florida Reef Resilience Program from 2006 to 2016",
       x = "Percent Bleaching",
       y = "Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

ggplot() +
  geom_histogram(#data = final_data3 |> filter(Percent_Bleaching != 0),
                 aes(x = Nonzer0C, after_stat(density)),
                 color = "#99c7c7", fill = "#bcdcdc",
                 bins = 100) +
  geom_density(#data = final_data3,
               aes(x = Nonzer0C),
               color = "#007C7C", 
               linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  #scale_x_continuous(limits = c(-0.5,1)) +
  labs(title = "Density Plot of Percent Bleaching from 2,394 Coral Reef Samples",
       subtitle = "Data was Collected by the Florida Reef Resilience Program from 2006 to 2016",
       x = "Percent Bleaching",
       y = "Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

ggplot() +
  geom_histogram(#data = final_data3 |> filter(Percent_Bleaching != 0),
    aes(x = log(Nonzer0B + var(Nonzer0B)/2), after_stat(density)),
    color = "#99c7c7", fill = "#bcdcdc",
    bins = 100) +
  geom_density(#data = final_data3,
    aes(x = log(Nonzer0B + var(Nonzer0B)/2)),
    color = "#007C7C", 
    linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  #scale_x_continuous(limits = c(-0.5,1)) +
  labs(title = "Density Plot of Percent Bleaching from 2,394 Coral Reef Samples",
       subtitle = "Data was Collected by the Florida Reef Resilience Program from 2006 to 2016",
       x = "Percent Bleaching",
       y = "Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

##### Back

ggplot() +
  geom_histogram(#data = final_data3 |> filter(Percent_Bleaching != 0),
    aes(x = exp(Nonzer0C), after_stat(density)),
    color = "#99c7c7", fill = "#bcdcdc",
    bins = 100) +
  geom_density(#data = final_data3,
    aes(x = exp(Nonzer0C)),
    color = "#007C7C", 
    linewidth = 1) +
  geom_density(#data = final_data3,
    aes(x = Nonzer0),
    color = "red", 
    linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  #scale_x_continuous(limits = c(-0.5,1)) +
  labs(title = "Density Plot of Percent Bleaching from 2,394 Coral Reef Samples",
       subtitle = "Data was Collected by the Florida Reef Resilience Program from 2006 to 2016",
       x = "Percent Bleaching",
       y = "Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

ggplot() +
  geom_histogram(#data = final_data3 |> filter(Percent_Bleaching != 0),
    aes(x = exp(Nonzer0), after_stat(density)),
    color = "#99c7c7", fill = "#bcdcdc",
    bins = 100) +
  geom_density(#data = final_data3,
    aes(x = Nonzer0),
    color = "#007C7C", 
    linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  #scale_x_continuous(limits = c(-0.5,1)) +
  labs(title = "Density Plot of Percent Bleaching from 2,394 Coral Reef Samples",
       subtitle = "Data was Collected by the Florida Reef Resilience Program from 2006 to 2016",
       x = "Percent Bleaching",
       y = "Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

ggplot() +
  geom_histogram(#data = final_data3 |> filter(Percent_Bleaching != 0),
    aes(x = Nonzer0B, after_stat(density)),
    color = "#99c7c7", fill = "#bcdcdc",
    bins = 100) +
  geom_density(#data = final_data3,
    aes(x = Nonzer0B),
    color = "#007C7C", 
    linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  #scale_x_continuous(limits = c(-0.5,1)) +
  labs(title = "Density Plot of Percent Bleaching from 2,394 Coral Reef Samples",
       subtitle = "Data was Collected by the Florida Reef Resilience Program from 2006 to 2016",
       x = "Percent Bleaching",
       y = "Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

ggplot() +
  geom_histogram(#data = final_data3 |> filter(Percent_Bleaching != 0),
    aes(x = Nonzer0C, after_stat(density)),
    color = "#99c7c7", fill = "#bcdcdc",
    bins = 100) +
  geom_density(#data = final_data3,
    aes(x = Nonzer0C),
    color = "#007C7C", 
    linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  #scale_x_continuous(limits = c(-0.5,1)) +
  labs(title = "Density Plot of Percent Bleaching from 2,394 Coral Reef Samples",
       subtitle = "Data was Collected by the Florida Reef Resilience Program from 2006 to 2016",
       x = "Percent Bleaching",
       y = "Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

ggplot() +
  geom_histogram(#data = final_data3 |> filter(Percent_Bleaching != 0),
    aes(x = log(Nonzer0B + var(Nonzer0B)/2), after_stat(density)),
    color = "#99c7c7", fill = "#bcdcdc",
    bins = 100) +
  geom_density(#data = final_data3,
    aes(x = log(Nonzer0B + var(Nonzer0B)/2)),
    color = "#007C7C", 
    linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  #scale_x_continuous(limits = c(-0.5,1)) +
  labs(title = "Density Plot of Percent Bleaching from 2,394 Coral Reef Samples",
       subtitle = "Data was Collected by the Florida Reef Resilience Program from 2006 to 2016",
       x = "Percent Bleaching",
       y = "Density") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )









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
                 aes(x = Percent_Bleaching_Log, after_stat(density)),
                 color = "black", fill = "blue",
                 bins = 100) +
  geom_density(data = final_data3,
               aes(x = Percent_Bleaching_Log),
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











