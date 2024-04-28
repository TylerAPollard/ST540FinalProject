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


byYearfun <- function(x){
  mean(x)
}
byYear_EDA <- ddply(final_data3, .(Date_Year), summarize,
      ClimSST = byYearfun(ClimSST),
      MedClimSST = median(ClimSST),
      SSTA = byYearfun(SSTA),
      SSTA_DHW = byYearfun(SSTA_DHW),
      TSA = byYearfun(TSA),
      TSA_DHW = byYearfun(TSA_DHW),
      Percent_Bleaching = byYearfun(Percent_Bleaching))

temp_data <- final_data3 |> select(
  ClimSST,
  SSTA,
  SSTA_DHW,
  TSA,
  TSA_DHW,
  Percent_Bleaching
)
ggpairs(temp_data)

lm_full <- lm(data = temp_data,
          Percent_Bleaching ~ .)
summary(lm_full)
anova(lm_full)
step(lm_full)

temp_data2U <- final_data3 |>
  filter(TSA_DHW >= 20)
temp_data2L <- final_data3 |>
  filter(TSA_DHW < 20)


ggplot() +
  # geom_histogram(data = final_tyler,
  #                aes(x = Percent_Bleaching, after_stat(density), fill = factor(Date_Year)), 
  #                color = "black", bins = 100) +
  geom_density(data = final_data3,
               aes(x = Percent_Bleaching, color = factor(Date_Year)), linewidth = 1) +
  theme_bw()

ggplot() +
  # geom_histogram(data = final_tyler,
  #                aes(x = Percent_Bleaching, after_stat(density), fill = factor(Date_Year)), 
  #                color = "black", bins = 100) +
  geom_density(data = final_data3,
               aes(x = log(Percent_Bleaching), color = factor(Date_Year)), linewidth = 1) +
  theme_bw()

unique(final_data3$City_Town_Name)
ddply(final_data3, .(Date_Year), summarize,
      Obs = length(Percent_Bleaching),
      Mean = mean(Percent_Bleaching))

world_coordinates <- map_data("county") 
ggplot() + 
  # geom_map() function takes world coordinates  
  # as input to plot world map 
  geom_map( 
    data = world_coordinates, map = world_coordinates, 
    aes(x = long, y = lat, map_id = region) 
  ) + 
  geom_point(
    data = final_data3,
    aes(x = Longitude_Degrees, y = Latitude_Degrees, 
        color = Percent_Bleaching)
  ) +
  xlim(c(-85,-77.5)) +
  ylim(c(23,32.5)) +
  scale_color_continuous(low = "green", high = "red") +
  facet_wrap(vars(City_Town_Name)) +
  theme_bw()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Model 3: Exponential regression model (Tyler) ----
## Modeled with Uninformative Gaussian Priors
Y3 <- final_tyler$Percent_Bleaching
mean(Y3)
sd(Y3)
zeroindex <- Y3 == 0
Y3zero <- Y3[zeroindex]
Y3nonzero <- Y3[!zeroindex]
# plot(density(Y3), xlim = c(0,100))
# plot(density(Y3zero), xlim = c(0,100))
# plot(density(Y3nonzero), xlim = c(0,100))

# Change 0's to 0.00001 and 1's to 0.99999 to avoid infinite density
# Log link
Y3B <- ifelse(Y3 == 0, 0.001, 
              ifelse(Y3 == 100, 99.999, Y3))
mean(Y3B)
sd(Y3B)
Y3B <- log(Y3B)
# plot(density(Y3B))

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

# Sqrt link
Y3D <- sqrt(Y3)
plot(density(Y3D))

Z3 <- ifelse(Y3 == 0, 0, 1)

X3 <- final_tyler |> 
  select(-c(
    Date,
    City_Town_Name,
    Latitude_Degrees,
    Longitude_Degrees,
    # Date_Year,
    # Exposure,
    # Turbidity,
    ClimSST,
    # SSTA,
    # TSA,
    # Windspeed,
    Percent_Bleaching
  ))
X3$Exposure <- ifelse(X3$Exposure == "Sheltered", 0, 1)
X3 <- scale(X3)

cor(X3, Y3)
ggpairs(X3)

n3 <- length(Y3)
p3 <- ncol(X3)

burn     <- 500
n.iter   <- 1000
thin     <- 5

model_string3A <- textConnection("model{
    # Likelihood
    for(i in 1:n){
      Y[i] ~ dnorm(mu[i], taue)
      mu[i] <- alphaY + inprod(X[i,], betaY[])
      
      Z[i] ~ dbern(q[i])
      logit(q[i]) <- alphaZ + inprod(X[i,], betaZ[])
    } 
    
    # Priors  
    for(j in 1:p){ 
      betaY[j] ~ dnorm(0, 0.01)
      betaZ[j] ~ dnorm(0, 0.01)
    } 
    alphaY ~ dnorm(0, 0.01)
    alphaZ ~ dnorm(0, 0.01)
    
    taue ~ dgamma(0.1,0.1)
    sigma <- 1/sqrt(taue)
    
    # # Posterior Predicitve Checks
    # for(i in 1:n){
    #   Zppc[i] ~ dbern(q[i])
    # }
    # D3[1] <- mean(Zppc[])
    # D3[2] <- sd(Zppc[])
}")

data3   <- list(Z=Z3,X=X3,n=n3,p=p3)
params3 <- c("alphaY", "betaY", 
             "alphaZ", "betaZ",
             "sigma"
             )

inits <- list(
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 52),
  list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 8)
)

## Run this all at once to get how long it took
# Start here
tic()
model3A <- jags.model(model_string3A, data=data3, inits = inits,
                     n.chains=2, quiet=FALSE)
#error persists. i've tried increasing burn in. next objective, is to review priors or perhaps try using a different sampler
update(model3A, burn, progress.bar = "none")
samples3A <- coda.samples(model3A, variable.names=params3, n.iter=n.iter, n.thin=thin, 
                         progress.bar="none")
toc()
summary3A <- summary(samples3A)
summary3A

save(
  final_data3,
  final_tyler,
  model3A,
  samples3A,
  summary3A,
  X3,
  Z3,
  n3,
  p3,
  file = "Model 3A Data.RData"
)

plot(samples3A)

# This reduces chance of crashing
dev.off()

stats3A <- summary3A$statistics[c(2397:2408),]
rownames(stats3A) <- c("Intercept", colnames(X3))
stats3A

quantiles3A <- summary3A$quantiles[c(2397:2408),]
rownames(quantiles3A) <- c("Intercept", colnames(X3))
quantiles3A

model_string3 <- textConnection("model{
    # Likelihood
    for(i in 1:n){
      Y[i] ~ dnorm(mu[i], taue) # Percent_Bleaching
      mu[i] <- alphaY + inprod(X[i,], betaY[])
      
      Z[i] ~ dbern(q[i]) # Whether zero or not
      q[i] <- alphaZ + inprod(X[i,], betaZ[])
    } 
    
    # Priors  
    for(j in 1:p){
      betaY[j] ~ dnorm(0, 0.01)
      betaZ[j] ~ dnorm(0, 0.01)
    }
    alphaY ~ dnorm(0, 0.01)
    alphaZ ~ dnorm(0, 0.01)
    
    taue ~ dgamma(0.1, 0.1)
    sigma 1/sqrt(taue)
    
}")

data3   <- list(Y = Y3,
                X = X3,
                Z = Z3,
                n = n3,
                p = p3
                )
params3 <- c("alphaY", "betaY",
             "alphaZ", "betaZ",
             "sigma"
             )

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
final_data3$Percent_Bleaching_log <- log(final_data3$Percent_Bleaching_Open)

## Create training/test index vector ----
## for CV and testing model performance and fit
set.seed(52)
trainIndex <- createDataPartition(final_data3$Percent_Bleaching,
                                  p = 0.75,
                                  list = FALSE)
trainIndex <- as.vector(trainIndex)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Model 3: Zero-Inflated Regression ----
## Programmed by Tyler
## Modeled Heirarchically

## Load Data ----
### Closed Support Data ----
## Y3 is original data with closed support [0,100]
Y3 <- final_data3$Percent_Bleaching
Z3 <- ifelse(Y3 == 0, 0, 1)

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
    #ClimSST,
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
      log(mu[i]) <- alphaY + inprod(X[i,], betaY[])
      
      Z[i] ~ dbern(q[i])
      logit(q[i]) <- alphaZ + inprod(X[i,], betaZ[])
    } 
    
    # Priors  
    for(j in 1:p){ 
      betaY[j] ~ dnorm(0, 0.01)
      betaZ[j] ~ dnorm(0, 0.01)
    } 
    alphaY ~ dnorm(0, 0.01)
    alphaZ ~ dnorm(0, 0.01)
    
    taue ~ dgamma(0.1,0.1)
    sigma <- 1/sqrt(taue)
    
    # # Posterior Predicitve Checks
    # for(i in 1:n){
    #   Zppc[i] ~ dbern(q[i])
    # }
    # D3[1] <- mean(Zppc[])
    # D3[2] <- sd(Zppc[])
}")

### Compile Model Inputs ----
data3   <- list(Y = Y3train,
                Z = Z3train,
                X = X3train,
                n = n3train,
                p = p3)
params3 <- c("alphaY", "betaY", 
             "alphaZ", "betaZ",
             "sigma")

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
paramNames <- rep(c("Intercept", colnames(X3)), each = 2)
paramNames2 <- paste0(c("Y_", "Z_"), paramNames)
paramNames3 <- c(paramNames2, "sigma")

stats3 <- summary3$statistics
rownames(stats3) <- paramNames3
stats3

#### Parameter 95% CIs ----
# CI = Credible Interval
quantiles3 <- summary3$quantiles
rownames(quantiles3) <- paramNames3
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
D0 <- c(
  min(Y3),
  max(Y3),
  max(Y3) - min(Y3),
  mean(Y3),
  sd(Y3)
)
names(D0) <- DPrintnames

# Chain 3 PPCs
D3A <- samples3[[3]][,3:5]
colnames(D3A) <- DPrintnames

# Chain 2 PPCs
D3B <- samples3[[2]][,3:5]
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

for(j in 3:5){
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



